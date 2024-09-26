#include "kelement_mlt.h"
#include "lodepng/lodepng.h"
#include "vec3.h"
#include <climits>
#include <iomanip>
#include <iostream>
#include <thread>
#include <vector>

struct Ray {
    Vec3 o, d;
    Ray(const Vec3 &o_, const Vec3 &d_) : o(o_), d(d_) {}
};

struct PathSample {
    int x, y;

    Vec3 F;
    double weight;

    PathSample() : x(-1), y(-1), F(Vec3(0.0)), weight(NAN) {}

    PathSample(const int x_, const int y_, const Vec3 &F_, const double weight_)
        : x(x_), y(y_), F(F_), weight(weight_) {}
};

enum class ReflectionType { diffuse, specular, Refractive }; // material types, used in radiance()

struct Sphere {
    double radius; // radius

    Vec3 position;
    Vec3 emission;
    Vec3 color;

    ReflectionType reflection_type;

    Sphere(double rad_, const Vec3 &p_, const Vec3 &e_, const Vec3 &c_,
           ReflectionType _reflection_type)
        : radius(rad_), position(p_), emission(e_), color(c_), reflection_type(_reflection_type) {}

    [[nodiscard]] double intersect(const Ray &r) const { // returns distance, 0 if nohit
        Vec3 op = position - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double epsilon = 1e-4;
        double b = op.dot(r.d);
        double determinant = b * b - op.dot(op) + radius * radius;
        if (determinant < 0) {
            return 0;
        }

        determinant = sqrt(determinant);

        if (double t = b - determinant; t > epsilon) {
            return t;
        }

        if (double t = b + determinant; t > epsilon) {
            return t;
        }

        return 0;
    }
};

std::vector<Sphere> spheres = {
    // Scene: radius, position, emission, color, material
    Sphere(1e5, Vec3(1e5 + 1, 40.8, 81.6), Vec3(0, 0, 0), Vec3(.75, .25, .25),
           ReflectionType::diffuse), // Left
    Sphere(1e5, Vec3(-1e5 + 99, 40.8, 81.6), Vec3(0, 0, 0), Vec3(.25, .25, .75),
           ReflectionType::diffuse), // Rght
    Sphere(1e5, Vec3(50, 40.8, 1e5), Vec3(0, 0, 0), Vec3(.75, .75, .75),
           ReflectionType::diffuse), // Back
    Sphere(1e5, Vec3(50, 40.8, -1e5 + 170), Vec3(0, 0, 0), Vec3(0, 0, 0),
           ReflectionType::diffuse), // Frnt
    Sphere(1e5, Vec3(50, 1e5, 81.6), Vec3(0, 0, 0), Vec3(.75, .75, .75),
           ReflectionType::diffuse), // Botm
    Sphere(1e5, Vec3(50, -1e5 + 81.6, 81.6), Vec3(0, 0, 0), Vec3(.75, .75, .75),
           ReflectionType::diffuse), // Top
    Sphere(16.5, Vec3(27, 16.5, 47), Vec3(0, 0, 0), Vec3(1, 1, 1) * .999,
           ReflectionType::specular), // Mirr
    Sphere(16.5, Vec3(73, 16.5, 78), Vec3(0, 0, 0), Vec3(1, 1, 1) * .999,
           ReflectionType::Refractive), // Glas
    Sphere(600, Vec3(50, 681.6 - .27, 81.6), Vec3(12, 12, 12), Vec3(0, 0, 0),
           ReflectionType::diffuse) // Lite
};

inline double clamp(double x, double low, double high) {
    return x < low ? low : x > high ? high : x;
}

inline int toInt(double x) {
    // with gamma correction
    return int(pow(clamp(x, 0, 1), 1 / 2.2) * 255 + .5);
}

inline int intersect(const Ray &r, double &t) {
    int id = -1;
    t = std::numeric_limits<double>::infinity();

    for (int i = 0; i < spheres.size(); ++i) {
        double d = spheres[i].intersect(r);
        if (d > 0 && d < t) {
            t = d;
            id = i;
        }
    }

    return id;
}

Vec3 trace(const Ray &camera_ray, KelemenMLT &mlt) {
    Vec3 radiance(0.0, 0.0, 0.0);
    Vec3 throughput(1.0, 1.0, 1.0);

    auto ray = camera_ray;

    for (int depth = 0; depth < 10; ++depth) {
        double t; // distance to intersection
        int hit_sphere_id = intersect(ray, t);
        if (hit_sphere_id < 0) {
            break;
        }

        const Sphere &obj = spheres[hit_sphere_id]; // the hit object
        Vec3 hit_point = ray.o + ray.d * t;
        Vec3 surface_normal = (hit_point - obj.position).norm(); // always face out
        Vec3 normal = surface_normal.dot(ray.d) < 0 ? surface_normal : surface_normal * -1;

        radiance += throughput * obj.emission;
        throughput *= obj.color;

        if (depth > 4) {
            // russian roulette
            double probability_russian_roulette = clamp(throughput.max_component_val(), 0.1, 0.95);

            if (mlt.next_sample() >= probability_russian_roulette) {
                // terminated
                break;
            }
            // survive and enhanced
            throughput *= (1.0 / probability_russian_roulette);
        }

        if (obj.reflection_type == ReflectionType::diffuse) { // Ideal DIFFUSE reflection
            double r1 = 2 * M_PI * mlt.next_sample();
            double r2 = mlt.next_sample();
            double r2s = sqrt(r2);
            Vec3 w = normal;
            Vec3 u = (fabs(w.x) > 0.1 ? Vec3(0, 1, 0) : Vec3(1, 0, 0)).cross(w).norm();
            Vec3 v = w.cross(u);
            Vec3 d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

            ray = Ray(hit_point, d);
            continue;
        }

        if (obj.reflection_type == ReflectionType::specular) { // Ideal SPECULAR reflection
            ray = Ray(hit_point, ray.d - surface_normal * 2 * surface_normal.dot(ray.d));
            continue;
        }

        Ray spawn_ray_reflect(hit_point,
                              ray.d - surface_normal * 2 *
                                          surface_normal.dot(ray.d)); // Ideal dielectric REFRACTION

        bool into = surface_normal.dot(normal) > 0; // Ray from outside going in?
        double nc = 1;
        double nt = 1.5;
        double nnt = into ? nc / nt : nt / nc;
        double ddn = ray.d.dot(normal);
        double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

        if (cos2t < 0) { // Total internal reflection
            ray = spawn_ray_reflect;
            continue;
        }

        Vec3 t_dir =
            (ray.d * nnt - surface_normal * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
        double a = nt - nc;
        double b = nt + nc;
        double R0 = a * a / (b * b);
        double c = 1 - (into ? -ddn : t_dir.dot(surface_normal));

        double Re = R0 + (1 - R0) * c * c * c * c * c;
        double Tr = 1 - Re;
        double probability_reflect = 0.25 + 0.5 * Re;

        double RP = Re / probability_reflect;
        double TP = Tr / (1 - probability_reflect);

        // refract or reflect
        if (mlt.next_sample() < probability_reflect) {
            // reflect
            ray = spawn_ray_reflect;
            throughput *= RP;
            continue;
        }

        // refract
        ray = Ray(hit_point, t_dir); // Ideal dielectric REFRACTION
        throughput *= TP;
    }

    return radiance;
}

PathSample generate_new_path(KelemenMLT &mlt, const int width, const int height) {
    const double weight = 4.0 * width * height;

    int x = mlt.next_sample() * width;
    x = x % width;

    int y = mlt.next_sample() * height;
    y = y % height;

    Ray cam(Vec3(50, 52, 295.6), Vec3(0, -0.042612, -1).norm()); // cam pos, dir
    Vec3 cx = Vec3(width * 0.5135 / height, 0, 0);
    Vec3 cy = cx.cross(cam.d).norm() * 0.5135;

    const double r1 = 2.0 * mlt.next_sample();
    const double r2 = 2.0 * mlt.next_sample();

    double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
    double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

    Vec3 d = cx * ((dx + x) / width - 0.5) + cy * ((dy + y - 20) / height - 0.5) + cam.d;

    auto c = trace(Ray(cam.o + d * 140, d.norm()), mlt);

    return PathSample(x, y, c, weight);
}

void render_with_selected_path(std::vector<Vec3> &image, KelemenMLT mlt, const int seed,
                               const PathSample &init_path, const double b,
                               const long long total_mutations, int width, int height) {
    RNG rng(seed);

    const double p_large = 0.5;
    int num_accept_mutation = 0;
    int num_reject_mutation = 0;

    auto old_path = init_path;

    for (long long i = 0; i < total_mutations; i++) {
        // この辺も全部論文と同じ（Next()）
        mlt.init_sample_idx();
        mlt.large_step = rng() < p_large;
        PathSample new_path = generate_new_path(mlt, width, height);

        double accept_prob = std::min(1.0, new_path.F.luminance() / old_path.F.luminance());
        const double new_path_weight =
            (accept_prob + mlt.large_step) / (new_path.F.luminance() / b + p_large);
        const double old_path_weight = (1.0 - accept_prob) / (old_path.F.luminance() / b + p_large);

        image[new_path.y * width + new_path.x] += new_path_weight * new_path.F;
        image[old_path.y * width + old_path.x] += old_path_weight * old_path.F;

        if (rng() < accept_prob) {
            num_accept_mutation++;
            old_path = new_path;
            if (mlt.large_step) {
                mlt.large_step_time = mlt.global_time;
            }
            mlt.global_time++;
            mlt.clear_backup_samples();
        } else {
            num_reject_mutation++;
            mlt.recover_samples();
        }
    }
}

void render_mlt(std::vector<Vec3> &pixels, int mutation_per_pixel, int width, int height) {
    auto start = std::chrono::system_clock::now();

    RNG rng(INT_MAX / 4);
    KelemenMLT mlt(INT_MAX);

    const int num_seed_paths = width * height;
    std::vector<PathSample> seed_paths(num_seed_paths);
    double sumI = 0.0;
    mlt.large_step = 1;
    for (int i = 0; i < num_seed_paths; i++) {
        mlt.init_sample_idx();
        PathSample sample = generate_new_path(mlt, width, height);

        mlt.global_time++;
        mlt.clear_backup_samples();

        sumI += sample.F.luminance();
        seed_paths[i] = sample;
    }

    // First pass is done, using importance sampling based on luminance values.
    // TODO: why is this selection important?
    int selected_path_id = -1;
    const double threshold = rng() * sumI;
    double accumulated_importance = 0.0;
    for (int i = 0; i < num_seed_paths; i++) {
        accumulated_importance += seed_paths[i].F.luminance();
        if (accumulated_importance >= threshold) {
            selected_path_id = i;
            break;
        }
    }

    const std::chrono::duration<double> duration_initial{std::chrono::system_clock::now() - start};
    std::cout << "building initial paths took " << std::fixed << std::setprecision(3)
              << duration_initial.count() << " seconds.\n"
              << std::flush;

    const double b = sumI / num_seed_paths;

    auto start_render = std::chrono::system_clock::now();

    auto seed = 42;

    const long long total_mutations = (long long)(width * height) * mutation_per_pixel;
    render_with_selected_path(pixels, mlt, seed, seed_paths[selected_path_id], b, total_mutations,
                              width, height);

    for (int i = 0; i < width * height; i++) {
        pixels[i] = pixels[i] / mutation_per_pixel;
    }
}

int main() {
    const int width = 1280;
    const int height = 960;

    for (auto mutation_per_pixel : {16, 32, 64, 128, 256, 512}) {
        auto pixels = std::vector<Vec3>(width * height);

        auto start = std::chrono::system_clock::now();

        render_mlt(pixels, mutation_per_pixel, width, height);

        const std::chrono::duration<double> duration{std::chrono::system_clock::now() - start};
        std::cout << "rendering (" << mutation_per_pixel << " spp) took " << std::fixed
                  << std::setprecision(3) << duration.count() << " seconds.\n"
                  << std::flush;

        std::vector<unsigned char> png_pixels(width * height * 4);

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int png_idx = y * width + x;
                int idx = (height - 1 - y) * width + x;

                png_pixels[4 * png_idx + 0] = toInt(pixels[idx].x);
                png_pixels[4 * png_idx + 1] = toInt(pixels[idx].y);
                png_pixels[4 * png_idx + 2] = toInt(pixels[idx].z);
                png_pixels[4 * png_idx + 3] = 255;
            }
        }

        std::string file_name = "smallpt_" + std::to_string(mutation_per_pixel) + ".png";

        if (unsigned error = lodepng::encode(file_name, png_pixels, width, height); error) {
            std::cerr << "lodepng::encoder error " << error << ": " << lodepng_error_text(error)
                      << std::endl;
            throw std::runtime_error("lodepng::encode() fail");
        }

        std::cout << "image saved to `" << file_name << "`\n\n" << std::flush;
    }
}
