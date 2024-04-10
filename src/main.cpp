#include <iostream>
#include <mutex>
#include <random>
#include <stack>
#include <thread>
#include <vector>

#include "lodepng/lodepng.h"

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}

    Vec3(double _x, double _y, double _z) {
        x = _x;
        y = _y;
        z = _z;
    }

    Vec3 operator+(const Vec3 &b) const { return Vec3(x + b.x, y + b.y, z + b.z); }

    void operator+=(const Vec3 &b) {
        x += b.x;
        y += b.y;
        z += b.z;
    }

    Vec3 operator-(const Vec3 &b) const { return Vec3(x - b.x, y - b.y, z - b.z); }

    Vec3 operator*(double b) const { return Vec3(x * b, y * b, z * b); }

    void operator*=(double b) {
        x *= b;
        y *= b;
        z *= b;
    }

    Vec3 operator*(const Vec3 &b) const { return Vec3(x * b.x, y * b.y, z * b.z); }

    void operator*=(const Vec3 &b) {
        x *= b.x;
        y *= b.y;
        z *= b.z;
    }

    Vec3 norm() const { return *this * (1 / sqrt(x * x + y * y + z * z)); }

    double dot(const Vec3 &b) const { return x * b.x + y * b.y + z * b.z; }

    Vec3 cross(const Vec3 &b) const {
        return Vec3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }

    double max_component_val() const { return x > y && x > z ? x : y > z ? y : z; }
};

struct Ray {
    Vec3 o, d;
    Ray(const Vec3 &o_, const Vec3 &d_) : o(o_), d(d_) {}
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

struct Sampler {
    Sampler() {
        uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::seed_seq ss{uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32)};
        rng.seed(ss);
        distribution_0_1 = std::uniform_real_distribution<double>(0, 1);
    }
    void seed(int seed_val) { rng.seed(seed_val); }

    double generate() { return distribution_0_1(rng); }

  private:
    std::mt19937_64 rng;
    std::uniform_real_distribution<double> distribution_0_1;
};

inline double clamp(double x, double low, double high) {
    return x < low ? low : x > high ? high : x;
}

inline int toInt(double x) { return int(pow(clamp(x, 0, 1), 1 / 2.2) * 255 + .5); }

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

Vec3 trace(const Ray &camera_ray, Sampler &sampler) {
    Vec3 radiance(0.0, 0.0, 0.0);
    Vec3 throughput(1.0, 1.0, 1.0);

    auto ray = camera_ray;

    for (int depth = 0;; ++depth) {

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

            if (sampler.generate() < probability_russian_roulette) {
                // survive and enhanced
                throughput *= (1.0 / probability_russian_roulette);
            } else {
                // terminated
                break;
            }
        }

        if (obj.reflection_type == ReflectionType::diffuse) { // Ideal DIFFUSE reflection
            double r1 = 2 * M_PI * sampler.generate();
            double r2 = sampler.generate();
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
        if (sampler.generate() < probability_reflect) {
            // reflect
            ray = spawn_ray_reflect;
            throughput *= RP;
            continue;
        }

        // refract
        ray = Ray(hit_point, t_dir); // Ideal dielectric REFRACTION
        throughput *= TP;
        continue;
    }

    return radiance;
}

void render(std::vector<Vec3> &pixels, std::stack<int> &job_list, std::mutex &mtx, int samples,
            int width, int height) {
    Ray cam(Vec3(50, 52, 295.6), Vec3(0, -0.042612, -1).norm()); // cam pos, dir
    Vec3 cx = Vec3(width * 0.5135 / height, 0, 0);
    Vec3 cy = cx.cross(cam.d).norm() * 0.5135;

    Sampler sampler;

    while (true) {
        mtx.lock();
        if (job_list.empty()) {
            mtx.unlock();
            return;
        }
        auto y = job_list.top();
        job_list.pop();
        mtx.unlock();

        for (int x = 0; x < width; x++) {
            int pixel_index = (height - y - 1) * width + x;
            sampler.seed(pixel_index);
            auto pixel_val = Vec3(0.0, 0.0, 0.0);

            for (int s = 0; s < samples; s++) {
                double r1 = 2 * sampler.generate();
                double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);

                double r2 = 2 * sampler.generate();
                double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

                Vec3 d = cx * (((0.5 + dx) / 2 + x) / width - 0.5) +
                         cy * (((0.5 + dy) / 2 + y) / height - 0.5) + cam.d;

                pixel_val = pixel_val + trace(Ray(cam.o + d * 140, d.norm()), sampler);
            } // Camera rays are pushed ^^^^^ forward to start in interior

            pixel_val = pixel_val * (1.0 / samples);

            pixels[pixel_index] =
                Vec3(clamp(pixel_val.x, 0, 1), clamp(pixel_val.y, 0, 1), clamp(pixel_val.z, 0, 1));
        }
    }
}

int main() {
    int width = 1024;
    int height = 768;

    int samples = 100;
    auto pixels = std::vector<Vec3>(width * height);

    std::stack<int> job_list;
    for (int y = 0; y < height; y++) {
        job_list.push(y);
    }

    std::mutex mtx;
    std::vector<std::thread> threads;
    for (int idx = 0; idx < std::thread::hardware_concurrency(); ++idx) {
        threads.push_back(std::thread(render, std::ref(pixels), std::ref(job_list), std::ref(mtx),
                                      samples, width, height));
    }

    for (auto &t : threads) {
        t.join();
    }

    std::string file_name = "small_pt_" + std::to_string(samples) + ".png";

    std::vector<unsigned char> png_pixels(width * height * 4);

    for (int i = 0; i < width * height; i++) {
        png_pixels[4 * i + 0] = toInt(pixels[i].x);
        png_pixels[4 * i + 1] = toInt(pixels[i].y);
        png_pixels[4 * i + 2] = toInt(pixels[i].z);
        png_pixels[4 * i + 3] = 255;
    }

    // Encode the image
    // if there's an error, display it
    if (unsigned error = lodepng::encode(file_name, png_pixels, width, height); error) {
        std::cerr << "lodepng::encoder error " << error << ": " << lodepng_error_text(error)
                  << std::endl;
        throw std::runtime_error("lodepng::encode() fail");
    }

    std::cout << "image saved to `" << file_name << "`\n\n" << std::flush;
}
