#include <iostream>
#include <mutex>
#include <random>
#include <stack>
#include <thread>
#include <vector>

struct Vec3 {
    double x, y, z;
    Vec3(double _x = 0, double _y = 0, double _z = 0) {
        x = _x;
        y = _y;
        z = _z;
    }

    Vec3 operator+(const Vec3 &b) const { return Vec3(x + b.x, y + b.y, z + b.z); }

    Vec3 operator-(const Vec3 &b) const { return Vec3(x - b.x, y - b.y, z - b.z); }

    Vec3 operator*(double b) const { return Vec3(x * b, y * b, z * b); }

    Vec3 pointwise_product(const Vec3 &b) const { return Vec3(x * b.x, y * b.y, z * b.z); }

    Vec3 &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }

    double dot(const Vec3 &b) const { return x * b.x + y * b.y + z * b.z; } // cross:

    Vec3 operator%(const Vec3 &b) {
        return Vec3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

struct Ray {
    Vec3 o, d;
    Ray(const Vec3 &o_, const Vec3 &d_) : o(o_), d(d_) {}
};

enum Refl_t { DIFF, SPEC, REFR }; // material types, used in radiance()

struct Sphere {
    double rad; // radius

    Vec3 p, e, c; // position, emission, color

    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)

    Sphere(double rad_, const Vec3 &p_, const Vec3 &e_, const Vec3 &c_, Refl_t refl_)
        : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

    double intersect(const Ray &r) const { // returns distance, 0 if nohit
        Vec3 op = p - r.o;                 // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        double t;
        double eps = 1e-4, b = op.dot(r.d);
        double det = b * b - op.dot(op) + rad * rad;
        if (det < 0) {
            return 0;
        }

        det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};

std::vector<Sphere> spheres = {
    // Scene: radius, position, emission, color, material
    Sphere(1e5, Vec3(1e5 + 1, 40.8, 81.6), Vec3(), Vec3(.75, .25, .25),
           DIFF), // Left
    Sphere(1e5, Vec3(-1e5 + 99, 40.8, 81.6), Vec3(), Vec3(.25, .25, .75),
           DIFF),                                                        // Rght
    Sphere(1e5, Vec3(50, 40.8, 1e5), Vec3(), Vec3(.75, .75, .75), DIFF), // Back
    Sphere(1e5, Vec3(50, 40.8, -1e5 + 170), Vec3(), Vec3(), DIFF),       // Frnt
    Sphere(1e5, Vec3(50, 1e5, 81.6), Vec3(), Vec3(.75, .75, .75), DIFF), // Botm
    Sphere(1e5, Vec3(50, -1e5 + 81.6, 81.6), Vec3(), Vec3(.75, .75, .75),
           DIFF),                                                         // Top
    Sphere(16.5, Vec3(27, 16.5, 47), Vec3(), Vec3(1, 1, 1) * .999, SPEC), // Mirr
    Sphere(16.5, Vec3(73, 16.5, 78), Vec3(), Vec3(1, 1, 1) * .999, REFR), // Glas
    Sphere(600, Vec3(50, 681.6 - .27, 81.6), Vec3(12, 12, 12), Vec3(),
           DIFF) // Lite
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

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }

inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

inline bool intersect(const Ray &r, double &t, int &id) {
    double inf = t = std::numeric_limits<double>::infinity();
    for (int i = int(spheres.size()); i--;) {
        double d;
        if ((d = spheres[i].intersect(r)) && d < t) {
            t = d;
            id = i;
        }
    }
    return t < inf;
}

Vec3 radiance(const Ray &r, int depth, Sampler &sampler) {
    double t;   // distance to intersection
    int id = 0; // id of intersected object
    if (!intersect(r, t, id)) {
        return Vec3(); // if miss, return black
    }

    const Sphere &obj = spheres[id]; // the hit object
    Vec3 x = r.o + r.d * t;
    Vec3 n = (x - obj.p).norm();
    Vec3 nl = n.dot(r.d) < 0 ? n : n * -1;
    Vec3 f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl

    if (++depth > 5) {
        if (sampler.generate() < p) {
            f = f * (1 / p);
        } else {
            return obj.e; // R.R.
        }
    }

    if (obj.refl == DIFF) { // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * sampler.generate();
        double r2 = sampler.generate();
        double r2s = sqrt(r2);
        Vec3 w = nl;
        Vec3 u = ((fabs(w.x) > .1 ? Vec3(0, 1) : Vec3(1)) % w).norm();
        Vec3 v = w % u;
        Vec3 d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.pointwise_product(radiance(Ray(x, d), depth, sampler));
    }

    if (obj.refl == SPEC) { // Ideal SPECULAR reflection
        return obj.e +
               f.pointwise_product(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, sampler));
    }

    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) { // Total internal reflection
        return obj.e + f.pointwise_product(radiance(reflRay, depth, sampler));
    }

    Vec3 tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P,
           TP = Tr / (1 - P);

    return obj.e + f.pointwise_product(depth > 2
                                           ? (sampler.generate() < P
                                                  ? // Russian roulette
                                                  radiance(reflRay, depth, sampler) * RP
                                                  : radiance(Ray(x, tdir), depth, sampler) * TP)
                                           : radiance(reflRay, depth, sampler) * Re +
                                                 radiance(Ray(x, tdir), depth, sampler) * Tr);
}

void render(std::vector<Vec3> &pixels, std::stack<int> &job_list, std::mutex &mtx, int samples,
            int width, int height) {
    Ray cam(Vec3(50, 52, 295.6), Vec3(0, -0.042612, -1).norm()); // cam pos, dir
    Vec3 cx = Vec3(width * .5135 / height);
    Vec3 cy = (cx % cam.d).norm() * .5135;
    Vec3 r;

    /*
    rng.seed(0);
    double currentRandomNumber = distribution_0_1(rng);
    */

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

        for (unsigned short x = 0; x < width; x++) { // Loop cols
            sampler.seed(y * width + x);

            for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++) { // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec3()) {                   // 2x2 subpixel cols
                    for (int s = 0; s < samples; s++) {
                        double r1 = 2 * sampler.generate();
                        double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);

                        double r2 = 2 * sampler.generate();
                        double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

                        Vec3 d = cx * (((sx + .5 + dx) / 2 + x) / width - .5) +
                                 cy * (((sy + .5 + dy) / 2 + y) / height - .5) + cam.d;

                        r = r +
                            radiance(Ray(cam.o + d * 140, d.norm()), 0, sampler) * (1. / samples);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    pixels[i] = pixels[i] + Vec3(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
            }
        }
    }
}

int main() {
    int width = 1024;
    int height = 768;

    int samples = 20;
    auto pixels = std::vector<Vec3>(width * height);

    std::stack<int> job_list;
    for (int y = 0; y < height; y++) {
        job_list.push(y);
    }

    std::mutex mtx;
    std::vector<std::thread> threads;
    for (int idx = 0; idx < std::thread::hardware_concurrency() / 2; ++idx) {
        threads.push_back(std::thread(render, std::ref(pixels), std::ref(job_list), std::ref(mtx),
                                      samples / 4, width, height));
    }

    for (auto &t : threads) {
        t.join();
    }

    std::string file_name = "cosine_sampling_" + std::to_string(samples / 4 * 4) + ".ppm";

    FILE *f = fopen(file_name.c_str(), "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", width, height, 255);
    for (int i = 0; i < width * height; i++) {
        fprintf(f, "%d %d %d ", toInt(pixels[i].x), toInt(pixels[i].y), toInt(pixels[i].z));
    }
    std::cout << "output to `" << file_name << "`\n\n" << std::flush;
}
