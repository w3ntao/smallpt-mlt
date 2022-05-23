#include <iostream>
#include <math.h>
#include <mutex>
#include <stack>
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <vector>

struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z; // position, also color (r,g,b)
    Vec(double _x = 0, double _y = 0, double _z = 0) {
        x = _x;
        y = _y;
        z = _z;
    }

    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }

    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }

    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }

    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }

    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }

    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; } // cross:

    Vec operator%(const Vec &b) {
        return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
    }
};

struct Ray {
    Vec o, d;
    Ray(const Vec &o_, const Vec &d_) : o(o_), d(d_) {}
};

enum Refl_t { DIFF, SPEC, REFR }; // material types, used in radiance()

struct Sphere {
    double rad; // radius

    Vec p, e, c; // position, emission, color

    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)

    Sphere(double rad_, const Vec &p_, const Vec &e_, const Vec &c_, Refl_t refl_)
        : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}

    double intersect(const Ray &r) const { // returns distance, 0 if nohit
        Vec op = p - r.o;                  // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
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
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25),
           DIFF), // Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75),
           DIFF),                                                     // Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF), // Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),       // Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF), // Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75),
           DIFF),                                                      // Top
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC), // Mirr
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR), // Glas
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(),
           DIFF) // Lite
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

Vec radiance(const Ray &r, int depth, unsigned short *Xi) {
    double t;   // distance to intersection
    int id = 0; // id of intersected object
    if (!intersect(r, t, id)) {
        return Vec(); // if miss, return black
    }

    const Sphere &obj = spheres[id]; // the hit object
    Vec x = r.o + r.d * t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    Vec f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // max refl

    if (++depth > 5) {
        if (erand48(Xi) < p) {
            f = f * (1 / p);
        } else {
            return obj.e; // R.R.
        }
    }

    if (obj.refl == DIFF) { // Ideal DIFFUSE reflection
        double r1 = 2 * M_PI * erand48(Xi);
        double r2 = erand48(Xi);
        double r2s = sqrt(r2);
        Vec w = nl;
        Vec u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
        Vec v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    }

    if (obj.refl == SPEC) { // Ideal SPECULAR reflection
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    }

    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTION
    bool into = n.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) { // Total internal reflection
        return obj.e + f.mult(radiance(reflRay, depth, Xi));
    }

    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P,
           TP = Tr / (1 - P);

    return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ? // Russian roulette
                                           radiance(reflRay, depth, Xi) * RP
                                                       : radiance(Ray(x, tdir), depth, Xi) * TP)
                                    : radiance(reflRay, depth, Xi) * Re +
                                          radiance(Ray(x, tdir), depth, Xi) * Tr);
}

void render(std::vector<Vec> &pixels, std::stack<int> &job_list, std::mutex &mtx, int samples,
            int width, int height) {
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
    Vec cx = Vec(width * .5135 / height);
    Vec cy = (cx % cam.d).norm() * .5135;
    Vec r;

    while (true) {
        mtx.lock();
        if (job_list.empty()) {
            mtx.unlock();
            return;
        }
        auto y = job_list.top();
        job_list.pop();
        mtx.unlock();

        for (unsigned short x = 0, Xi[3] = {0, 0, y * y * y}; x < width; x++) { // Loop cols
            for (int sy = 0, i = (height - y - 1) * width + x; sy < 2; sy++) {  // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec()) {                     // 2x2 subpixel cols
                    for (int s = 0; s < samples; s++) {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / width - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / height - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samples);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    pixels[i] = pixels[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
            }
        }
    }
}

int main() {
    int width = 1024;
    int height = 768;

    for (auto samples : {10, 20, 30, 40, 50}) {
        auto pixels = std::vector<Vec>(width * height);

        std::stack<int> job_list;
        for (int y = 0; y < height; y++) {
            job_list.push(y);
        }

        std::mutex mtx;
        std::vector<std::thread> threads;
        for (int idx = 0; idx < std::thread::hardware_concurrency() / 2; ++idx) {
            threads.push_back(std::thread(render, std::ref(pixels), std::ref(job_list),
                                          std::ref(mtx), samples / 4, width, height));
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
}
