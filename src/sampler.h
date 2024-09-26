#pragma once

#include <chrono>
#include <random>

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
