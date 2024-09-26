#pragma once

#include <random>
#include <stack>
#include <vector>

class RNG {
  public:
    RNG(int seed_val) : distribution(std::uniform_real_distribution<double>(0.0, 1.0)) {
        generator.seed(seed_val);
    }

    void seed(uint val) { generator.seed(val); }

    double operator()() { return distribution(generator); }

  private:
    std::mt19937 generator;
    std::uniform_real_distribution<double> distribution;
};

struct PrimarySample {
    int modify_time;
    double value;

    PrimarySample() {
        modify_time = 0;
        value = NAN;
    }
};

struct PSSMLTSampler {
  private:
    RNG rng;
    int sample_idx;

    std::vector<PrimarySample> samples;
    std::stack<PrimarySample> backup_samples;

    inline double mutate_small_step(const double x) {
        const double s1 = 1.0 / 512.0;
        const double s2 = 1.0 / 16.0;
        const double r = rng();
        const double dx = s1 / (s1 / s2 + std::abs(2.0 * r - 1.0)) - s1 / (s1 / s2 + 1.0);
        if (r < 0.5) {
            const double x1 = x + dx;
            return (x1 < 1.0) ? x1 : x1 - 1.0;
        } else {
            const double x1 = x - dx;
            return (x1 < 0.0) ? x1 + 1.0 : x1;
        }
    }

  public:
    int global_time;
    int large_step;
    int large_step_time;

    PSSMLTSampler(int seed_val) : rng(seed_val) {
        global_time = 0;
        large_step = 0;
        large_step_time = 0;
        sample_idx = 0;
        samples.resize(128);
    }

    void seed(int val) { rng.seed(val); }

    void init_sample_idx() { sample_idx = 0; }

    void clear_backup_samples() { backup_samples = std::stack<PrimarySample>(); }

    void recover_samples() {
        int idx = sample_idx - 1;
        while (!backup_samples.empty()) {
            // recover sample stack
            samples[idx--] = backup_samples.top();
            backup_samples.pop();
        }
    }

    inline double next_sample() {
        if (samples.size() <= sample_idx) {
            samples.resize(samples.size() * 1.5);
        }

        if (samples[sample_idx].modify_time < global_time) {
            if (large_step > 0) {
                backup_samples.push(samples[sample_idx]);
                samples[sample_idx].modify_time = global_time;
                samples[sample_idx].value = rng();
            } else {
                // small step

                if (samples[sample_idx].modify_time < large_step_time) {
                    samples[sample_idx].modify_time = large_step_time;
                    samples[sample_idx].value = rng();
                }

                for (int idx = 0; idx < global_time - 1 - samples[sample_idx].modify_time; ++idx) {
                    // TODO: rewrite this part to mutate multiple steps in one function
                    samples[sample_idx].value = mutate_small_step(samples[sample_idx].value);
                }

                backup_samples.push(samples[sample_idx]);
                samples[sample_idx].value = mutate_small_step(samples[sample_idx].value);
                samples[sample_idx].modify_time = global_time;
            }
        }

        sample_idx++;
        return samples[sample_idx - 1].value;
    }
};
