#include <iostream>
#include <iomanip>
#include <vector>
#include <future>
#include <thread>
#include <cmath>
#include <chrono>
#include "Neutron.h"
#include "CrossSections.h"

int main() {
    constexpr unsigned int NUM_PARTICLES = 500000;

    unsigned int num_cores = std::thread::hardware_concurrency();
    if (num_cores == 0) num_cores = 4; //fallback

    const unsigned int particles_per_core = NUM_PARTICLES / num_cores;

    CrossSections continuous_xsec;

    {
        std::cout << "Running simulation with " << NUM_PARTICLES << " particles in analog mode...\n";

        std::vector<Fission_Neutron> fission_bank;
        std::vector<Fission_Neutron> fission_bank_new;
        Global_Tallies tallies;

        const auto start_time = std::chrono::high_resolution_clock::now();

        std::vector<std::future<Semilocal_Results>> futures;

        for (int i = 0; i < num_cores; ++i) {
            futures.push_back(std::async(std::launch::async, [&continuous_xsec, &fission_bank, particles_per_core]() {
                Semilocal_Results thread_results;
                for (unsigned int j = 0; j < particles_per_core; ++j) {
                    Local_Results results;
                    runNeutronAnalog(&continuous_xsec, &fission_bank, &results);

                    thread_results.N += results.N;
                    thread_results.k_inf += results.k_inf;
                    thread_results.k_inf_2 += std::pow(results.k_inf, 2.0f);
                }
                return thread_results;
            }));
        }

        for (auto& f : futures) {
            const Semilocal_Results results = f.get();

            tallies.N += results.N;
            tallies.k_inf += results.k_inf;
            tallies.k_inf_2 += results.k_inf_2;
        }

        fission_bank = fission_bank_new;
        fission_bank_new.clear();

        const auto end_time = std::chrono::high_resolution_clock::now();
        const auto delta_time = end_time - start_time;
        const auto dt = static_cast<double>(delta_time.count());

        const auto crit_var = static_cast<float>(std::abs(std::pow(tallies.k_inf / tallies.N, 2.0) - tallies.k_inf_2 / tallies.N) / (tallies.N - 1));
        std::cout << std::fixed << std::setprecision(5) << "k_inf = " << tallies.k_inf / tallies.N << " +/- " << std::sqrt(crit_var) << "\n";
        std::cout << "Elapsed time: " << dt / 1e6 << " ms \n";
        std::cout << "FOM: " << 1.0 / (dt / 1e9 * crit_var) << "\n";
    }
    std::cout << "\n";
    {
        std::cout << "Running simulation with " << NUM_PARTICLES << " particles in implicit mode...\n";

        std::vector<Fission_Neutron> fission_bank;
        std::vector<Fission_Neutron> fission_bank_new;
        Global_Tallies tallies;

        const auto start_time = std::chrono::high_resolution_clock::now();

        std::vector<std::future<Semilocal_Results>> futures;

        for (int i = 0; i < num_cores; ++i) {
            futures.push_back(std::async(std::launch::async, [&continuous_xsec, &fission_bank, particles_per_core]() {
                Semilocal_Results thread_results;
                for (unsigned int j = 0; j < particles_per_core; ++j) {
                    Local_Results results;
                    runNeutronImplicit(&continuous_xsec, &fission_bank, &results);

                    thread_results.N += results.N;
                    thread_results.k_inf += results.k_inf;
                    thread_results.k_inf_2 += std::pow(results.k_inf, 2.0f);
                }
                return thread_results;
            }));
        }

        for (auto& f : futures) {
            const Semilocal_Results results = f.get();

            tallies.N += results.N;
            tallies.k_inf += results.k_inf;
            tallies.k_inf_2 += results.k_inf_2;
        }

        fission_bank = fission_bank_new;
        fission_bank_new.clear();

        const auto end_time = std::chrono::high_resolution_clock::now();
        const auto delta_time = end_time - start_time;
        const auto dt = static_cast<double>(delta_time.count());

        const auto crit_var = static_cast<float>(std::abs(std::pow(tallies.k_inf / tallies.N, 2.0) - tallies.k_inf_2 / tallies.N) / (tallies.N - 1));
        std::cout << std::fixed << std::setprecision(5) << "k_inf = " << tallies.k_inf / tallies.N << " +/- " << std::sqrt(crit_var) << "\n";
        std::cout << "Elapsed time: " << dt / 1e6 << " ms \n";
        std::cout << "FOM: " << 1.0 / (dt / 1e9 * crit_var) << "\n";
    }

    return !ExportCrossSectionsToCSV(&continuous_xsec, "continuous_xsec.csv");
}
