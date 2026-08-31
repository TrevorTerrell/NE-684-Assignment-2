//
// Created by Trevor on 8/28/2026.
//

#pragma once
#import <cmath>
#import <vector>
#import "RandomManager.h"
#import "CrossSections.h"

#ifndef NE684A2_NEUTRON_H
#define NE684A2_NEUTRON_H

#define HYPER_WEIGHT_THRESH 1e-5
#define HYPER_ROUNDS 10

#define FINE_FLUX_GROUPS 1000

struct Global_Tallies {
    float N = 0.0f;
    float k_inf = 0.0f;
    float k_inf_2 = 0.0f;
    std::vector<double> group_flux = {0.0f, 0.0f, 0.0f};
    std::vector<double> group_flux_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr_2 = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> flux;
    std::vector<double> flux_2;
};

//Fission table will be a vector of pointers to Fission_Neutron(s)
//for asynchronicity, a neutron will compute its contribution to every tally then apply said contribution itself.

struct Fission_Neutron {
    float weight = -1.0f;
    // This would also hold position, but in infinite homogenous, the transport problem simplifies to 0D
};

struct Semilocal_Results {
    float N = 0.0f;
    float k_inf = 0.0f;
    float k_inf_2 = 0.0f;
    std::vector<double> group_flux = {0.0f, 0.0f, 0.0f};
    std::vector<double> group_flux_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr_2 = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr_2 = {0.0f, 0.0f, 0.0f};
    std::vector<double> flux;
    std::vector<double> flux_2;
};

struct Local_Results {
    float N = 0.0f;
    float k_inf = 0.0f;
    std::vector<double> group_flux = {0.0f, 0.0f, 0.0f};
    std::vector<double> scatter_rr = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    std::vector<double> capture_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> fission_rr = {0.0f, 0.0f, 0.0f};
    std::vector<double> flux;
};

class analog_neutron {
public:
    explicit analog_neutron(const Fission_Neutron n) {
        weight = n.weight;
        energy = RandomManager::getRandomFrac() * (ENERGY_MAX - ENERGY_MIN) + ENERGY_MIN;
    }

    void simulate(CrossSections *crossSections, Local_Results *results) {
        bool alive = true;
        while (alive) {
            std::vector<double> XSec = crossSections->getCrossSections(energy);

            //Test collision type
            float collision = static_cast<float>(RandomManager::getRandomFrac()) * XSec[7]; //ranges from 0-XSec_total
                //only the first 7 are actual cross section values
            uint32_t i;
            for (i = 0; i < 7; ++i) {
                    //strictly less than b/c 0.0 from the random number generator is inclusive (1.0 is exclusive)
                if (collision < XSec[i]) {
                    break;
                }
                collision -= XSec[i];
            }
            // Found collision in cross section i
            switch (i % 3) {
                case 0:
                    //Scatter
                    energy = RandomManager::getRandomFrac() * (energy - ENERGY_MIN) + ENERGY_MIN;
                    // apply to local tallies
                    break;
                case 1:
                    //Capture
                    // apply to local tallies
                    alive = false;
                    break;
                default:
                    //Fission
                    const float nu = crossSections->neutronsFromFission[i / 3] * weight; // int division to floor
                    //apply to local tallies
                    results->k_inf += nu;
                    alive = false;
                    break;
            }
        }
    }

private:
    double energy;
    float weight;
};

class implicit_neutron {
public:
    explicit implicit_neutron(const Fission_Neutron n) {
        weight = n.weight;
        energy = RandomManager::getRandomFrac() * (ENERGY_MAX - ENERGY_MIN) +ENERGY_MIN;
    }

    void simulate(CrossSections *crossSections, Local_Results *results) {
        bool alive = true;
        std::vector<double> XSec(8);
        std::vector<double> XSec_sums(3);
        double avg_nu = 0.0;
        while (alive) {
            XSec = crossSections->getCrossSections(energy);
            XSec_sums = {
                XSec[0] + XSec[3] + XSec[6],
                XSec[1] + XSec[4],
                XSec[2] + XSec[5],
            };
            avg_nu = (crossSections->neutronsFromFission[0] * XSec[2] +
                crossSections->neutronsFromFission[0] * XSec[5]) / XSec_sums[2];

            unsigned int g{};
            for (g = 0; g < energy_groups.size(); ++g) {
                if (energy >= energy_groups[g]) {
                    break;
                }
            }
            //std::cout << "g: " << g << "\n";

            //collision estimator of flux
            results->group_flux[g] += weight / XSec[7];
            fine_energy = static_cast<unsigned int>((FINE_FLUX_GROUPS - 1) * std::log(ENERGY_MAX / energy) / std::log(ENERGY_MAX / ENERGY_MIN));
            //std::cout << "flux_g: " << fine_energy << "\n";
            results->flux[fine_energy] += weight / XSec[7];

            //Fission Contribution
            results->k_inf += static_cast<float>(weight * avg_nu * XSec_sums[2] / XSec[7]);
            results->fission_rr[g] += weight * XSec_sums[2] / XSec[7];

            //Capture Contribution
            results->capture_rr[g] += weight * XSec_sums[1] / XSec[7];

            //Scatter Contribution
            new_energy = RandomManager::getRandomFrac() * (energy - ENERGY_MIN) + ENERGY_MIN;
            unsigned int g_new{};
            for (g_new = 0; g_new < energy_groups.size(); ++g_new) {
                if (new_energy >= energy_groups[g_new]) {
                    break;
                }
            }
            //std::cout << "g_new: " << g_new << "\n";

            g = g * (results->scatter_rr.size() - g + 1) / 2 + (g_new - g);
            results->scatter_rr[g] += weight * XSec_sums[0] / XSec[7];

            weight *= XSec_sums[0] / XSec[7];
            energy = new_energy;

            //Russian Roulette
            if (weight <= HYPER_WEIGHT_THRESH) {
                if (1.0 - RandomManager::getRandomFrac() <= 1.0 / HYPER_ROUNDS) {
                    alive = false;
                }
                // neutron lives or is dead and weight doesn't matter
                weight *= HYPER_ROUNDS;
            }
        }
    }

private:
    double energy;
    unsigned int fine_energy = 0;
    double new_energy{};
    double weight;

    const std::vector<double> energy_groups = {1e2, 1.0, ENERGY_MIN};
};

void runNeutronAnalog(CrossSections *crossSections, const std::vector<Fission_Neutron> *fission_bank, Local_Results *results) {
    Fission_Neutron neutron{.weight = 1.0f};
    if (!fission_bank->empty()) {
        const int index = static_cast<int>(RandomManager::getRandomFrac() * static_cast<double>(fission_bank->size()));
        neutron = fission_bank->at(index);
    }

    results->N += neutron.weight;

    analog_neutron n(neutron);
    n.simulate(crossSections, results);
}

void runNeutronImplicit(CrossSections *crossSections, const std::vector<Fission_Neutron> *fission_bank, Local_Results *results) {
    Fission_Neutron neutron{.weight = 1.0f};
    if (!fission_bank->empty()) {
        const int index = static_cast<int>(RandomManager::getRandomFrac() * static_cast<double>(fission_bank->size()));
        neutron = fission_bank->at(index);
    }

    results->N += neutron.weight;

    implicit_neutron n(neutron);
    n.simulate(crossSections, results);
}

bool exportTallies(Global_Tallies *tallies, const std::string& filename) {
    std::cout <<std::defaultfloat;

    std::cout << "\nFission Cross Sections:\n";
    for (unsigned int g = 0; g < tallies->group_flux.size(); ++g) {
        const auto fiss_xsec = tallies->fission_rr[g] / tallies->group_flux[g];
        const auto fiss_var = static_cast<float>(std::abs(std::pow(tallies->fission_rr[g] / tallies->group_flux[g], 2.0) - tallies->fission_rr_2[g] / tallies->group_flux[g]) / (tallies->N - 1));
        std::cout << g << ":\t" << fiss_xsec << " +/- " << std::sqrt(fiss_var) << "\n";
    }
    std::cout << "\nCapture Cross Sections:\n";
    for (unsigned int g = 0; g < tallies->group_flux.size(); ++g) {
        const auto cap_xsec = tallies->capture_rr[g] / tallies->group_flux[g];
        const auto cap_var = static_cast<float>(std::abs(std::pow(tallies->capture_rr[g] / tallies->group_flux[g], 2.0) - tallies->capture_rr_2[g] / tallies->group_flux[g]) / (tallies->N - 1));
        std::cout << g << ":\t" << cap_xsec << " +/- " << std::sqrt(cap_var) << "\n";
    }
    std::cout << "\nScatter Cross Sections:\n";
    for (int g = 0; g < tallies->scatter_rr.size(); ++g) {
        int e_g = 0;
        if (g > 2)
            e_g = 1;
        if (g > 4)
            e_g = 2;

        const auto scat_xsec = tallies->scatter_rr[g] / tallies->group_flux[e_g];
        const auto scat_var = static_cast<float>(std::abs(std::pow(tallies->scatter_rr[g] / tallies->group_flux[e_g], 2.0) - tallies->scatter_rr_2[g] / tallies->group_flux[e_g]) / (tallies->N - 1));

        std::cout << e_g << "->" << g - e_g * (7 - e_g) / 2 + e_g << ":\t" << scat_xsec << " +/- " << std::sqrt(scat_var) << "\n";
    }
    std::cout << "\nFew-Group Flux:\n";
    const std::vector<double> energy_groups = {ENERGY_MAX, 1e2, 1.0, ENERGY_MIN};
    for (int g = 0; g < tallies->group_flux.size(); ++g) {
        tallies->group_flux[g] /= (energy_groups[g] - energy_groups[g + 1]);
        tallies->group_flux_2[g] /= (energy_groups[g] - energy_groups[g + 1]);

        const auto flux_var = static_cast<float>(std::abs(std::pow(tallies->group_flux[g] / tallies->N, 2.0) - tallies->group_flux_2[g] / tallies->N) / (tallies->N - 1));

        std::cout << g << ":\t" << tallies->group_flux[g] / tallies->N << " +/- " << std::sqrt(flux_var) << "\n";
    }

    std::vector<std::vector<float>> flux(FINE_FLUX_GROUPS);
    std::vector<float> single_group_flux(3);
    double energy;
    double flux_var;
    double delta_energy;
    for (int g = 0; g < FINE_FLUX_GROUPS; ++g) {
        energy = std::log(ENERGY_MAX) - g * std::log(ENERGY_MAX / ENERGY_MIN) / (FINE_FLUX_GROUPS - 1);
        delta_energy = std::exp(energy - std::log(ENERGY_MAX / ENERGY_MIN) / (FINE_FLUX_GROUPS - 1)) - energy;
        energy = std::exp(energy);
        flux_var = std::abs(std::pow(tallies->flux[g] / tallies->N, 2.0) - tallies->flux_2[g] / tallies->N) / (tallies->N - 1);
        single_group_flux[0] = static_cast<float>(energy);
        single_group_flux[1] = static_cast<float>(tallies->flux[g] / delta_energy);
        single_group_flux[2] = static_cast<float>(std::sqrt(flux_var / delta_energy));

        flux[g] = single_group_flux;
    }

    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Unable to open file " << filename << "\n";
        return false;
    }
    outfile.clear();

    for (const auto &row : flux) {
        for (size_t i = 0; i < row.size(); ++i) {
            outfile << row[i];

            if (i < row.size() - 1)
                outfile << ",";

        }
        outfile << "\n";
    }

    outfile.close();

    return true;
}

#endif //NE684A2_NEUTRON_H
