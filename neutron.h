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
// k_eff optimal: 1e-3, 1
// general: 1e-5?, 10?
#define HYPER_WEIGHT_THRESH 1e-3
#define HYPER_ROUNDS INFINITY

#define FINE_FLUX_GROUPS 1000

/**
 * @brief Holds the tally data lumped from every thread.
 */
struct Global_Tallies {
    float N = 0.0f;
    float k_inf = 0.0f;
    float k_inf_2 = 0.0f;
};

//Fission table will be a vector of pointers to Fission_Neutron(s)
//for asynchronicity, a neutron will compute its contribution to every tally then apply said contribution itself.
/**
 * @brief Holds the necessary information to store a neutron in the fission bank.
 */
struct Fission_Neutron {
    float weight = -1.0f;
    // This would also hold position, but in infinite homogenous, the transport problem simplifies to 0D
};

/**
 * @brief Holds the tally information collected per-thread.
 */
struct Semilocal_Results {
    float N = 0.0f;
    float k_inf = 0.0f;
    float k_inf_2 = 0.0f;
};

/**
 * @brief Holds the tally contribution from a single neutron.
 */
struct Local_Results {
    float N = 0.0f;
    float k_inf = 0.0f;
};

/**
 * @breif Neutron object used for Monte Carlo simulation.
 */
class neutron {
public:
    /**
     * @brief Creates a neutron object from a banked neutron.
     * @param n A neutron from the fission bank.
     */
    explicit neutron(const Fission_Neutron n) {
        weight = n.weight;
        energy = RandomManager::getRandomFrac() * (ENERGY_MAX - ENERGY_MIN) +ENERGY_MIN;
    }

    /**
     * @breif Simulates this neutron via Implicit Capture with Russian Roulette.
     * @param crossSections A set of continuous energy macroscopic cross sections.
     * @param results A mutable struct of this neutrons tally contributions.
     */
    void simulate_implicit(const CrossSections *crossSections, Local_Results *results) {
        bool alive = true;
        std::vector<double> XSec(8);
        std::vector<double> XSec_sums(3);
        double avg_nu = 0.0;
        unsigned int g{};
        while (alive) {
            XSec = crossSections->getCrossSections(energy);
            XSec_sums = {
                XSec[0] + XSec[3] + XSec[6],
                XSec[1] + XSec[4],
                XSec[2] + XSec[5],
            };
            avg_nu = crossSections->neutronsFromFission[0] * XSec[2] +
                crossSections->neutronsFromFission[0] * XSec[5];

            //Fission Contribution
            results->k_inf += static_cast<float>(weight * avg_nu / XSec[7]);

            //Capture Contribution

            //Scatter Contribution
            weight *= XSec_sums[0] / XSec[7];

            //Russian Roulette
            if (weight <= HYPER_WEIGHT_THRESH) {
                if (1.0 - RandomManager::getRandomFrac() >= 1.0 / HYPER_ROUNDS) {
                    alive = false;
                }
                // neutron lives or is dead and weight doesn't matter
                weight *= HYPER_ROUNDS;
            }

            // doing this here to optimize the number of operations
            energy = RandomManager::getRandomFrac() * (energy - ENERGY_MIN) + ENERGY_MIN;
        }
    }

    /**
     * @breif Simulates this neutron via Analog Monte Carlo.
     * @param crossSections A set of continuous energy macroscopic cross sections.
     * @param results A mutable struct of this neutrons tally contributions.
     */
    void simulate_analog(const CrossSections *crossSections, Local_Results *results) {
        bool alive = true;
        std::vector<double> XSec(8);
        double collision;
        unsigned int i;

        unsigned int g{};
        while (alive) {
            XSec = crossSections->getCrossSections(energy);

            //determining the type of interaction
            collision = RandomManager::getRandomFrac() * XSec[7];
            for (i = 0; i < 7; ++i) {
                if (collision < XSec[i])
                    break;
                collision -= XSec[i];
            }

            //given cross section vector structure, the modulo can determine the type of interaction
            switch (i % 3) {
                case 0:
                    //Scatter
                    energy = RandomManager::getRandomFrac() * (energy - ENERGY_MIN) + ENERGY_MIN;
                    break;
                case 1:
                    //Capture
                    alive = false;
                    break;
                default:
                    const double nu = crossSections->neutronsFromFission[i / 3];
                    results->k_inf += static_cast<float>(weight * nu);
                    alive = false;
                    break;
            }
        }
    }

private:
    double energy;
    double weight;

    const std::vector<double> energy_groups = {1e2, 1.0, ENERGY_MIN};
};

/**
 * @brief Handles the entire lifetime of a single neutron. Simulates with Analog Monte Carlo.
 * @param crossSections A set of continuous energy macroscopic cross sections.
 * @param fission_bank
 * @param results A mutable struct of this neutrons tally contributions.
 */
inline void runNeutronAnalog(const CrossSections *crossSections, const std::vector<Fission_Neutron> *fission_bank, Local_Results *results) {
    Fission_Neutron banked_neutron{.weight = 1.0f};
    if (!fission_bank->empty()) {
        const int index = static_cast<int>(RandomManager::getRandomFrac() * static_cast<double>(fission_bank->size()));
        banked_neutron = fission_bank->at(index);
    }

    results->N += banked_neutron.weight;

    neutron n(banked_neutron);
    n.simulate_analog(crossSections, results);
}

/**
 * @brief Handles the entire lifetime of a single neutron. Simulates with Implicit Capture.
 * @param crossSections A set of continuous energy macroscopic cross sections.
 * @param fission_bank
 * @param results A mutable struct of this neutrons tally contributions.
 */
inline void runNeutronImplicit(const CrossSections *crossSections, const std::vector<Fission_Neutron> *fission_bank, Local_Results *results) {
    Fission_Neutron banked_neutron{.weight = 1.0f};
    if (!fission_bank->empty()) {
        const int index = static_cast<int>(RandomManager::getRandomFrac() * static_cast<double>(fission_bank->size()));
        banked_neutron = fission_bank->at(index);
    }

    results->N += banked_neutron.weight;

    neutron n(banked_neutron);
    n.simulate_implicit(crossSections, results);
}


#endif //NE684A2_NEUTRON_H
