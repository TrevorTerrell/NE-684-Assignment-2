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

#define HYPER_WEIGHT_THRESH 1e-3
#define HYPER_ROUNDS 1

struct Global_Tallies {
    float N = 0.0f;
    float k_inf = 0.0f;
    float k_inf_2 = 0.0f;
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
    std::vector<Fission_Neutron> neutrons;
};

struct Local_Results {
    float N = 0.0f;
    float k_inf = 0.0f;
    std::vector<Fission_Neutron> neutrons;
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
            std::vector<float> XSec = crossSections->getCrossSections(energy);
            dist = -std::logf(static_cast<float>(1.0 - RandomManager::getRandomFrac())) / XSec[7];

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
    float dist = 0.0f;
};

class implicit_neutron {
public:
    explicit implicit_neutron(const Fission_Neutron n) {
        weight = n.weight;
        energy = RandomManager::getRandomFrac() * (ENERGY_MAX - ENERGY_MIN) + ENERGY_MIN;
    }

    void simulate(CrossSections *crossSections, Local_Results *results) {
        bool alive = true;
        std::vector<float> XSec(8);
        std::vector<float> XSec_sums(3);
        float avg_nu = 0.0f;
        while (alive) {
            XSec = crossSections->getCrossSections(energy);
            XSec_sums = {
                XSec[0] + XSec[3] + XSec[6],
                XSec[1] + XSec[4],
                XSec[2] + XSec[5],
            };
            avg_nu = (crossSections->neutronsFromFission[0] * XSec[2] +
                crossSections->neutronsFromFission[0] * XSec[5]) / XSec_sums[2];

            //Fission Contribution
            results->k_inf += weight * avg_nu * XSec_sums[2] / XSec[7];
            //results->neutrons.push_back(Fission_Neutron(weight * avg_nu * XSec_sums[2] / XSec[7]));

            //Capture Contribution
                //none

            //Scatter Contribution
            energy = RandomManager::getRandomFrac() * (energy - ENERGY_MIN) + ENERGY_MIN;
            weight *= XSec_sums[0] / XSec[7];

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
    float weight;
    float dist = 0.0f;
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

#endif //NE684A2_NEUTRON_H
