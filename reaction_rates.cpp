#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "qr_algorithm.hpp"


const int SPECIES_COUNT = 7;
const int REACTION_COUNT = 7;
const double T = 3000.0;  // Temperature in K
const double RHO = 1000.0;  // Density in kg/m^3
const double R = 8.31446;  // Universal gas constant in J/(mol·K)

// Molecular masses in g/mol
const double M[SPECIES_COUNT] = {
    1.00784,    // H
    15.999,     // O
    17.00734,   // OH
    2.01568,    // H2
    31.998,     // O2
    18.01528,   // H2O
    28.014      // N2
};

class ReactionSystem {
private:
    std::vector<std::vector<double>> nu_forward;  // Forward stoichiometric coefficients
    std::vector<std::vector<double>> nu_backward; // Backward stoichiometric coefficients
    std::vector<std::vector<double>> kinetic_params; // A, n, E for forward and backward reactions

public:
    ReactionSystem() {
        nu_forward.resize(REACTION_COUNT, std::vector<double>(SPECIES_COUNT));
        nu_backward.resize(REACTION_COUNT, std::vector<double>(SPECIES_COUNT));
        kinetic_params.resize(REACTION_COUNT, std::vector<double>(6)); // 3 params for forward and 3 for backward
    }

    void readData(const std::string& forward_file,
                  const std::string& backward_file,
                  const std::string& kinetics_file) {
        // Read forward coefficients
        std::ifstream fin_forward(forward_file);
        for (int i = 0; i < REACTION_COUNT; ++i) {
            for (int j = 0; j < SPECIES_COUNT; ++j) {
                fin_forward >> nu_forward[i][j];
            }
        }

        // Read backward coefficients
        std::ifstream fin_backward(backward_file);
        for (int i = 0; i < REACTION_COUNT; ++i) {
            for (int j = 0; j < SPECIES_COUNT; ++j) {
                fin_backward >> nu_backward[i][j];
            }
        }

        // Read kinetic parameters
        std::ifstream fin_kinetics(kinetics_file);
        for (int i = 0; i < REACTION_COUNT; ++i) {
            for (int j = 0; j < 6; ++j) {
                fin_kinetics >> kinetic_params[i][j];
            }
        }
    }

    // Calculate reaction rate coefficient k
    double calculateK(double A, double n, double E) {
        return A * pow(T, n) * exp(-E /(R* T));
    }

    // Calculate gamma values (kg/mole)
    std::vector<double> calculateGamma(const std::vector<double>& rho) {
        std::vector<double> gamma(SPECIES_COUNT);
        for (int i = 0; i < SPECIES_COUNT; ++i) {
            gamma[i] = rho[i] * (M[i] / 1000.0);  // Convert g/mol to kg/mol and multiply by density
        }
        return gamma;
    }

    // Calculate the Jacobian matrix
    std::vector<std::vector<double>> calculateJacobian(const std::vector<double>& gamma) {
        std::vector<std::vector<double>> jacobian(SPECIES_COUNT, std::vector<double>(SPECIES_COUNT, 0.0));

        for (int i = 0; i < SPECIES_COUNT; ++i) {
            for (int j = 0; j < SPECIES_COUNT; ++j) {
                double sum = 0.0;

                for (int r = 0; r < REACTION_COUNT; ++r) {
                    // Calculate forward and backward rate coefficients
                    double kf = calculateK(kinetic_params[r][0],
                                        kinetic_params[r][1],
                                        kinetic_params[r][2]);
                    double kb = calculateK(kinetic_params[r][3],
                                        kinetic_params[r][4],
                                        kinetic_params[r][5]);

                    // Calculate forward and backward rates
                    double forward_rate = kf;
                    double backward_rate = kb;

                    for (int s = 0; s < SPECIES_COUNT; ++s) {
                        if (nu_forward[r][s] > 0) {
                            forward_rate *= pow(gamma[s], nu_forward[r][s]);
                        }
                        if (nu_backward[r][s] > 0) {
                            backward_rate *= pow(gamma[s], nu_backward[r][s]);
                        }
                    }

                    // Calculate partial derivatives
                    if (nu_forward[r][j] > 0) {
                        sum += (nu_backward[r][i] - nu_forward[r][i]) *
                               nu_forward[r][j] * forward_rate / gamma[j];
                    }
                    if (nu_backward[r][j] > 0) {
                        sum += (nu_forward[r][i] - nu_backward[r][i]) *
                               nu_backward[r][j] * backward_rate / gamma[j];
                    }
                }

                jacobian[i][j] = sum;
            }
        }

        return jacobian;
    }
};

int test() {
    std::vector<std::vector<double>> A = {
        {6, 8, 2, 6, 34, 13},
        {54, 76, 12, 43, 23, 12},
        {11, 65, 34, 98, 9, 25},
        {27, 59, 73, 51, 49, 23},
        {17, 39, 49, 27, 74, 94},
        {69, 44, 23, 47, 26, 39}
    };
    std::vector<std::complex<double>> eigenvalues = find_eigenvalues(A);
    for (int i = 0; i < 6; ++i) {
        std::cout << "lambda " << i + 1 << " = " << std::fixed << std::setprecision(6)
                 << eigenvalues[i].real() << " + " << eigenvalues[i].imag() << "i" << std::endl;
    }
    return 0;
}

int main() {
    //set locale ru
    std::setlocale(LC_ALL, "ru");
    ReactionSystem system;

    // Read data from files
    system.readData("forward.txt", "backward.txt", "kinetics.txt");

    // Initial densities
    std::vector<double> rho(SPECIES_COUNT, RHO/SPECIES_COUNT);  // Distribute density equally among species

    // Calculate gamma values
    auto gamma = system.calculateGamma(rho);
    for (int i = 0; i < SPECIES_COUNT; ++i) {
        std::cout << "Gamma[" << i << "]: " << gamma[i] << std::endl;
    }

    // Calculate Jacobian
    auto jacobian = system.calculateJacobian(gamma);
    
    // Print Jacobian matrix
    std::cout << "Jacobian matrix dW_i/dgamma_j:" << std::endl;
    bool has_valid_values = false;
    for (int i = 0; i < SPECIES_COUNT; ++i) {
        for (int j = 0; j < SPECIES_COUNT; ++j) {
            std::cout.precision(4);
            std::cout << std::setw(10) << jacobian[i][j] << "\t";
            if (std::isfinite(jacobian[i][j])) {
                has_valid_values = true;
            }
        }
        std::cout << std::endl;
    }

    if (!has_valid_values) {
        std::cout << "\nWarning: Jacobian matrix contains no valid finite values!" << std::endl;
        return 1;
    }

    // Calculate eigenvalues using QR algorithm
    std::cout << "\nEigenvalues of the Jacobian matrix:" << std::endl;
    auto eigenvalues = find_eigenvalues(jacobian);
    for (int i = 0; i < SPECIES_COUNT; ++i) {
        std::cout << "lambda" << i + 1 << " = " << std::fixed << std::setprecision(6)
                 << eigenvalues[i].real() << " + " << eigenvalues[i].imag() << "i" << std::endl;
    }

    // find min and max abs eigenvalues
    double min_abs_eigenvalue = std::abs(eigenvalues[0].real());
    double max_abs_eigenvalue = std::abs(eigenvalues[0].real());
    for (int i = 1; i < SPECIES_COUNT; ++i) {
        if(std::abs(eigenvalues[i].real()) != 0){
            min_abs_eigenvalue = std::min(min_abs_eigenvalue, std::abs(eigenvalues[i].real()));
        }

        max_abs_eigenvalue = std::max(max_abs_eigenvalue, std::abs(eigenvalues[i].real()));
    }
    std::cout << "Stiffness: " << max_abs_eigenvalue / min_abs_eigenvalue << std::endl;
    return 0;
}
