#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "qr_algorithm.hpp"
#include <bits/stdc++.h>
#include <complex>

int SPECIES_COUNT = 7;
const int REACTION_COUNT = 7;
double T = 3000.0;  // Temperature in K
const double RHO = 1000.0;  // Density in kg/m^3
const double R = 8.31446;  // Universal gas constant in J/(mol·K)

// Molecular masses in g/mol
const double M[7] = {
    2.01568 * 10e-3,    // H2
    31.998 * 10e-3,     // O2
    18.01528 * 10e-3,   // H2O
    17.00734 * 10e-3,   // OH
    1.00784 * 10e-3,    // H
    15.999 * 10e-3,     // O
    28.014 * 10e-3      // N2
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
    std::pair<vector<double>,vector<double>> calcKmatr()
    {
        vector<double> Kb = {};
        vector<double> Kf = {};
        for (int r = 0; r < REACTION_COUNT; ++r) {
            // Calculate forward and backward rate coefficients
            double kf = calculateK(kinetic_params[r][0],
                                kinetic_params[r][1],
                                kinetic_params[r][2]);
            double kb = calculateK(kinetic_params[r][3],
                                kinetic_params[r][4],
                                kinetic_params[r][5]);
            Kf.push_back(kf);
            Kb.push_back(kb);

        }
        return std::make_pair(Kf, Kb);
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
    // std::vector<std::complex<double>> eigenvalues = find_eigenvalues(A);
    // for (int i = 0; i < 6; ++i) {
    //     std::cout << "lambda " << i + 1 << " = " << std::fixed << std::setprecision(6)
    //              << eigenvalues[i].real() << " + " << eigenvalues[i].imag() << "i" << std::endl;
    // }
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


    //if (!has_valid_values) {
    //  std::cout << "\nWarning: Jacobian matrix contains no valid finite values!" << std::endl;
    //     return 1;
        //}


//ВОТ ЗДЕСЬ ТИПА ГАММЫ ОК ДА ОК ДА ДА ИЛИ НЕТ СПРАШИВАЮ ОК ДА ИЛИ СУКА НЕТ ДА ИЛИ НЕТ
    auto g = vector<double> {2./7. * M[0], 2./7. * M[1], 0, 0, 0, 4./7.* M[6]}; //system.calculateGamma(rho);
    auto [Kf, Kr] = system.calcKmatr();



    for(int i = 0; i < g.size(); i++)
        cout << g[i] << '\t';

    std::vector<std::pair<double,double>> res = std::vector<std::pair<double, double>>();
    for(double i = 300; i <= 6000; i+=300)
    {
        T = i;
//ВОТ ЗДЕСЬ ОГРОМНЫЙ ХУЙ
        std::vector<std::vector<double>> jacobian = std::vector<std::vector<double>>
        {
        {
            -Kf[0] * g[1] - Kf[2] * g[3] - Kf[3] * g[4],
            -Kf[0] * g[0],
            Kr[2] * g[4],
            Kr[0] * 2 * g[3] - Kf[2] * g[0] + Kr[3] * g[4],
            Kr[2] * g[2] + Kr[3] * g[3] + Kf[6] * 2 * g[4] * g[6],
            -Kf[3] * g[0],
            -Kr[6] * g[0] + Kf[6] * g[4] * g[4]
        },
        {
            -Kf[0] * g[1],
            -Kf[0] * g[0] - Kf[1] * g[4],
            0,
            2 * Kr[0] * g[3] + Kr[1] * g[4],
            -Kf[1] * g[1],
            Kr[1] * g[3],
            0
        },
        {
            Kf[2] * g[3],
            0,
            -Kr[2] * g[4] - Kr[4] * g[5] - Kr[5] * g[6],
            Kf[2] * g[0] + Kr[4] * 2 * g[3] + Kf[5] * g[4] * g[6],
            -Kr[2] * g[2] + Kf[5] * g[3] * g[6],
            -Kr[4] * g[2],
            Kf[5] * g[6]
        },
        {
            2 * Kf[0] * g[1] - Kr[3] * g[4] + Kr[5] * g[5],
            2 * Kf[0] * g[0],
            Kr[2] * g[4] + 2 * Kr[5] * g[5] + Kr[6] * g[6],
            -4 * Kr[0] * g[3] - Kr[4] * g[5] - Kf[0] * g[0] - Kr[3] * g[4] - 4 * Kr[5] * g[5] - Kf[5] * g[4] * g[6],
            Kf[1] * g[1] + Kr[2] * g[2] - Kr[4] * g[3] - Kf[5] * g[3] * g[6],
            -Kr[4] * g[3] + Kf[0] * g[0] + 2 * Kr[2] * g[2],
            Kr[6] * g[2]
        },
        {
            Kr[3] * g[4] + Kr[5] * g[5] + 2 * Kr[6] * g[6],
            -Kf[1] * g[4],
            -Kr[3] * g[4] + Kr[6] * g[6],
            Kr[5] * g[5] + Kr[0] * g[0] - Kr[3] * g[4] - Kf[3] * g[4] * g[6],
            -Kf[1] * g[1] - Kf[2] * g[2] - Kr[4] * g[3] - Kr[5] * g[3] * g[6] - 4 * Kr[5] * g[5] * g[6],
            Kr[4] * g[3] + Kr[0] * g[0],
            0
        },
        {
            -Kr[5] * g[5],
            Kr[3] * g[4],
            -Kr[5] * g[5],
            -Kr[5] * g[5] + Kr[3] * g[4] + 2 * Kr[0] * g[3],
            Kf[1] * g[1] + Kr[4] * g[3],
            -Kr[4] * g[3] - Kf[0] * g[0] - Kr[2] * g[2],
            0
        },
        {
            0,
            0,
            0,
            0,
            0,
            0,
            0
        }
        };
        std::cout << '\n';
        vector<vector<double>> jacob2 = vector<vector<double>>(6,vector<double>(6,0));
        for(int i = 0; i < 7; i++)
        {
            for(int j = 0; j < 7; j++)
            {
                std::cout << scientific << showpoint << setprecision(2) << jacobian[i][j] << "\t";
            }
            std::cout << '\n';
        }
        std::cout << "T: " << T << '\n';
        for(int k = 1; k <= 6; k++)
            std::cout << "gamma " << k << "= " << '\n';
        std::cout << "coef: " << '\n';

        //ТЫ МЕНЯ УВАЖАЕШЬ ИЛИ НЕТ УВАЖАЕШЬ ИЛИ НЕТ

        // for (int i = 0; i < SPECIES_COUNT; ++i) {
        //     std::cout << "Gamma[" << i << "]: " << gamma[i] << std::endl;
        // }

        // Calculate Jacobian
        //auto jacobian = system.calculateJacobian(gamma);
        
        // Print Jacobian matrix
        //std::cout << "Jacobian matrix dW_i/dgamma_j:" << std::endl;
       
        //Нужно вычеркнуть третий столбец и третью строку
        // std::cout << "Якобиан (c вычеркнутым столбцом и строкой): \n";
        // auto jacob2 = std::vector<std::vector<double>>(jacobian.size() - 1, std::vector<double>(jacobian.size() - 1, 0));
        // int l = 0;
        // int k = 0;
        // for (int i = 0; i < SPECIES_COUNT; ++i) {
        //     l = 0;
        //     for (int j = 0; j < SPECIES_COUNT; ++j) {
        //         if(i == 2)
        //             i += 1;
        //         if(j == 2)
        //             j += 1;
        //         if (std::isfinite(jacobian[i][j])) {
        //             has_valid_values = true;
        //         }
        //         if(std::fabs(jacobian[i][j]) < 1e-5)
        //             jacobian[i][j] = 0;
        //         jacob2[k][l] = jacobian[i][j];
        //         std::cout << std::scientific << std::setprecision(2) << std::showpoint << jacobian[i][j] << "\t";
        //         l++;
        //     }
        //     std::cout << '\n';
        //     k++;
        //     //  std::cout << std::endl;
        // }
        // Calculate eigenvalues using QR algorithm
        //std::cout << "\nEigenvalues of the Jacobian matrix:" << std::endl;
        // auto [eigenvalues, A] = get_eigens(jacobian);

        // for (int i = 0; i < SPECIES_COUNT; ++i) {
        //     std::cout << "lambda_" << i + 1 << std::scientific << std::setprecision(2) << std::showpoint << " = " << std::fixed << std::setprecision(6)
        //           << eigenvalues[i].real() << " + " << eigenvalues[i].imag() << "i" << std::endl;
        // }

        // // find min and max abs eigenvalues
        // double min_abs_eigenvalue = std::fabs(eigenvalues[0].real());
        // double max_abs_eigenvalue = std::fabs(eigenvalues[0].real());
        // for (int i = 1; i < SPECIES_COUNT; ++i) {
        //     if(std::abs(eigenvalues[i].real()) != 0){
        //         min_abs_eigenvalue = std::min(min_abs_eigenvalue, std::abs(eigenvalues[i].real()));
        //     }

        //     max_abs_eigenvalue = std::max(max_abs_eigenvalue, std::abs(eigenvalues[i].real()));
        // }
        // res.push_back({i, max_abs_eigenvalue / min_abs_eigenvalue});

        // std::cout << "T: "<< i << " " << "Stiffness: " << max_abs_eigenvalue / min_abs_eigenvalue << std::endl;
        // std::cout << "----------------------------------------------------------------------------------" << std::endl;
    }
    // std::cout << res.size() << std::endl;
    // for(int i = 0; i < res.size(); i++)
    //     std::cout << res[i].first << " ";
    // std::cout << std::endl;
    // for(int i = 0; i < res.size(); i++)
    //     std::cout << res[i].second << " ";
     return 0;
}
