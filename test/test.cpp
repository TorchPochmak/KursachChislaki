

#include <iostream>

#include "qr.hpp"
#include "../parser/include/parser.h"


int main (int argc, char **argv)
{
    // std::string path_reaction = "C:\\Users\\Asus\\Desktop\\Coding\\NumericalsCourseWork\\KursachChislaki\\test\\reaction.txt";
    // std::string path_coefs = "C:\\Users\\Asus\\Desktop\\Coding\\NumericalsCourseWork\\KursachChislaki\\test\\mass.txt";

    std::string path_reaction = "X:\\Coding\\CPP\\KursachChislaki\\test\\reaction.txt";
    std::string path_coefs = "X:\\Coding\\CPP\\KursachChislaki\\test\\mass.txt";
    std::string path_fuel = "X:\\Coding\\CPP\\KursachChislaki\\test\\fuel.txt";

    Calculation::reaction_data reactions = Calculation::reaction_data();

    reactions.
        parse_reactions(path_reaction).
        parse_coefficients(path_coefs).
        parse_fuel(path_fuel);

    for (int T = 300; T <= 6000; T+=300)
    {
        auto result = Calculation::get_matrix(reactions, T);

        auto sz = result.size();
        for (int i = 0; i < sz; i++)
        {
            for (int j = 0; j < sz; j++)
            {
                std::cout << result[i][j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        auto [eigenvalues, A] = get_eigens(result);

        for (int i = 0; i < reactions._cnt_components; ++i) {
            std::cout << "lambda_" << i + 1 << std::scientific << std::setprecision(2) << std::showpoint << " = " << std::fixed << std::setprecision(6)
                  << eigenvalues[i].real() << " + " << eigenvalues[i].imag() << "i" << std::endl;
        }   
    }
}