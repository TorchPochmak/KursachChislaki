

#include <iostream>

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

    auto result = Calculation::get_matrix(reactions, 300);

}