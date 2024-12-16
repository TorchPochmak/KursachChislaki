

#include <iostream>

#include "../parser/include/parser.h"


int main (int argc, char **argv)
{
    std::string path_reaction = "C:\\Users\\Asus\\Desktop\\Coding\\NumericalsCourseWork\\KursachChislaki\\test\\reaction.txt";
    std::string path_coefs = "C:\\Users\\Asus\\Desktop\\Coding\\NumericalsCourseWork\\KursachChislaki\\test\\mass.txt";

    Calculation::reaction_data data(path_reaction, path_coefs);

    auto result = Calculation::get_chemical_rate(300, data);

    std::cout << "gfy";

}