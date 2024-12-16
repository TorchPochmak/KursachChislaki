

#include <fstream>
#include <iostream>
#include <sstream>

#include "parser.h"


std::vector< std::vector< std::complex<double>>> parsing_tool::parse (const std::string &path) noexcept
{
    std::ifstream input(path);

    if (!input.is_open())
    {
        std::cout << "File seems to be abscent or corrupted.\n";
        return {};
    }

    std::vector< std::vector< std::complex<double>>> result({});
    std::string buff;

    while (std::getline(input, buff))
    {
        std::vector< std::complex<double>> row_result({});
        std::istringstream buff_row(buff);
        double value = 0;

        while (buff_row >> value)
        {
            row_result.push_back(std::complex<double>(value, 0));
        }
        
        result.push_back(row_result);
    }

    return result;
}


