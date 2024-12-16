

#include <fstream>
#include <iostream>
#include <sstream>

#include "../include/parser.h"



Calculation::reaction_data ::reaction_data  (const std::string &path_reactions,
    const std::string &path_coefs) :
        _reactions({}),
        _mass({}),
        _A({}),
        _n({}),
        _E({})
{
    std::ifstream input_reactions(path_reactions);

    if (!input_reactions.is_open())
    {
        std::cout << "File seems to be abscent or corrupted.\n";
    }

    std::string buff;

    while (std::getline(input_reactions, buff))
    {
        std::vector< std::complex<double>> row_result({});
        std::istringstream buff_row(buff);
        double value = 0;

        while (buff_row >> value)
        {
            row_result.push_back(std::complex<double>(value, 0));
        }
        
        _reactions.push_back(row_result);
    }

    std::ifstream input_coefs(path_coefs);

    if (!input_coefs.is_open())
    {
        std::cout << "File seems to be abscent or corrupted.\n";
        return;
    }

    std::vector< std::vector<double>> buff_result({});
    while (std::getline(input_coefs, buff))
    {
        std::vector<double> line_result({});
        double value = 0;

        std::istringstream line_buff(buff);

        while (line_buff >> value)
        {
            line_result.push_back(value);
        }

        buff_result.push_back(line_result);
    }

    _mass = buff_result[0];
    _A = buff_result[1];
    _n = buff_result[2];
    _E = buff_result[3];
};



std::vector< std::pair< double, double>> Calculation::get_chemical_rate (const double &T, const reaction_data &data) noexcept
{
    std::vector< std::pair<double, double>> result({});

    int dim = data._reactions[0].size() / 2;
    for (int i = 0; i < dim; i++)
    {
        std::vector<double> gamma = get_gamma((data._reactions)[i], data._mass);
        
        double comp_f = 1;
        double comp_r = 1;
        for (int j = 0; j < dim; j++)
        {
            comp_f *= pow(ro * gamma[j], data._reactions[i][j].real());
            comp_r *= pow(ro * gamma[j], data._reactions[i][j + dim].real());
        }

        result.push_back(std::make_pair(get_k(T, data._A[i], data._n[i], data._E[i]) * comp_f, 
            get_k(T, data._A[i + dim], data._n[i + dim], data._E[i + dim]) * comp_r));
    }

    return result;
};


double Calculation::get_k (const double &T,
    const double &A, 
    const double &n,
    const double &E) noexcept
{
    return (A * pow(T, n) * exp(-(E/T)));
};


std::vector<double> Calculation::get_gamma (const std::vector< std::complex<double>> &reaction,
        const std::vector<double> &mass) noexcept
{
    std::vector<double> result({});

    std::vector<double> parts({});
    int dim = reaction.size() / 2;

    double total = 0;
    for (int i = 0; i < dim; i++)
    {
        total += reaction[i].real();
    }

    double sum = 0;
    for (int i = 0; i < dim; i++)
    {
        parts.push_back(reaction[i].real() / total);
        sum += parts[i] * mass[i];
    }

    for (int i = 0; i < dim; i++)
    {
        result.push_back(parts[i] / sum);
    }

    return result;
};


