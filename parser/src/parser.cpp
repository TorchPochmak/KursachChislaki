

#include <fstream>
#include <iostream>
#include <sstream>

#include "../include/parser.h"



Calculation::reaction_data::reaction_data  () :
        _reactions({}),
        _mass({}),
        _fuel_fraction({}),
        _A({}),
        _n({}),
        _E({})
{   };


Calculation::reaction_data &Calculation::reaction_data::parse_reactions (const std::string &path) noexcept
{
    if (_reactions.size() != 0)
    {
        return *this;
    }

    std::ifstream input(path);
    if (!input.is_open())
    {
        std::cout << "File seems to be abscent or corrupted.\n";
    }

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
        _reactions.push_back(row_result);
    }

    return *this;
};


Calculation::reaction_data &Calculation::reaction_data::parse_coefficients (const std::string &path) noexcept
{
    std::ifstream input(path);
    if (!input.is_open())
    {
        std::cout << "File seems to be abscent or corrupted.\n";
        return *this;
    }

    std::string buff;
    std::vector< std::vector<double>> buff_result({});
    while (std::getline(input, buff))
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
    _cnt_components = buff_result[0].size();

    _A = buff_result[1];
    _n = buff_result[2];
    _E = buff_result[3];

    return *this;
};


Calculation::reaction_data &Calculation::reaction_data::parse_fuel (const std::string &path) noexcept
{
    std::ifstream input(path);

    if (!(input.is_open()))
    {
        std::cout << "File seems to abscent or corructed.\n";
        return *this;
    }

    double overall_mass = 0;
    double overall_cnt = 0;

    std::string buff;
    std::getline(input, buff);
    std::istringstream fuel(buff);

    double fraction = 0;
    while (fuel >> fraction)
    {
        overall_cnt += fraction;
        _fuel_fraction.push_back(fraction);
    }

    if (_fuel_fraction.size() != _cnt_components)
    {
        _fuel_fraction = {};
        return *this;
    }

    auto iter = _mass.begin();
    for (auto &elem : _fuel_fraction)
    {
        overall_mass += elem * (*iter);
        ++iter;
    }
    //overall_mass /= overall_cnt;

    for (auto &elem : _fuel_fraction)
    {
        elem = elem / overall_mass; // / overall_cnt;
    }

    return *this;
};



std::vector<double> Calculation::get_mole_count (const reaction_data &data, const double T) noexcept
{
    std::vector<double> result({});
    auto reaction_rate = get_reaction_rates(data, T);
    int sz = reaction_rate.size();

    for (int el = 0; el < data._cnt_components; el++)
    {
        double per_element = 0;
        for (int r = 0; r < sz; r++)
        {
            per_element += (data._reactions[r][el].real() - data._reactions[r + data._cnt_components][el].real()) *
                reaction_rate[r].second - reaction_rate[r].first;
        }
        result.push_back(per_element);
    }

    return result;
};


std::vector< std::pair<double, double>> Calculation::get_reaction_rates (const reaction_data &data,
    const double T) noexcept
{
    std::vector< std::pair<double, double>> result({});

    for (int i = 0; i < data._cnt_components; i++)
    {
        std::pair<double, double> per_reaction = std::make_pair(
            (get_reaction_density_comp_forward(data._cnt_components, data._reactions[i], data._fuel_fraction) *
                get_K_forward(data, T, i)),
            (get_reaction_density_comp_reverse(data._cnt_components, data._reactions[i], data._fuel_fraction) *
                get_K_reverse(data, T, i)));
        result.push_back(per_reaction);
    }

    return result;
};


double Calculation::get_reaction_density_comp_forward (
    const int &elem_cnt,
    const std::vector<std::complex<double>> &reaction, 
    const std::vector<double> &fuel_fractions) noexcept
{
    return get_reaction_density_comp (elem_cnt,
        reaction,
        fuel_fractions);
};


double Calculation::get_reaction_density_comp_reverse (
    const int &elem_cnt,
    const std::vector<std::complex<double>> &reaction, 
    const std::vector<double> &fuel_fractions) noexcept
{
    //std::vector< std::complex<double>> reaction_reverse = ;
    return get_reaction_density_comp (elem_cnt,
        {reaction.begin() + elem_cnt, reaction.end() - 1},
        fuel_fractions);
};


double Calculation::get_reaction_density_comp (
        const int &elem_cnt,
        const std::vector<std::complex<double>> &reaction, 
        const std::vector<double> &fuel_fractions) noexcept
{
    double result = 1;

    for (int i = 0; i < elem_cnt; i++)
    {
        result *= pow(ro * fuel_fractions[i], reaction[i].real());
    }

    return result;
};


double Calculation::get_K_forward (const reaction_data &data, const double &T, const int &r) noexcept
{
    return get_K(data._A[r], data._E[r], data._n[r], T);
};


double Calculation::get_K_reverse (const reaction_data &data, const double &T, const int &r) noexcept
{
    int shifted_r = r + data._cnt_components;
    return get_K(data._A[shifted_r], data._E[shifted_r], data._n[shifted_r], T);
};


double Calculation::get_K (const double &A,
        const double &E,
        const double &n,
        const double &T) noexcept
{
    return (A * pow(T, n) * exp(-(E / T)));
};

