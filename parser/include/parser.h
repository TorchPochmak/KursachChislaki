
#include <string>
#include <vector>
#include <complex>


namespace Calculation
{
    const double ro = 1;

    class reaction_data  final
    {
    public:

        std::vector< std::vector<double>> _reactions;

        std::vector<double> _mass;
        int _cnt_components;

        std::vector<double> _fuel_fraction;
        std::vector<double> _A;
        std::vector<double> _n;
        std::vector<double> _E;

    public:

        reaction_data ();

    public:

        reaction_data &parse_reactions (const std::string &path) noexcept;
        reaction_data &parse_coefficients (const std::string &path) noexcept;
        reaction_data &parse_fuel (const std::string &path) noexcept;
    };

    std::vector< std::vector<double>> get_matrix (const reaction_data &data, const double T) noexcept;

    std::vector<double> get_mole_count (const reaction_data &data, const double T) noexcept;

    std::vector< std::pair<double, double>> get_reaction_rates (const reaction_data &data,
        const double T) noexcept;

    double get_reaction_density_comp_forward (
        const int &elem_cnt,
        const std::vector<double> &reaction, 
        const std::vector<double> &fuel_fractions) noexcept;

    double get_reaction_density_comp_reverse (
        const int &elem_cnt,
        const std::vector<double> &reaction, 
        const std::vector<double> &fuel_fractions) noexcept;

    double get_reaction_density_comp (
        const int &elem_cnt,
        const std::vector<double> &reaction, 
        const std::vector<double> &fuel_fractions) noexcept;



    double get_K_forward (const reaction_data &data, const double &T, const int &r) noexcept;
    double get_K_reverse (const reaction_data &data, const double &T, const int &r) noexcept;

    double get_K (const double &A,
        const double &E,
        const double &n,
        const double &T) noexcept;
    
  
}; // namespace Calculation

