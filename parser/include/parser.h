
#include <string>
#include <vector>
#include <complex>


namespace Calculation
{
    const double ro = 1;

    class reaction_data  final
    {
    public:

        std::vector< std::vector< std::complex<double>>> _reactions;
        std::vector<double> _mass;
        std::vector<double> _A;
        std::vector<double> _n;
        std::vector<double> _E;

    public:

        reaction_data (const std::string &path_reactions,
            const std::string &path_coefs);
    };

    std::vector< std::pair< double, double>> get_chemical_rate (const double &T, const reaction_data &data) noexcept;

    double get_k (const double &T,
        const double &A, 
        const double &n,
        const double &E) noexcept;

    std::vector<double> get_gamma (const std::vector< std::complex<double>> &reaction,
        const std::vector<double> &mass) noexcept;

    
  
}; // namespace Calculation

