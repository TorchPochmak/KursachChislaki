

#include <vector>
#include <string>

namespace matrixes
{

    class vector;
    class square_matrix;

    class vector final
    {
    private:

        std::vector< double> _data;
        short _dimensionality;

    public:

        vector ();
        vector (std::vector< double> data);
        vector (short dimensionality);

    public:

        ~vector() = default;

        vector (const matrixes::vector &other);
        vector (matrixes::vector &&other);

        vector &operator= (const matrixes::vector &other);
        vector &operator= (matrixes::vector &&other);
    
    public:

        matrixes::vector &turn_into_identity () noexcept;
        matrixes::vector get_normalized () noexcept;
        matrixes::vector normalize () noexcept;

    public:

        void show_vector () const noexcept;

        double get_norm () const noexcept;

        std::vector< double> get_vector () const noexcept;
        std::vector< double> &get_vector (int) noexcept;

        short get_dimensionality () const noexcept;
        short &get_dimensionality (int) noexcept;

    public:

        bool is_orthogonal (const vector &other, double precision) const noexcept;

        double vec_composition (const vector &other) const noexcept;

    public:

        double operator*= (const vector &other) const noexcept;
        double operator* (const vector &other) const noexcept;

        vector &operator*= (const matrixes::square_matrix &other) noexcept;
        vector &operator* (const matrixes::square_matrix &other) noexcept;

        vector operator/= (const double div) noexcept;
        vector operator/ (const double &div) const noexcept;

    };

    class square_matrix final
    {

    private:
        
        std::vector< std::vector< double>> _data;
        short _dimensionality;

    public:

        square_matrix ();
        square_matrix (const std::vector< std::vector< double>> &matrix);
        square_matrix (const short &dimensionality);
        square_matrix (const std::string &path);

        matrixes::square_matrix qr_alg () const noexcept;   

    public:

        ~square_matrix () = default;    

        square_matrix (const square_matrix &other) noexcept;
        square_matrix (square_matrix &&other) noexcept;

        square_matrix &operator= (const square_matrix &other) noexcept;
        square_matrix &operator= (square_matrix &&other) noexcept;

    public:

        square_matrix &turn_into_identity_matrix () noexcept;
        square_matrix &transpose () noexcept;

        double get_det () const noexcept;

    public:

        short &get_dimensionality_c () noexcept;
        short get_dimensionality () const noexcept;
        void set_dimensionality (short dimensionality) noexcept;

        std::vector< std::vector< double>> &get_matrix_c () noexcept;
        std::vector< std::vector< double>> get_matrix () const noexcept;
        void set_matrix (std::vector< std::vector< double>> data) noexcept;

        const square_matrix &show_matrix (const double &width = 4) const noexcept;
        const square_matrix &dump_matrix (const std::string o_path) const noexcept;

    public:

        square_matrix &operator*= (const square_matrix &other) noexcept;
        square_matrix &operator* (const square_matrix &other) noexcept;

        matrixes::vector operator*= (const matrixes::vector &other) const noexcept;
        matrixes::vector operator* (const matrixes::vector &other) const noexcept;
    };


}
