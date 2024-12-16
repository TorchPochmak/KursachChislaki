
#include "../include/Matrix.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>


matrixes::vector::vector () :
    _data({}),
    _dimensionality(0)
{   };


matrixes::vector::vector (std::vector< double> data) :
    _data(data),
    _dimensionality(data.size())
{   };


matrixes::vector::vector (short dimensionality) :
    _dimensionality(dimensionality)
{
    std::vector< double> buff({});
    buff.resize(dimensionality);
    std::fill(buff.begin(), buff.end(), 0);
    _data = buff;
};


matrixes::vector::vector (const matrixes::vector &other) :
    _data(other.get_vector()),
    _dimensionality(other.get_dimensionality())
{   };


matrixes::vector::vector (matrixes::vector &&other) :
    _data(other.get_vector()),
    _dimensionality(other.get_dimensionality())
{   };


matrixes::vector &matrixes::vector::operator= (const matrixes::vector &other) 
{
    if (this != &other)
    {
        this->get_vector(1) = other.get_vector();
        this->get_dimensionality(1) = other.get_dimensionality();
    }
    return *this;
};


matrixes::vector &matrixes::vector::operator= (matrixes::vector &&other)
{
    if (this != &other)
    {
        this->get_vector(1) = other.get_vector();
        other.get_vector(1) = {};

        this->get_dimensionality(1) = other.get_dimensionality();
        other.get_dimensionality(1) = 0;
    }
    return *this;
};


matrixes::vector &matrixes::vector::turn_into_identity () noexcept
{
    for (auto &element : _data)
    {
        element = 1;
    }
    return *this;
};


matrixes::vector matrixes::vector::get_normalized () noexcept
{
    matrixes::vector buff(*this / this->get_norm());
    return buff;
}


matrixes::vector matrixes::vector::normalize () noexcept
{
    return *this /= this->get_norm();
}



void matrixes::vector::show_vector () const noexcept
{
    std::cout << "transposed vector: (";
    for (auto element : _data)
    {
        std::cout << " " << element;
    }
    std::cout << ")\n";
};


double matrixes::vector::get_norm () const noexcept
{
    double buff = (-1.0) * std::numeric_limits<double>::max();

    for (auto el : _data)
    {
        buff = std::max(buff, el);
    }

    return buff;
};


std::vector< double> matrixes::vector::get_vector () const noexcept
{
    return _data;
};


std::vector< double> &matrixes::vector::get_vector (int) noexcept
{
    return _data;
};

short matrixes::vector::get_dimensionality () const noexcept
{
    return _dimensionality;
};


short &matrixes::vector::get_dimensionality (int) noexcept
{
    return _dimensionality;
};


bool matrixes::vector::is_orthogonal (const matrixes::vector &other, double precision) const noexcept
{
    return ((*this * other) + std::numeric_limits<double>::epsilon() <= precision);
};


double matrixes::vector::operator*= (const matrixes::vector &other) const noexcept
{
    double result = 0;
    for (int i = 0; i < this->_dimensionality; i++)
    {
        result += (this->_data)[i] * (other._data)[i];
    }
    return result;
};


double matrixes::vector::operator* (const matrixes::vector &other) const noexcept
{
    return (*this *= other);
};


matrixes::vector &matrixes::vector::operator*= (const matrixes::square_matrix &other) noexcept
{

    std::vector< double> buff_result({});

    short dim = other.get_dimensionality();
    for (int column = 0; column < dim; column++)
    {
        double sum = 0;
        for (int elem = 0; elem < dim; elem++)
        {
            sum += (this->_data)[elem] * (other.get_matrix())[elem][column];
        }
        buff_result.push_back(sum);
    }

    this->get_vector(1) = buff_result;
    return *this;
};

matrixes::vector &matrixes::vector::operator* (const matrixes::square_matrix &other) noexcept
{
    return (*this *= other);
};


matrixes::vector matrixes::vector::operator/= (const double div) noexcept
{   
    for (auto &element : this->get_vector(1))
    {
        element /= div;
    }
    return *this;
};


matrixes::vector matrixes::vector::operator/ (const double &div) const noexcept
{
    std::vector< double> result({});

    for (auto element : this->get_vector())
    {
        result.push_back(element / div);
    }

    matrixes::vector buff(result);
    return buff;
};


matrixes::square_matrix::square_matrix () :
    _data({}),
    _dimensionality(0)
{   };


matrixes::square_matrix::square_matrix (const std::vector< std::vector< double>> &matrix) :
    _data(matrix),
    _dimensionality(matrix.size())
{   };


matrixes::square_matrix::square_matrix (const short &dimensionality) :
    _data({}),
    _dimensionality(dimensionality)
{
    for (int i = 0; i < dimensionality; i++)
    {
        std::vector< double> row;
        row.resize(dimensionality);
        std::fill(row.begin(), row.end(), 0);

        _data.push_back(row);
    }
};


matrixes::square_matrix::square_matrix (const std::string &path)
{
    std::ifstream input_data(path);
    if (!input_data.is_open())
    {
        std::cout << "File seems to be abscent or corrupted." << std::endl;
    }

    std::string buff;
    short dimensionality = 0;

    while (std::getline(input_data, buff))
    {
        std::istringstream buff_row(buff);
        std::vector< double> row = {};
        double buff_coef = 0;

        while (buff_row >> buff_coef)
        {
            row.push_back(buff_coef);
        }
        ++dimensionality;
        _data.push_back(row);
    }

    _dimensionality = dimensionality;
};


matrixes::square_matrix::square_matrix (const square_matrix &other) noexcept :
    _data(other.get_matrix()),
    _dimensionality(other.get_dimensionality())
{   };


matrixes::square_matrix::square_matrix (square_matrix &&other) noexcept :
    _data(other.get_matrix()),
    _dimensionality(other.get_dimensionality())
{   };


matrixes::square_matrix matrixes::square_matrix::qr_alg () const noexcept
{

    


};


matrixes::square_matrix &matrixes::square_matrix::operator= (const square_matrix &other) noexcept
{   
    if (this != &other)
    {
        this->set_matrix(other.get_matrix());
        this->set_dimensionality(other.get_dimensionality());
    }
    return *this;
};


matrixes::square_matrix &matrixes::square_matrix::operator= (square_matrix &&other) noexcept
{
    if (this != &other)
    {
        this->set_matrix(other.get_matrix());
        other.set_matrix({});

        this->set_dimensionality(other.get_dimensionality());
        other.set_dimensionality(0);
    }
    return *this;
};


matrixes::square_matrix &matrixes::square_matrix::turn_into_identity_matrix () noexcept
{
    for (int i = 0; i < _dimensionality; i ++)
    {
        _data[i][i] = 1;
    }

    return *this;
};


matrixes::square_matrix &matrixes::square_matrix::transpose () noexcept
{
    for (int row = 0; row < _dimensionality; row++)
    {
        for (int column = 0; column <= row; column++)
        {
            if (row != column)
            {
                std::swap(_data[row][column], _data[column][row]);
            }
        }
    }

    return *this;
};


double matrixes::square_matrix::get_det () const noexcept
{
    double result = 0;
    return (_data[0][0] * _data[1][1] - _data[0][1] * _data[1][0]);
};


const matrixes::square_matrix &matrixes::square_matrix::show_matrix (const double &w) const noexcept
{
    for (auto row : _data)
    {
        for (auto element : row)
        {
            std::cout << std::setw(w) << element << "| ";
        }
        std::cout << "\n";
    }
    std::cout << std::endl;
    
    return *this;
};


const matrixes::square_matrix &matrixes::square_matrix::dump_matrix (const std::string o_path) const noexcept
{
    std::ofstream output(o_path);

    if (!output.is_open())
    {
        std::cout << "File [ " << o_path << " ] can't be reached\n";
    }

    for (auto row : _data)
    {
        for (auto element : row)
        {
            output << element << " ";
        }
        output << "\n";
    }

    return *this;
};


short &matrixes::square_matrix::get_dimensionality_c () noexcept
{
    return _dimensionality;
};


std::vector< std::vector< double>> &matrixes::square_matrix::get_matrix_c () noexcept
{
    return _data;
};


short matrixes::square_matrix::get_dimensionality () const noexcept
{
    return _dimensionality;
};

void matrixes::square_matrix::set_dimensionality (short dimensionality) noexcept
{
    _dimensionality = dimensionality;
};


std::vector< std::vector< double>> matrixes::square_matrix::get_matrix () const noexcept
{
    return _data;
};

void matrixes::square_matrix::set_matrix (std::vector< std::vector< double>> data) noexcept
{
    _data = data;
};


matrixes::square_matrix &matrixes::square_matrix::operator*= (const matrixes::square_matrix &other) noexcept
{
    if (this != &other)
    {
        short dim = this->get_dimensionality();
        for (int row = 0; row < dim; row++)
        {
            std::vector< double> resulted_row({});
            for (int column = 0; column < dim; column++)
            {
                double sum = 0;
                for (int element_pos = 0; element_pos < dim; element_pos++)
                {
                    sum += (this->get_matrix())[row][element_pos] * (other.get_matrix())[element_pos][column];
                }
                resulted_row.push_back(sum);
            }
            (this->get_matrix_c())[row] = resulted_row;
        }
    }
    return *this;
};


matrixes::square_matrix &matrixes::square_matrix::operator* (const matrixes::square_matrix &other) noexcept
{
    return (*this) *= other;
};


matrixes::vector matrixes::square_matrix::operator*= (const matrixes::vector &other) const noexcept
{
    std::vector< double> result_row({});

    short dim = this->get_dimensionality();
    for (int row = 0; row < dim; row++)
    {
        
        double sum = 0;
        for (int element = 0; element < dim; element++)
        {
            sum += (this->get_matrix())[row][element] * (other.get_vector())[element];
        }
        result_row.push_back(sum);
    }
    matrixes::vector buff(result_row);
    return buff;
};


matrixes::vector matrixes::square_matrix::operator* (const matrixes::vector &other) const noexcept
{
    return (*this *= other);
};
