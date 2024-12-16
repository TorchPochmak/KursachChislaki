#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <complex>


const double EPS = 1e-5;
using namespace std;
class Matrix {
public:
	Matrix(const std::vector<std::vector<double>>& data) : data(data) {
		initialise();
	}

	Matrix(const Matrix& matrix) : data(matrix.data) {
		initialise();
	}

	Matrix(const uint64_t& row_count, const uint64_t& column_count) {
		data.assign(row_count, std::vector<double>(column_count, 0.0));
		initialise();
	}

	Matrix add(const Matrix& other) {
		if ( other.row_count != row_count || other.column_count != column_count ) {
			throw std::invalid_argument("Incompatible matrix dimensions.");
		}

		Matrix result(data);
		for ( int32_t i = 0; i < row_count; ++i ) {
			for ( int32_t j = 0; j < column_count; ++j ) {
				result.data[i][j] = data[i][j] + other.data[i][j];
			}
		}
		return result;
	}

	Matrix multiply(const Matrix& other) {
		if ( column_count != other.row_count ) {
			throw std::invalid_argument("Incompatible matrix dimensions.");
		}

		Matrix result(row_count, other.column_count);
		for ( int32_t i = 0; i < row_count; ++i ) {
			for ( int32_t j = 0; j < other.column_count; ++j ) {
				for ( int32_t k = 0; k < row_count; k++ ) {
					result.data[i][j] += data[i][k] * other.data[k][j];
				}
			}
		}
		return result;
	}

	Matrix transpose() {
		Matrix result(column_count, row_count);
		for ( int32_t i = 0; i < row_count; ++i ) {
			for ( int32_t j = 0; j < column_count; ++j ) {
				result.data[j][i] = data[i][j];
			}
		}
		return result;
	}

	Matrix minor(const int32_t& index) {
		Matrix result(row_count, column_count);
		for ( int32_t i = 0; i < index; ++i ) {
			result.set_entry(i, i, 1.0);
		}

		for ( int32_t i = index; i < row_count; ++i ) {
			for ( int32_t j = index; j < column_count; ++j ) {
				result.set_entry(i, j, data[i][j]);
			}
		}
		return result;
	}

	Matrix column(const int32_t& index) {
		Matrix result(row_count, 1);
		for ( int32_t i = 0; i < row_count; ++i ) {
			result.set_entry(i, 0, data[i][index]);
		}
		return result;
	}

	Matrix scalarMultiply(const double& value) {
		if ( column_count != 1 ) {
			throw std::invalid_argument("Incompatible matrix dimension.");
		}

		Matrix result(row_count, column_count);
		for ( int32_t i = 0; i < row_count; ++i ) {
			result.data[i][0] = data[i][0] * value;
		}
		return result;
	}

	Matrix unit() {
		if ( column_count != 1 ) {
			throw std::invalid_argument("Incompatible matrix dimensions.");
		}

		const double the_magnitude = magnitude();
		Matrix result(row_count, column_count);
		for ( int32_t i = 0; i < row_count; ++i ) {
			result.data[i][0] = data[i][0] / the_magnitude;
		}
		return result;
	}

	double magnitude() {
		if ( column_count != 1 ) {
			throw std::invalid_argument("Incompatible matrix dimensions.");
		}

		double norm = 0.0;
		for ( int32_t i = 0; i < row_count; ++i ) {
			norm += data[i][0] * data[i][0];
		}
		return std::sqrt(norm);
	}

	void display(const std::string& title) {
		std::cout << title << std::endl;
		for ( int32_t i = 0; i < row_count; ++i ) {
			for ( int32_t j = 0; j < column_count; ++j ) {
				std::cout << std::setw(9) << std::fixed << std::setprecision(4) << data[i][j];
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	double get_entry(const int32_t& row, const int32_t& col) {
		return data[row][col];
	}

	void set_entry(const int32_t& row, const int32_t& col, const double& value) {
		data[row][col] = value;
	}

	int32_t get_row_count() {
		return row_count;
	}

	int32_t get_column_count() {
		return column_count;
	}

	void initialise() {
		row_count = data.size();
		column_count = data[0].size();
	}

	int32_t row_count;
	int32_t column_count;
	std::vector<std::vector<double>> data;
};

typedef std::pair<Matrix, Matrix> matrix_pair;

Matrix householder_factor(Matrix vector) {
	if ( vector.get_column_count() != 1 ) {
		throw std::invalid_argument("Incompatible matrix dimensions.");
	}

	const int32_t size = vector.data.size();
	Matrix result(size, size);
	for ( int32_t i = 0; i < size; ++i ) {
		for ( int32_t j = 0; j < size; ++j ) {
			result.set_entry(i, j, -2 * vector.get_entry(i, 0) * vector.get_entry(j, 0));
		}
	}

	for ( int32_t i = 0; i < size; ++i ) {
		result.set_entry(i, i, result.get_entry(i, i) + 1.0);
	}
	return result;
}

matrix_pair householder(Matrix matrix) {
	const int32_t row_count = matrix.get_row_count();
	const int32_t column_count = matrix.get_column_count();
	std::vector<Matrix> versions_of_Q;
	Matrix z(matrix);
	Matrix z1(row_count, column_count);

	for ( int32_t k = 0; k < column_count && k < row_count - 1; ++k ) {
		Matrix vectorE(row_count, 1);
		z1 = z.minor(k);
		Matrix vectorX = z1.column(k);
		double magnitudeX = vectorX.magnitude();
		if ( matrix.get_entry(k, k) > 0 ) {
			magnitudeX = -magnitudeX;
		}

		for ( int32_t i = 0; i < vectorE.data.size(); ++i ) {
			vectorE.set_entry(i, 0, ( i == k ) ? 1 : 0);
		}
		vectorE = vectorE.scalarMultiply(magnitudeX).add(vectorX).unit();
		versions_of_Q.emplace_back(householder_factor(vectorE));
		z = versions_of_Q[k].multiply(z1);
	}

	Matrix Q = versions_of_Q[0];
	for ( int32_t i = 1; i < column_count && i < row_count - 1; ++i ) {
		Q = versions_of_Q[i].multiply(Q);
	}

	Matrix R = Q.multiply(matrix);
	Q = Q.transpose();
	return matrix_pair(R, Q);
}

Matrix solve_upper_triangular(Matrix r, Matrix b) {
	const int32_t column_count = r.get_column_count();
	Matrix result(column_count, 1);

	for ( int32_t k = column_count - 1; k >= 0; --k ) {
		double total = 0.0;
		for ( int32_t j = k + 1; j < column_count; ++j ) {
			total += r.get_entry(k, j) * result.get_entry(j, 0);
		}
		result.set_entry(k, 0, ( b.get_entry(k, 0) - total ) / r.get_entry(k, k));
	}
	return result;
}

Matrix least_squares(Matrix vandermonde, Matrix b) {
	matrix_pair pair = householder(vandermonde);
	return solve_upper_triangular(pair.first, pair.second.transpose().multiply(b));
}

Matrix fit_polynomial(Matrix x, Matrix y, const int32_t& polynomial_degree) {
	Matrix vandermonde(x.get_column_count(), polynomial_degree + 1);
	for ( int32_t i = 0; i < x.get_column_count(); ++i ) {
		for ( int32_t j = 0; j < polynomial_degree + 1; ++j ) {
			vandermonde.set_entry(i, j, std::pow(x.get_entry(0, i), j));
		}
	}
	return least_squares(vandermonde, y.transpose());
}
bool check_valid(Matrix A)
{
    int n = A.data.size();
    for(int i = 0; i < n; ++i)
    {
        if(i != n - 1)
        {
            double cur = A.data[i][i];
            double down = A.data[i+1][i];

            for(int j = i + 2; j < n; j++)
            {
                if(fabs(A.data[j][i]) > EPS)
                    return false;
            }
            if(fabs(down) < EPS)
            {
                continue;
            }
            else
            {
                for(int j = i + 2; j < n; j++)
                {
                    if(fabs(A.data[j][i + 1]) > EPS)
                        return false;
                }
                i++;
            }
        }
    }
    return true;
}

vector<complex<double>> get_eigenvalues(Matrix A)
{
    vector<complex<double>> result = vector<complex<double>>();
    for(int i = 0; i < A.data.size(); i++)
    {
        if(i == A.data.size() - 1)
        {
            result.push_back(complex<double>(A.get_entry(i, i), 0));
            return result;
        }
        double cur = A.get_entry(i,i);//a
        double down = A.get_entry(i + 1, i);//c
        double right = A.get_entry(i, i + 1);//b
        double down_right = A.get_entry(i + 1, i + 1);//d
        if(fabs(down) < EPS)
        {
            result.push_back(complex<double>(cur, 0));
        }
        else
        {
            double a = 1;
            double b = - cur - down_right;
            double c = cur * down_right - down * right;
            complex<double> discr = sqrt(complex<double>(b*b - 4 * a * c));
            complex<double> lambda1 = (-b + discr) / (2 * a);
            complex<double> lambda2 = (-b - discr) / (2 * a);
            result.push_back(lambda1);
            result.push_back(lambda2);
            i++;
        }
    }
}
int main() {
    int n;
    cin >> n;
	auto data = vector<vector<double>>(n, vector<double>(n,0));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cin >> data[i][j];
        }
    }

	// Task 1
	Matrix A(data);
	A.display("Initial matrix A:");
    matrix_pair pair = householder(A);
    Matrix Q = pair.second;
    Matrix R = pair.first;
    for(int i = 0; ; i++)
    {
        pair = householder(A);
        Q = pair.second;
        R = pair.first;      
        A = Matrix(R.multiply(Q));
       // if(i % 100 == 0)
        //{
            //Q.display("Matrix Q:");
            //R.display("Matrix R:");
        
        if(check_valid(A))
        {
            cout << "Число итераций: " << i << '\n';
            break;
        }
    }
    A.display("Matrix A:");
    vector<complex<double>> gol = get_eigenvalues(A);
    for(int i = 0; i < gol.size(); ++i)
    {
        cout << gol[i].real() << "+" << gol[i].imag() << "i" << '\n';
    }
}