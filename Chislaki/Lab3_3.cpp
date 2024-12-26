#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

class LeastSquares {
private:
    std::vector<double> x, y; // Табличные данные
    int n;                   // Количество точек

    // Метод для вычисления суммы квадратов ошибок
    double calculateError(const std::vector<double>& coefficients, int degree) const {
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            double approx = 0.0;
            for (int j = 0; j <= degree; ++j) {
                approx += coefficients[j] * std::pow(x[i], j);
            }
            error += std::pow(y[i] - approx, 2);
        }
        return error;
    }

    // Решение системы линейных уравнений методом Гаусса
    std::vector<double> solveSystem(std::vector<std::vector<double>>& matrix) {
        int size = matrix.size();
        std::vector<double> result(size);

        for (int i = 0; i < size; ++i) {
            // Прямой ход
            for (int j = i + 1; j < size; ++j) {
                double factor = matrix[j][i] / matrix[i][i];
                for (int k = i; k <= size; ++k) {
                    matrix[j][k] -= factor * matrix[i][k];
                }
            }
        }

        // Обратный ход
        for (int i = size - 1; i >= 0; --i) {
            result[i] = matrix[i][size] / matrix[i][i];
            for (int j = i - 1; j >= 0; --j) {
                matrix[j][size] -= matrix[j][i] * result[i];
            }
        }

        return result;
    }

public:
    LeastSquares(const std::vector<double>& x, const std::vector<double>& y)
        : x(x), y(y), n(x.size()) {}

    // Метод для вычисления коэффициентов многочлена заданной степени
    std::vector<double> fit(int degree) {
        int size = degree + 1;
        std::vector<std::vector<double>> matrix(size, std::vector<double>(size + 1, 0.0));

        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                for (int k = 0; k < n; ++k) {
                    matrix[i][j] += std::pow(x[k], i + j);
                }
            }
            for (int k = 0; k < n; ++k) {
                matrix[i][size] += std::pow(x[k], i) * y[k];
            }
        }

        return solveSystem(matrix);
    }

    // Метод для вывода коэффициентов и ошибки
    void printResult(int degree) {
        std::vector<double> coefficients = fit(degree);
        double error = calculateError(coefficients, degree);

        std::cout << "Аппроксимирующий многочлен степени " << degree << ":\n";
        std::cout << "f(x) = ";
        for (int i = 0; i <= degree; ++i) {
            if (i > 0 && coefficients[i] >= 0) {
                std::cout << " + ";
            }
            std::cout << coefficients[i];
            if (i > 0) {
                std::cout << " * x^" << i;
            }
        }
        std::cout << "\nСумма квадратов ошибок: " << error << "\n\n";
    }
};



void solve()
{
    // let's go kid
     // Табличные данные
    std::vector<double> x = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
    std::vector<double> y = {1.0, 1.0032, 1.0512, 1.2592, 1.8192, 3.0};

    // Создаём объект для вычислений
    LeastSquares ls(x, y);

    // Аппроксимация многочленом 1-й степени
    ls.printResult(1);

    // Аппроксимация многочленом 2-й степени
    ls.printResult(2);

}
int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    solve();
}