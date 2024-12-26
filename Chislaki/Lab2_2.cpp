#include <cmath>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef unsigned long long ull;
typedef long double ld;
typedef pair<ll, ll> pll;

const ll INF = 1e18;
const int INF_int = 1e9;
const ll MOD = 1e9 + 7;
const double PI = acos(-1.0);

double __speedstart__;
#define SPEEDTEST_START __speedstart__ = clock();
#define SPEEDTEST_STOP cout << "Process finished in: " << (double)(clock() - __speedstart__)/CLOCKS_PER_SEC << " seconds";
#define fastboi ios::sync_with_stdio(false); cin.tie(0);

using matrix = vector<vector<double>>;

matrix operator-(const matrix& mat1, const matrix& mat2) {

    if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
        throw std::invalid_argument("Размеры матриц не совпадают.");
    }

    matrix result = mat1; 
    for (std::size_t i = 0; i < mat1.size(); ++i) {
        for (std::size_t j = 0; j < mat1[i].size(); ++j) {
            result[i][j] -= mat2[i][j]; 
        }
    }
    return result;
}
matrix operator&(const matrix& mat1, const matrix& mat2) {

    if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
        throw std::invalid_argument("Размеры матриц не совпадают.");
    }

    matrix result = mat1; 
    for (std::size_t i = 0; i < mat1.size(); ++i) {
        for (std::size_t j = 0; j < mat1[i].size(); ++j) {
            result[i][j] = ((mat2[i][j] + mat1[i][j]) >= 2); 
        }
    }
    return result;
}
matrix operator+(const matrix& mat1, const matrix& mat2) {

    if (mat1.size() != mat2.size() || mat1[0].size() != mat2[0].size()) {
        throw std::invalid_argument("Размеры матриц не совпадают.");
    }

    matrix result = mat1; 
    for (std::size_t i = 0; i < mat1.size(); ++i) {
        for (std::size_t j = 0; j < mat1[i].size(); ++j) {
            result[i][j] += mat2[i][j]; 
        }
    }
    return result;
}
matrix equations_system (matrix x)
{
    matrix f(2, vector<double>(1));

    f[0][0] = pow(x[0][0],2) + pow(x[1][0],2) - 4;
    f[1][0] = x[0][0] - exp(x[1][0]) + 2;

    return f;
}

matrix Jacobi_matrix (matrix x)
{
    matrix J(2, vector<double>(2));

    J[0][0] = 2 * x[0][0];
    J[0][1] = 2 * x[1][0];

    J[1][0] = 1;
    J[1][1] = -exp(x[1][0]);

    return J;
}
double norm(const matrix& a)
{
    double o = -INF_int;
    for(int i = 0; i < a.size(); i++)
    {
        double sm = 0;
        for(int j = 0; j < a[i].size(); j++)
        {
            sm += abs(a[i][j]);
        }
        o = max(sm, o);
    }
    return o;
}
matrix equivalent_equations_system (matrix x)
{
    matrix f(2, vector<double>(1));

    f[0][0] = sqrt(4 - x[1][0] * x[1][0]);
    f[1][0] = log(x[0][0] + 2);

    return f;
}

matrix derivative_equivalent_system (double x_1, double x_2)
{
    matrix J(2, vector<double>(2));

    J[0][0] = 0;
    J[0][1] = -x_2 / (sqrt(4 - x_2 * x_2));

    J[1][0] = 1 / (x_1 + 2);
    J[1][1] = 0;

    return J;
}

matrix inverse_matrix(matrix mat) {
    if (mat.size() != 2 || mat[0].size() != 2 || mat[1].size() != 2) {
        throw std::invalid_argument("Матрица должна быть размером 2x2.");
    }

    // Вычисляем детерминант
    double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    
    if (det == 0) {
        throw std::runtime_error("Матрица не имеет обратной (детерминант равен 0).");
    }

    // Вычисляем обратную матрицу
    double invDet = 1.0 / det;
    matrix invMat = {
        { mat[1][1] * invDet, -mat[0][1] * invDet },
        { -mat[1][0] * invDet, mat[0][0] * invDet }
    };

    return invMat;
}

matrix mult_matr_vec(matrix& a, vector<double>& b)
{
    int a_r = a.size();
    int b_r = b.size();

    int a_c = a[0].size();
    int b_c = 1;

    if(a_c != b_r)
    {
        throw invalid_argument("ERROR: ошибка умножения матриц - не совместны");
    }

    matrix res(a_r, vector<double>(b_c, 0));
    for(int i = 0; i < a_r; i++)
    {
        for(int j = 0; j < b_c; j++)
        {
            for(int k = 0; k < a_c; k++)
            {
                res[i][j] += a[i][k] * b[k];
            }
        }
    }
    return res;
}

matrix mult(matrix& mat1, matrix& mat2) {
    size_t rows1 = mat1.size();
    size_t cols1 = mat1.empty() ? 0 : mat1[0].size();
    size_t rows2 = mat2.size();
    size_t cols2 = mat2.empty() ? 0 : mat2[0].size();

    if (cols1 != rows2) {
        throw std::invalid_argument("Число столбцов первой матрицы должно быть равно числу строк второй матрицы.");
    }

    matrix result(rows1, std::vector<double>(cols2, 0.0));

    // Выполняем умножение матриц
    for (size_t i = 0; i < rows1; ++i) {
        for (size_t j = 0; j < cols2; ++j) {
            for (size_t k = 0; k < cols1; ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return result;
}

matrix Newton_method (const std::function<matrix(matrix)>& f,
                      const std::function<matrix(matrix)>& J,
                      matrix x, double EPS)
{
    matrix prev_x;

    do
    {
        prev_x = x;
        matrix h = J(x);
        matrix h1 = inverse_matrix(h);
        matrix h2= f(x);
        x = x - mult(h1,h2);

    } while (norm((x - prev_x)) > EPS);

    return x;
}

matrix simple_iteration_method (const std::function<matrix(matrix)>& f,
                                const std::function<matrix(double, double)>& d_f,
                                matrix x, double EPS, double a_1, double b_1, double a_2, double b_2)
{
    double q, q_1;

    matrix d(2, vector<double>(2));
    d = d_f(a_1, b_1);

    q = norm(d);

    d = d_f(a_1, b_2);
    q_1 = norm(d);
    if (q_1 > q) q = q_1;

    d = d_f(a_2, b_1);
    q_1 = norm(d);
    if (q_1 > q) q = q_1;

    d = d_f(a_2, b_2);
    q_1 = norm(d);
    if (q_1 > q) q = q_1;

    if (q >= 1)
    {
        throw std::runtime_error("Simple iteration method does not converge");
    }

    matrix prev_x;

    do
    {
        prev_x = x;

        x = f(x);

    } while (q / (1 - q) * norm((x - prev_x)) > EPS);

    return x;
}
int solve()
{
    double EPS = 1E-3;

    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    matrix x_1(2, vector<double>(1));
    x_1[0][0] = 0.25;
    x_1[1][0] = 0.25;
    matrix x_2 = x_1;

    std::cout << "Newton method" << std::endl;

    try
    {
        x_1 = Newton_method(equations_system, Jacobi_matrix, x_1, EPS);
    }
    catch (std::runtime_error &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;

        return 1;
    }

    for (int i = 0; i < x_1.size(); ++i)
    {
        std::cout << "x" << i + 1 << " = " << x_1[i][0] << std::endl;
    }
    cout << "Проверка первого: " << x_1[0][0] * x_1[0][0] + x_1[1][0] * x_1[1][0] - 4 << '\n';
    cout << "Проверка второго: " << x_1[0][0] - exp(x_1[1][0]) + 2 << '\n';
    std::cout << "Simple iteration method" << std::endl;

    try
    {
        x_2 = simple_iteration_method(equivalent_equations_system, derivative_equivalent_system, x_2, EPS, 0, 0.5, 0, 0.5);
    }
    catch (std::runtime_error &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;

        return 1;
    }

    for (int i = 0; i < x_2.size(); ++i)
    {
        std::cout << "x" << i + 1 << " = " << x_2[i][0] << std::endl;
    }
    return 0;
}
void print(matrix& x)
{
    int n = x.size();
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            cout << x[i][j] << ' ';
        }
        cout << '\n';
    }
}
vector<vector<double>> transpose(const vector<vector<double>>& matrix) {
    if (matrix.empty()) {
        return {};
    }
    
    // Определяем размер результирующей транспонированной матрицы
    size_t numRows = matrix.size();
    size_t numCols = matrix[0].size();
    
    // Создаем матрицу для результата с переставленными размерностями
    vector<vector<double>> transposed(numCols, vector<double>(numRows));
    
    // Выполняем транспонирование
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }
    
    return transposed;
}


void solve2()
{
    int n;
    cin >> n;
    vector<vector<string>> matr(n, vector<string>(n, ""));
    for(int i = 0;i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cin >> matr[i][j]; 
        }
    }
    vector<vector<string>> matr2(n, vector<string>(n, ""));
    for(int i = 0;i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cin >> matr2[i][j]; 
        }
    }
    vector<vector<string>> result(n, std::vector<string>(n, ""));

    // Выполняем умножение матриц
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            for (size_t k = 0; k < n; ++k) {
                result[i][j] += matr[i][k] + "*" + matr2[k][j] + "+";
            }
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            result[i][j] = "(" + result[i][j] + ")";
        }
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            cout << result[i][j] << '\t';
        }
        cout << '\n';
    }
    // for(int i = 1; i < n; i++)
    // {
    //     matrix temp = matr;
    //     for(int j = 1; j < i; j++)
    //     {
    //         temp = mult(temp, matr);
    //     }
    //     result = result + temp;
    // }
    // for(int i = 0;i < n; i++)
    // {
    //     for(int j = 0; j < n; j++)
    //     {
    //         result[i][j] = result[i][j] > 0;
    //     }
    // }
    //result = result & transpose(result);
    // int gol = 0;
    // for(int i = 1; i < 4; i++)
    // {
    //     result = mult(matr, result);
    // }
    // for(int i = 0; i < n; i++)
    // {
    //     for(int j = 0; j < n;j++)
    //     {
    //         if(i != j)
    //             gol += result[i][j];
    //     }
    // }

    // matrix a;

    // matrix b;
    // print(result);
    // cout << gol;

}
int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    fastboi

    SPEEDTEST_START
    //int r = solve();
    solve2();
    SPEEDTEST_STOP
    return 0;
}
