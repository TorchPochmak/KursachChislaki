
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


void print(vector<double>& x)
{
    for(int i = 0; i < x.size(); i++)
    {
        cout << x[i] << ' ';
    }
    cout << '\n';
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


std::vector<double> Gauss_solve(matrix A) {
    int n = A.size();

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; ++i) {
        // Поиск максимального элемента в колонке для избежания деления на ноль
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        // Перестановка строк
        std::swap(A[i], A[maxRow]);

        // Проверка на единичный элемент в диагонали после перестановки
        if (std::abs(A[i][i]) < 1e-9) {
            throw std::runtime_error("Матрица вырождена или система имеет бесконечно много решений");
        }

        // Приведение диагонального элемента к 1 и обнуление ниже этой строки
        for (int k = i + 1; k < n; ++k) {
            double coeff = A[k][i] / A[i][i];
            for (int j = i; j <= n; ++j) {
                A[k][j] -= coeff * A[i][j];
            }
        }
    }

    // Обратный ход
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = A[i][n] / A[i][i];
        for (int k = i - 1; k >= 0; --k) {
            A[k][n] -= A[k][i] * x[i];
        }
    }

    return x;
}

std::vector<double> find_c (const std::vector<double>& f,
                            const std::vector<double>& h)
{
    matrix c(3, vector<double>(4));

    c[0][0] = 2 * (h[0] + h[1]);
    c[0][1] = h[1];
    c[0][2] = 0;
    c[0][3] = 3 * ((f[2] - f[1]) / h[1] - (f[1] - f[0]) / h[0]);

    c[1][0] = h[1];
    c[1][1] = 2 * (h[1] + h[2]);
    c[1][2] = h[2];
    c[1][3] = 3 * ((f[3] - f[2]) / h[2] - (f[2] - f[1]) / h[1]);

    c[2][0] = 0;
    c[2][1] = h[2];
    c[2][2] = 2 * (h[2] + h[3]);
    c[2][3] = 3 * ((f[4] - f[3]) / h[3] - (f[3] - f[2]) / h[2]);

    std::vector<double> res(4);
    res[0] = 0;

    vector<double> c_i = Gauss_solve(c);
    for (int i = 0; i < 3; ++ i)
    {
        res[i + 1] = c_i[i];
    }

    return res;
}

matrix cubic_spline (const std::vector<double>& x,
                     const std::vector<double>& f)
{
    matrix table(4, vector<double>(4));

    std::vector<double> h(4);
    for (int i = 0; i < 4; ++i)
    {
        h[i] = x[i + 1] - x[i];
    }

    //c
    table[2] = find_c(f, h);

    //a
    for (int i = 0; i < 4; ++i)
    {
        table[0][i] = f[i];
    }

    //b
    for (int i = 0; i < 3; ++i)
    {
        table[1][i] = (f[i + 1] - f[i]) / h[i] - h[i] * (table[2][i + 1] + 2 * table[2][i]) / 3;
    }
    table[1][3] = (f[4] - f[3]) / h[3] - 2 * h[3] * table[2][3] / 3;

    //d
    for (int i = 0; i < 3; ++i)
    {
        table[3][i] = (table[2][i + 1] - table[2][i]) / (3 * h[i]);
    }
    table[3][3] = -table[2][3] / (3 * h[3]);

    return table;
}

double cubic_spline_result (matrix& table, double value,
                            const std::vector<double>& x)
{
    if (value < x[0] || value > x[x.size() - 1])
    {
        throw std::runtime_error("Cubic spline is not defined on this area");
    }

    int i = 0;
    while (value >= x[i])
    {
        ++i;
    }
    --i;

    return table[0][i] + table[1][i] * (value - x[i]) + table[2][i] * pow((value - x[i]), 2)
            + table[3][i] * pow((value - x[i]), 3);
}

void cubic_spline_polynomials (matrix& table, const std::vector<double>& x)
{
    for (int i = 0; i < x.size() - 1; ++i)
    {
        std::cout << x[i] << " - " << x[i + 1] << ": ";
        std::cout << table[0][i] << " + " << table[1][i] << " * (x - " << x[i] << ") + " << table[2][i]
        << " * (x - " << x[i] << ")^2" << " + " << table[3][i] << " * (x - " << x[i] << ")^3" << std::endl;
    }
}

void solve()
{
    // let's go kid
    std::vector<double> x(5);
    std::vector<double> f(5);
    double value;

    for (int i = 0; i < 5; ++i)
    {
        cin >> x[i];
    }

    for (int i = 0; i < 5; ++i)
    {
        cin >> f[i];
    }

    cin >> value;

    for (int i = 0; i < 5; ++ i)
    {
        std::cout << x[i] << "  " << f[i] << std::endl;
    }

    matrix res = cubic_spline(x, f);
    print(res);

    std::cout << "Result: " << cubic_spline_result(res, value, x) << std::endl;

    cubic_spline_polynomials(res, x);
}
int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    fastboi

    SPEEDTEST_START
    solve();
    SPEEDTEST_STOP
}