#include <cmath>
#include <functional>
#include <iomanip>
#include <algorithm>
#include <bits/stdc++.h>
#include <vector>
using namespace std;
double f (double x, double y, double z)
{
    return z;
}

double g (double x, double y, double z)
{
    //return ((x + 1) * z - y) / x;
    return (4*x*z-(4*x*x-2)*y);
}

double func (double x)
{
    return (1+x)*exp(x*x);
}

vector<vector<double>> Euler_method (double x_1, double x_2, double y, double z, double h)
{
    int n = (x_2 - x_1) / h;
    vector<vector<double>> res(n+1, vector<double>(4));

    double x = x_1;
    res[0][0] = x;
    res[0][1] = y;

    res[0][2] = func(x);
    res[0][3] = std::abs(res[0][2] - res[0][1]);
    
    double next_x;
    double next_y;

    for (int i = 0; i < n; ++i)
    {
        next_x = x + h;

        next_y = y + h * f(x + h / 2,
                           y + h / 2 * f(x, y, z),
                           z + h / 2 * g(x, y, z));

        z += h * g(x + h / 2,
                   y + h / 2 * f(x, y, z),
                   z + h / 2 * g(x, y, z));

        x = next_x;
        y = next_y;

        res[i + 1][0] = x;
        res[i + 1][1] = y;

        res[i + 1][2] = func(x);
        res[i + 1][3] = std::abs(res[i + 1][2] - res[i + 1][1]);
    }

    return res;
}

void K_and_L (double h, double x, double y, double z,
              std::vector<double> &K, std::vector<double> &L)
{
    K[0] = h * f(x, y, z);
    L[0] = h * g(x, y, z);

    for (int i = 1; i <= 2; ++i)
    {
        K[i] = h * f(x + 0.5 * h, y + 0.5 * K[i - 1], z + 0.5 * L[i - 1]);
        L[i] = h * g(x + 0.5 * h, y + 0.5 * K[i - 1], z + 0.5 * L[i - 1]);
    }

    K[3] = h * f(x + h, y + K[2], z + L[2]);
    L[3] = h * g(x + h, y + K[2], z + L[2]);
}

vector<vector<double>> Runge_Kutta_method (double x_1, double x_2, double y, double z, double h)
{
    int n = (x_2 - x_1) / h;
    vector<vector<double>> res(n+1, vector<double>(5));

    std::vector<double> K(4);
    std::vector<double> L(4);

    double x = x_1;
    res[0][0] = x;
    res[0][1] = y;
    res[0][2] = z;

    res[0][3] = func(x);
    res[0][4] = std::abs(res[0][3] - res[0][1]);

    for (int i = 0; i < n; ++i)
    {
        K_and_L(h, x, y, z, K, L);

        x += h;
        y += (K[0] + 2 * K[1] + 2 * K[2] + K[3]) / 6;
        z += (L[0] + 2 * L[1] + 2 * L[2] + L[3]) / 6;

        res[i + 1][0] = x;
        res[i + 1][1] = y;
        res[i + 1][2] = z;

        res[i + 1][3] = func(x);
        res[i + 1][4] = std::abs(res[i + 1][3] - res[i + 1][1]);
    }

    return res;
}

vector<vector<double>> Adams_method (double x_1, double x_2, double y, double z, double h, vector<vector<double>> &M)
{
    int n = (x_2 - x_1) / h;
    vector<vector<double>> res(n+1, vector<double>(5));

    for (int i = 0; i < 4; ++i)
    {
        res[i][0] = M[i][0];
        res[i][1] = M[i][1];
        res[i][2] = M[i][2];

        res[i][3] = M[i][3];
        res[i][4] = M[i][4];
    }

    double x = M[3][0];
    y = M[3][1];
    double p_y, p_z;

    for (int i = 4; i <= n; ++i)
    {
        x += h;
        p_y = y + h / 24 * (55 * f(res[i - 1][0], res[i - 1][1], res[i - 1][2]) -
                       59 * f(res[i - 2][0], res[i - 2][1], res[i - 2][2]) +
                       37 * f(res[i - 3][0], res[i - 3][1], res[i - 3][2]) -
                        9 * f(res[i - 4][0], res[i - 4][1], res[i - 4][2]));

        p_z = z + h / 24 * (55 * g(res[i - 1][0], res[i - 1][1], res[i - 1][2]) -
                       59 * g(res[i - 2][0], res[i - 2][1], res[i - 2][2]) +
                       37 * g(res[i - 3][0], res[i - 3][1], res[i - 3][2]) -
                        9 * g(res[i - 4][0], res[i - 4][1], res[i - 4][2]));

        y += h / 24 * ( 9 * f(x, p_y, p_z) +
                       19 * f(res[i - 1][0], res[i - 1][1], res[i - 1][2]) -
                        5 * f(res[i - 2][0], res[i - 2][1], res[i - 2][2]) +
                            f(res[i - 3][0], res[i - 3][1], res[i - 3][2]));

        z += h / 24 * ( 9 * g(x, p_y, p_z) +
                       19 * g(res[i - 1][0], res[i - 1][1], res[i - 1][2]) -
                        5 * g(res[i - 2][0], res[i - 2][1], res[i - 2][2]) +
                            g(res[i - 3][0], res[i - 3][1], res[i - 3][2]));


        res[i][0] = x;
        res[i][1] = y;
        res[i][2] = z;

        res[i][3] = func(x);
        res[i][4] = std::abs(res[i][3] - res[i][1]);
    }

    return res;
}


void print(vector<vector<double>>& x)
{
    int n = x.size();
    int m = x[0].size();
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < m; ++j)
        {
            cout << x[i][j] << ' ';
        }
        cout << '\n';
    }
}

void solve()
{
    // let's go kid
    double x_1 = 0;
    double x_2 = 1;
    double y = 1;
    double z = 1;
    double h = 0.1;

    vector<vector<double>> res, res_a;

    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    std::cout << "Euler method" << std::endl;
    res = Euler_method(x_1, x_2, y, z, h);
    print(res);

    //x y z r p

    std::cout << "Runge-Kutta method" << std::endl;
    res = Runge_Kutta_method(x_1, x_2, y, z, h);
    print(res);

    std::cout << "Adams method" << std::endl;
    res_a = Adams_method(x_1, x_2, y, z, h, res);
    print(res_a);

    cout << "Answer expected:\n";
    for(double i = 0; i < 1; i += h)
    {
        cout << i << " " << func(i);
        cout << '\n';
    }
}
int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    solve();
}