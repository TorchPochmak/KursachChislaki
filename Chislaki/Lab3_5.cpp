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

double funct (const double x)
{
    return 1 / (3 * x * x + 4 * x + 2);
    //return x / std::pow(3*x + 4, 2);
}

double rectangle_method (const std::function<double(double)>& func,
                         double x_0, double x_1, double h)
{
    double F = 0;
    size_t N = std::round((x_1 - x_0) / h);
    double x = x_0;
    double next_x = x + h;

    for (int i = 0; i < N; ++i)
    {
        F += func((x + next_x) / 2);
        x = next_x;
        next_x += h;
    }

    return F * h;
}

double trapezoid_method (const std::function<double(double)>& func,
                         double x_0, double x_1, double h)
{
    size_t N = std::round((x_1 - x_0) / h);
    double x = x_0;

    double F = func(x) / 2;
    x += h;

    for (int i = 0; i < N - 1; ++i)
    {
        F += func(x);
        x += h;
    }

    F += func(x) / 2;

    return F * h;
}

double Simpson_method (const std::function<double(double)>& func,
                       double x_0, double x_1, double h)
{
    size_t N = std::round((x_1 - x_0) / h);
    double x = x_0;

    double F = func(x);
    x += h;

    for (int i = 0; i < N - 1; ++i)
    {
        if (i % 2 == 0)
        {
            F += 4 * func(x);
        }
        else
        {
            F += 2 * func(x);
        }

        x += h;
    }

    F += func(x);

    return F * h / 3;
}

double Runge_Romberg_Richardson_method (double F_h, double F_kh, double h_2, double h_1, int p)
{
    double k = h_1 / h_2;

    return F_h + (F_h - F_kh) / (pow(k, p) - 1);
}

void solve()
{
    // let's go kid
    std::cout << std::fixed;
    std::cout << std::setprecision(8);

    double exact_value = 1.88742;

    double x_0 = -2;
    double x_1 = 2;
    double h_1 = 1;
    double h_2 = 0.5;

//    double x_0 = -1;
//    double x_1 = 1;
//    double h_1 = 0.5;
//    double h_2 = 0.25;

    double rectangle = rectangle_method(funct, x_0, x_1, h_1);
    double trapezoid = trapezoid_method(funct, x_0, x_1, h_1);
    double Simpson = Simpson_method(funct, x_0, x_1, h_1);

    std::cout << "The rectangle method for h_1: " << rectangle << std::endl;
    std::cout << "The trapezoid method for h_1: " << trapezoid << std::endl;
    std::cout << "The Simpson method for h_1: " << Simpson << std::endl;

    std::cout << std::endl;

    double rectangle_2 = rectangle_method(funct, x_0, x_1, h_2);
    double trapezoid_2 = trapezoid_method(funct, x_0, x_1, h_2);
    double Simpson_2 = Simpson_method(funct, x_0, x_1, h_2);

    std::cout << "The rectangle method for h_2: " << rectangle_2 << std::endl;
    std::cout << "The trapezoid method for h_2: " << trapezoid_2 << std::endl;
    std::cout << "The Simpson method for h_2: " << Simpson_2 << std::endl;

    std::cout << std::endl;

    std::cout << "The Runge-Romberg-Richardson method" << std::endl;

    double RRR_r = Runge_Romberg_Richardson_method(rectangle_2, rectangle, h_2, h_1, 2);
    double RRR_t = Runge_Romberg_Richardson_method(trapezoid_2, trapezoid, h_2, h_1, 2);
    double RRR_s = Runge_Romberg_Richardson_method(Simpson_2, Simpson, h_2, h_1, 4);

    std::cout << "The rectangle method: " << RRR_r << std::endl;
    std::cout << "Error rate: " << std::abs(exact_value - RRR_r) << std::endl;
    std::cout << "The trapezoid method: " << RRR_t << std::endl;
    std::cout << "Error rate: " << std::abs(exact_value - RRR_t) << std::endl;
    std::cout << "The Simpson method: " << RRR_s << std::endl;
    std::cout << "Error rate: " << std::abs(exact_value - RRR_s) << std::endl;
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