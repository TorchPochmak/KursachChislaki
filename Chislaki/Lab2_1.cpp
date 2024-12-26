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

//2^x+x^2-2=0
double function_main (double x)
{
    return pow(2, x) + x * x - 2;
}

double derivative_function (double x)
{
    return log(2) * pow(2, x) + 2 * x;
}

double derivative_2_function (double x)
{
    return pow(log(2), 2) * pow(2, x) + 2;
}

double equivalent_function (double x)
{
    //return (pow(4, x) - 2) / 5;
    return sqrt(2-pow(2,x));
}

double derivative_equivalent_function (double x)
{
    //return log(4) * pow(4, x) / 5;
    return -log(2)*pow(2,x-1) / (sqrt(2-pow(2,x)));
}

double dichotomy_method (const function<double(double)>& func, double a, double b, double EPS)
{
    double x = (a + b) / 2;
    EPS *= 2;

    while (abs(func(x)) > EPS)
    {
        double f_x = func(x);

        if (func(a) * f_x < 0)
        {
            b = x;
            x = (a + b) / 2;
        }
        else if (func(b) * f_x < 0)
        {
            a = x;
            x = (a + b) / 2;
        }
        else
        {
            throw runtime_error("There is no root in this range");
        }
    }

    return x;
}

double Newton_method (const function<double(double)>& func, 
                      const function<double(double)>& d_func,
                      const function<double(double)>& d2_func,
                      double a, double b, double EPS)
{
    double x = a;
    double prev_x;

    if (func(a) * func(b) >= 0)
    {
        throw runtime_error("Newton method does not converge");
    }

    while (func(x) * d2_func(x) <= 0)
    {
        x = (x + b) / 2;
    }

    do
    {
        prev_x = x;

        x -= func(x) / d_func(x);

    } while (abs(x - prev_x) > EPS);

    return x;
}

double secant_method (const function<double(double)>& func,
                      const function<double(double)>& d2_func,
                      double a, double b, double EPS)
{
    double x = a;
    double prev_x;

    while (func(x) * d2_func(x) <= 0)
    {
        x = (x + b) / 2;
    }

    prev_x = x;
    x = (prev_x + b) / 2;

    while (abs(x - prev_x) > EPS)
    {
        double tmp = x;
        x -= func(x) * (x - prev_x) / (func(x) - func(prev_x));

        prev_x = tmp;
    }

    return x;
}

double simple_iteration_method (const function<double(double)>& func,
                                const function<double(double)>& d_func,
                                double a, double b, double EPS)
{
    double q = d_func(a);

    if (func(a) < a || func(b) > b || q >= 1)
    {
        throw runtime_error("Simple iteration method does not converge");
    }

    double x = (a + b) / 2;
    double x_prev;

    do
    {
        x_prev = x;

                x = func(x);

    } while (abs(x - x_prev) > EPS);

    return x;
}

int solve()
{
    double EPS = 1E-3;
    cout << fixed;
    cout << setprecision(8);

    double x_1;

    try
    {
        x_1 = dichotomy_method(function_main, -1, 3, EPS);
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what() <<'\n';

        return 1;
    }

    cout << "Dichotomy method: " << x_1 <<'\n';

    double x_2;

    try
    {
        x_2 = Newton_method(function_main, derivative_function, derivative_2_function, -1, 3, EPS);
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what() <<'\n';

        return 1;
    }

    cout << "Newton method: " << x_2 <<'\n';

    double x_3;

    try
    {
        x_3 = secant_method(function_main, derivative_2_function,-1, 3, EPS);
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what() <<'\n';

        return 1;
    }

    cout << "Secant method: " << x_3 <<'\n';

    double x_4 = 0;

    try
    {
        x_4 = simple_iteration_method(equivalent_function, derivative_equivalent_function, 0.5, 0.71, EPS);
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what() <<'\n';

        return 1;
    }

    cout << "Simple iteration method: " << x_4 << '\n';
    return 0;
    // let's go kid
}
int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    fastboi

    SPEEDTEST_START
    int res = solve();
    SPEEDTEST_STOP
    return res;
}