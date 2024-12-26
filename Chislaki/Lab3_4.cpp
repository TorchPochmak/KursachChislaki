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


size_t find_x (const std::vector<double>& x, double value)
{
    double EPS = 1E-10;

    for (size_t i = 0; i < x.size(); ++i)
    {
        if (std::abs(x[i] - value) < EPS)
        {
            return i;
        }
    }

    return -1;
}

double first_order_of_accuracy_derivative (double x_1, double x_2, double y_1, double y_2)
{
    return (y_2 - y_1) / (x_2 - x_1);
}

double derivative (const std::vector<double>& x,
                   const std::vector<double>& y,
                   double value)
{
    size_t ind = find_x(x, value);

    if (ind == -1)
    {
        throw std::runtime_error("The function is not defined at this point");
    }

    if (ind == 0)
    {
        return first_order_of_accuracy_derivative(x[0], x[1], y[0], y[1]);
    }
    if (ind == x.size() - 1)
    {
        return first_order_of_accuracy_derivative(x[ind - 1], x[ind], y[ind - 1], y[ind]);
    }

    double left = first_order_of_accuracy_derivative(x[ind - 1], x[ind], y[ind - 1], y[ind]);
    double right = first_order_of_accuracy_derivative(x[ind], x[ind + 1], y[ind], y[ind + 1]);

    return left + (right - left) / (x[ind + 1] - x[ind - 1]) * (2 * value - x[ind - 1] - x[ind]);
}

double second_derivative (const std::vector<double>& x,
                          const std::vector<double>& y,
                          double value)
{
    size_t ind = find_x(x, value);

    if (ind == -1 || ind == 0 || ind == x.size() - 1)
    {
        throw std::runtime_error("It is impossible to find the second derivative");
    }

    double left = first_order_of_accuracy_derivative(x[ind - 1], x[ind], y[ind - 1], y[ind]);
    double right = first_order_of_accuracy_derivative(x[ind], x[ind + 1], y[ind], y[ind + 1]);

    return 2 * (right - left) / (x[ind + 1] - x[ind - 1]);
}


void solve()
{
    // let's go kid
    std::vector<double> x(5);
    for(int i = 0; i < 5; i++)
    {
        cin >> x[i];
    }
    std::vector<double> y(5);
    for(int i = 0; i < 5; i++)
    {
        cin >> y[i];
    }
    double value;
    cin >> value;

    
    std::cout << "The first derivative: " << derivative(x, y, value) << std::endl;
    std::cout << "The second derivative: " << second_derivative(x, y, value) << std::endl;
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
