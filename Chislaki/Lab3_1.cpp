#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include <cmath>

double function (double x)
{
    return sqrt(x);
}

std::vector<double> Lagrange_polynomial (const std::function<double(double)>& func,
                                         const std::vector<double> &x)
{
    std::vector<double> L(x.size());

    for (int i = 0; i < x.size(); ++i)
    {
        L[i] = func(x[i]);
        for (int j = 0; j < x.size(); ++j)
        {
            if (i != j)
            {               
                L[i] /= x[i] - x[j];
            }
        }
    }

    return L;
}

double Lagrange_polynomial_result (const std::vector<double> &L,
                                   const std::vector<double> &x,
                                   double value)
{
    double res = 0;

    for (int i = 0; i < L.size(); ++i)
    {
        double tmp = L[i];

        for (int j = 0; j < x.size(); ++j)
        {
            if (i != j)
            {
                tmp *= value - x[j];
            }
        }

        res += tmp;
    }

    return res;
}

void print_Lagrange_polynomial (const std::vector<double> &L,
                                const std::vector<double> &x)
{
    for (int i = 0; i < L.size(); ++i)
    {
        if (i != 0 && L[i] >= 0) std::cout << "+";
        std::cout << L[i];

        for (int j = 0; j < x.size(); ++j)
        {
            if (i != j)
            {
                std::cout << "(x";
                if (x[j] >= 0) std::cout << "+";
                std::cout << x[j] << ")";
            }
        }
    }

    std::cout << std::endl;
}

std::vector<double> separated_differences (const std::function<double(double)>& func,
                                           const std::vector<double> &x)
{
    int size = (1 + x.size()) * x.size() / 2;
    std::vector<double> f(size);

    for (int i = 0; i < x.size(); ++i)
    {
        f[i] = func(x[i]);
    }

    int k = x.size();
    int p = 0;
    int e = x.size();
    for (int j = 1; j < x.size(); ++j)
    {
        for (int i = 0; i < x.size() - j; ++i)
        {
            f[e] = (f[p + i] - f[p + i + 1]) / (x[i] - x[i + j]);
            ++e;
        }

        p += k;
        --k;
    }

    return f;
}

std::vector<double> Newton_polynomial (const std::function<double(double)>& func,
                                       const std::vector<double> &x)
{
    std::vector<double> P(x.size());

    std::vector<double> f = separated_differences(func, x);

    int p = 0;
    int k = x.size();
    for (int i = 0; i < P.size(); ++i)
    {
        P[i] = f[p];

        p += k;
        --k;
    }

    return P;
}

double Newton_polynomial_result (const std::vector<double> &P,
                                 const std::vector<double> &x,
                                 double value)
{
    double res = 0;

    for (int i = 0; i < P.size(); ++i)
    {
        double tmp = P[i];

        for (int j = 0; j < i; ++j)
        {
            tmp *= value - x[j];
        }

        res += tmp;
    }

    return res;
}

void print_Newton_polynomial (const std::vector<double> &P,
                              const std::vector<double> &x)
{
    for (int i = 0; i < P.size(); ++i)
    {
        if (i != 0 && P[i] >= 0) std::cout << "+";
        std::cout << P[i];

        for (int j = 0; j < i; ++j)
        {
            std::cout << "(x";
            if (x[j] >= 0) std::cout << "+";
            std::cout << x[j] << ")";
        }
    }

    std::cout << std::endl;
}


void solve()
{
    // let's go kid
    std::vector<double> a = {0, 1.7, 3.4, 5.1};
    std::vector<double> b = {0,1.7,4.0,5.1};
    double X = 3;
    double f = function(X);

    std::cout << "Lagrange polynomial" << std::endl;
    std::cout << "a)" << std::endl;

    std::vector<double> L = Lagrange_polynomial(function, a);
    print_Lagrange_polynomial(L, a);

    double La = Lagrange_polynomial_result(L, a, X);
    std::cout << "Lagrange polynomial result: " << La << std::endl;
    std::cout << "Function result: " << f << std::endl;
    std::cout << "Error rate: " << std::abs(La - f) << std::endl;

    std::cout << "b)" << std::endl;

    std::vector<double> L_2 = Lagrange_polynomial(function, b);
    print_Lagrange_polynomial(L_2, b);

    double Lb = Lagrange_polynomial_result(L_2, b, X);
    std::cout << "Lagrange polynomial result: " << Lb << std::endl;
    std::cout << "Function result: " << f << std::endl;
    std::cout << "Error rate: " << std::abs(Lb - f) << std::endl;

    std::cout << "Newton polynomial" << std::endl;
    std::cout << "a)" << std::endl;

    std::vector<double> P = Newton_polynomial(function, a);
    print_Newton_polynomial(P, a);

    double Pa = Newton_polynomial_result(P, a, X);
    std::cout << "Newton polynomial result: " << Pa << std::endl;
    std::cout << "Function result: " << f << std::endl;
    std::cout << "Error rate: " << std::abs(Pa - f) << std::endl;

    std::cout << "b)" << std::endl;

    std::vector<double> P_2 = Newton_polynomial(function, b);
    print_Newton_polynomial(P_2, b);

    double Pb = Newton_polynomial_result(P_2, b, X);
    std::cout << "Newton polynomial result: " << Pb << std::endl;
    std::cout << "Function result: " << f << std::endl;
    std::cout << "Error rate: " << std::abs(Pb - f) << std::endl;
}
int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);


    solve();

}