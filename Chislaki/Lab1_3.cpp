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
double EPS = 1e-6;

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

bool diag_dominant(matrix& matr)
{
    size_t n = matr.size();
    for (size_t i = 0; i < n; ++i) {
        double sum_d = std::abs(matr[i][i]);
        double sum = 0.0;

        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                sum += std::abs(matr[i][j]);
            }
        }

        if (sum_d <= sum) {
            return false;
        }
    }
    return true;
}

// норма для точности |cur - prev|
double norm(vector<double>& x, vector<double>& prev_x)
{
    double norm = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        norm = std::max(norm, std::abs(x[i] - prev_x[i]));
    }
    return norm;
}

vector<double> solve_simple_iter(matrix& matr, vector<double>& b, int n, int& iter_total)
{
    iter_total = 0;

    vector<double> cur(n, 0);
    vector<double> prev(n, 0);

    do
    {
        prev = cur;
        for(int i = 0; i < n; ++i)
        {
            double sm = 0;
            for(int j = 0; j < n; j++)
            {
                if(i != j)
                {
                    sm += matr[i][j] * prev[j];
                }
            }
            if(abs(matr[i][i]) < EPS)
            {
                throw logic_error("ERROR: деление на ноль");
            } 
            cur[i] = (b[i] - sm) / matr[i][i];
        }
        iter_total++;
        cout << fixed << setprecision(6) << "Итоговая норма: " << norm(cur, prev) << '\n';
    } while (norm(cur, prev) >= EPS);
    return cur;
}

vector<double> solve_seidel(matrix& matr, vector<double>& b, int n, int& iter_total)
{
    vector<double> cur(n, 0);
    vector<double> prev(n, 0);

    do
    {
        prev = cur;
        for(int i = 0; i < n; ++i)
        {
            double sm = 0;
            for(int j = 0; j < i; ++j)
            {
                sm += matr[i][j] * cur[j];
            }
            for(int j = i + 1; j < n; j++)
            {
                sm += matr[i][j] * prev[j];
            }
            if(abs(matr[i][i]) < EPS)
            {
                throw logic_error("ERROR: деление на ноль");
            } 
            cur[i] = (b[i] - sm) / matr[i][i];
        }
        iter_total++;
        //cout << fixed << setprecision(6) << "Итоговая норма: " << norm(cur, prev) << '\n';
    } while (norm(cur, prev) >= EPS);
    return cur;   
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

vector<double> to_vector(matrix a)
{
    vector<double> res = {};
    if(a.size() == 1)
    {
        for(int i = 0; i < a[0].size(); i++)
            res.push_back(a[0][i]);
    }
    else if(a[0].size() == 1)
    {
        for(int i = 0; i < a.size(); i++)
            res.push_back(a[i][0]);
    }
    return res;
}

void solve()
{
    // let's go kid
    int n;
    cin >> n;
    matrix matr(n, vector<double>(n, 0));
    vector<double> b(n, 0);
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < n; ++j)
        {
            cin >> matr[i][j];
        }
    }
    for(int i = 0; i < n; ++i)
    {
        cin >> b[i];
    }
    if(!diag_dominant(matr))
    {
        cout << "Матрица не диаг. доминирующая\n";
        return;
    }

    int iter_total_simple = 0;
    int iter_total_seidel = 0;

    vector<double> result_simple(n, 0);
    vector<double> result_seidel(n, 0);

    try
    {
        result_simple = solve_simple_iter(matr,b, n, iter_total_simple);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << '\n';
        return;
    }

    try
    {
        result_seidel = solve_seidel(matr, b, n, iter_total_seidel);
    }
    catch(const std::exception& e)
    {
        std::cout << e.what() << '\n';
        return;
    }
    
    std::cout << "> Решение методом простых итераций:";
    print(result_simple);
    std::cout << "> Количество итераций: " << iter_total_simple << '\n';

    std::cout << "\n> Решение методом Зейделя:" << '\n';
    print(result_seidel);
    std::cout << "> Количество итераций: " << iter_total_seidel << '\n';

    vector<double> check_input_iters = vector<double>(n,0);
    vector<double> check_input_seidel = vector<double>(n,0);
    
    try
    {
       check_input_iters = to_vector(mult_matr_vec(matr, result_simple));
    }
    catch(std::exception& e)
    {
        cout << e.what() << '\n';
    }

    try
    {
       check_input_seidel = to_vector(mult_matr_vec(matr, result_seidel));
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    cout << "Проверка метода простых итераций" << '\n';
    print(check_input_iters);
    cout << "Проврека метода Зейделя" << '\n';
    print(check_input_seidel);

    return;
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