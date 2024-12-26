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

matrix identityMatrix(int n)
{
    vector<vector<double>> result = vector<vector<double>>(n, vector<double>(n, 0));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            result[i][j] = 1;
        }
    }
    return result;
}

void print(const matrix& vec) 
{
    for (int i = 0; i < vec.size(); i++) 
    {
        for (int j = 0; j < vec[i].size(); j++) 
        {
            cout << setw(6) << fixed << setprecision(2) << vec[i][j] << ' ';
        }
        cout << '\n';
    }
}

matrix addMatrixes(matrix a, matrix b, int n, int m)
{
    matrix result = matrix(n, vector<double>(m, 0));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            result[i][j] = a[i][j] + b[i][j];
        }
    }
}

matrix multMatrixes(matrix a, int n1, int m1, matrix b, int n2, int m2)
{

}

void solve()
{
    // let's go kid

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