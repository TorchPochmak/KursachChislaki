#include <bits/stdc++.h>
using namespace std;

#define fastboi ios::sync_with_stdio(false); cin.tie(0);

void check(
    vector<double>& a,
    vector<double>& b,
    vector<double>& c)
{
    try
    {    
        int n = a.size() - 1;
        for(int i = 1; i <= n; i++)
        {
            //|b| >= |a| + |c|
            if(std::abs(b[i]) < std::abs(a[i]) + std::abs(c[i]))
                throw std::runtime_error("Error: Нет диагонального преобладания");
        }
        for(int i = 2; i <= n - 1; i++)
        {
            if(a[i] == 0 || b[i] == 0 || c[i] == 0)
                throw std::runtime_error("Error: Не трехдиагональная");
        }
        if(b[1] == 0 || c[1] == 0 || a[n] == 0 || b[n] == 0)
            throw std::runtime_error("Error: Не трехдиагональная");
    }
    catch(const std::exception& e)
    {
        throw;
    }
}
std::pair<vector<double>, vector<double>> getPairForward(
    vector<double>& a,
    vector<double>& b,
    vector<double>& c,
    vector<double>& d,
    vector<double>& y)
{
    try
    {        
        int n = a.size() - 1;
        auto p = vector<double>(n + 1, 0);
        auto q = vector<double>(n + 1, 0);

        if(b[1] == 0)
            throw std::runtime_error("Error: Деление на 0");

        p[1] = -c[1]/b[1];
        q[1] = d[1]/b[1];

        for(int i = 1; i <= n; i++)
        {
            y[i] = (b[i] + a[i] * p[i - 1]);
            if(y[i] == 0)
                throw std::runtime_error("Error: Деление на 0");
            p[i] = -c[i] / y[i];
            q[i] = (d[i] - a[i] * q[i - 1]) / y[i];
        }
        return std::make_pair(p, q);
    }
    catch(const std::exception& e)
    {
        throw;
    }
}
vector<double> getResultReverse(
    vector<double>& p,
    vector<double>& q
)
{
    int n = p.size() - 1;
    auto x = vector<double>(n + 1, 0);
    x[n] = q[n];
    for(int i = n - 1; i >= 1; i--)
    {
        x[i] = p[i] * x[i + 1] + q[i];
    }
    return x;
}
void check_result(
    vector<double>& a,
    vector<double>& b,
    vector<double>& c,
    vector<double>& x,
    vector<double>& d
)
{
    x.resize(x.size() + 1); //чтобы последний элемент проверить не парясь
    cout << "\n Проверка результата: \n";
    for(int i = 1; i <= a.size() - 1; i++)
    {
        cout << "Строка " << i << "\n";
        cout << "x = " << a[i] * x[i - 1]  + b[i] * x[i] + c[i] * x[i + 1] << '\n';
        cout << "d = " << d[i] << '\n';
    }
}

double det(
    vector<double>& a,
    vector<double>& b,
    vector<double>& c,
    vector<double>& y
)
{
    double d = 1;
    for(int i = 1; i <= a.size() - 1; i++)
        d *= y[i];
    return d;
}
void solve()
{
    int n;
    cin >> n;
    if(n < 3)
        throw std::runtime_error("Error: N < 3");
    //векторы диагоналей a, b, c - b - главная
    //нулевые не трогаем, строки i = 1..n
    auto a = vector<double>(n + 1, 0);
    auto b = vector<double>(n + 1, 0);
    auto c = vector<double>(n + 1, 0);

    auto d = vector<double>(n + 1, 0);
    auto y = vector<double>(n + 1, 0);

    cin >> b[1] >> c[1] >> d[1];
    for(int i = 2; i <= n - 1; i++)
        cin >> a[i] >> b[i] >> c[i] >> d[i];
    cin >> a[n] >> b[n] >> d[n];
    try
    {
        check(a, b, c);
        auto pq_pair = getPairForward(a, b, c, d, y);
        auto x = getResultReverse(pq_pair.first, pq_pair.second);
        cout << "Столбец результата x: \n";
        for(int i = 1; i <= n; i++)
            cout << x[i] << '\n';
        cout << "Детерминант: " << det(a, b, c, y);

        //Проверка
        check_result(a,b,c,x,d);
    }
    catch(const std::exception& e)
    {
        throw;
    }
}
int main()
{
   
    cout << "Выберите режим \n 1) Чтение из файла input.txt \n 2) Из терминала \n";
    int choose = 0;
    while(true)
    {
        cout << "Введите ваш выбор (1 или 2)\n";
        cin >> choose;
        if(choose < 1 || choose > 2)
            cout << "Попробуйте снова\n";
        else
        {
            if(choose == 1)
            {
                fastboi
                freopen("input.txt", "r", stdin);
            }
            else
                cout << "Введите сначала размерность матрицы, а затем коэффициенты в каждой строке \n";
            break;
        }
    }
    try
    {
        solve();
    }
    catch(const std::exception& e)
    {
        cout << e.what() << '\n';
    }
}