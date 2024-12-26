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

using matrix =vector<vector<double>>;

double EPS = 1e-6;
int MAX_ITER = 100000;

vector<double> matVecMultiply(const matrix& matrix, const vector<double>& x) {
    int n = matrix.size();
    vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i] += matrix[i][j] * x[j];
        }
    }
    return result;
}

void print(vector<double>& x)
{
    for(int i = 0; i < x.size(); i++)
    {
        cout << x[i] << ' ';
    }
    cout << '\n';
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

void norm(vector<double>& v) {
    double norm = 0.0;
    for (double elem : v) {
        norm += elem * elem;
    }
    norm = sqrt(norm);
    for (double& elem : v) {
        elem /= norm;
    }
}

// Функция поворота для обнуления элемента matrix[p][q]
void rotate(matrix& matr, matrix& V, int p, int q) {
    // Если уже ноль, ничего не делаем
    if (matr[p][q] == 0) {
        return;
    }
    int n = matr.size(); // Размер матрицы
    // Вычисляем угол theta для поворота
    double theta = 0.5 * atan2(2 * matr[p][q], matr[q][q] - matr[p][p]);
    // Коррекция угла для близких диагональных значений
    if (abs(matr[p][p] - matr[q][q]) < EPS) {
        theta = cos(-1) / 4; // Устанавливаем на 45 градусов
    }
    // Вычисляем косинус и синус от угла theta
    double cosTheta = cos(theta);
    double sinTheta = sin(theta);

    // Сохраняем текущие значения диагональных элементов и элемента вне диагонали
    double app = matr[p][p], aqq = matr[q][q], apq = matr[p][q];
    
    // Обновляем элементы matrix[p][p] и matrix[q][q] после вращения
    matr[p][p] = cosTheta * cosTheta * app + sinTheta * sinTheta * aqq - 2 * sinTheta * cosTheta * apq;
    matr[q][q] = sinTheta * sinTheta * app + cosTheta * cosTheta * aqq + 2 * sinTheta * cosTheta * apq;
    
    // Обнуляем элемента вне диагонали matrix[p][q] и matrix[q][p]
    matr[p][q] = matr[q][p] = 0;

    // Обновляем другие элементы, затронутые вращением
    for (int i = 0; i < n; i++) {
        if (i != p && i != q) {
            double aip = matr[i][p], aiq = matr[i][q];
            matr[i][p] = matr[p][i] = cosTheta * aip - sinTheta * aiq;
            matr[i][q] = matr[q][i] = sinTheta * aip + cosTheta * aiq;
        }
    }

    // Обновляем матрицу собственных векторов V
    for (int i = 0; i < n; i++) {
        double vip = V[i][p], viq = V[i][q];
        V[i][p] = cosTheta * vip - sinTheta * viq;
        V[i][q] = sinTheta * vip + cosTheta * viq;
    }
}
// Функция для нахождения собственных значений и векторов матрицы методом Якоби
void jacobiEigenvalMethod(matrix& matr, vector<double>& eigenvals, matrix& eigenvecs) 
{
    int n = matr.size(); // Размер матрицы

    //единичкая
    eigenvecs.assign(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; i++) 
    {
        eigenvecs[i][i] = 1.0;
    }

    // Основной цикл итераций метода Якоби
    for (int k = 0; k < MAX_ITER; k++) 
    {
        int p = 0, q = 1;
        double max = 0;
        
        // Находим наибольший по модулю элемент вне диагонали
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) 
            {
                if (abs(matr[i][j]) > max) 
                {
                    max = abs(matr[i][j]);
                    p = i;
                    q = j;
                }
            }
        }
        
        // Если все элементы вне диагонали меньше заданного порога, завершаем
        if (max < EPS) 
        {
            break; 
        }
        
        // Выполняем ротацию для обнуления matrix[p][q]
        rotate(matr, eigenvecs, p, q);
    }

    // Извлекаем собственные значения из диагональных элементов матрицы
    for (int i = 0; i < n; i++) 
    {
        eigenvals[i] = matr[i][i];
    }
}

// Функция для нахождения наибольшего собственного значения и вектора методом Лейбница (степенное возведение)
pair<double, vector<double>> powerMethod(matrix& matr, int maxIterations = MAX_ITER, double eps = EPS) 
{
    int n = matr.size();
    vector<double> x(n, 1.0);  
    double eigenval = 0.0;

    for (int k = 0; k < maxIterations; k++) {
        vector<double> y = matVecMultiply(matr, x);
        
        double tmpEigenval = y[0] / x[0];
        
        norm(y);
        
        if (fabs(tmpEigenval - eigenval) < eps) {
            return {tmpEigenval, y};
        }

        x = y;

        eigenval = tmpEigenval;
    }

    return {eigenval, x};
}

void solve()
{
    // let's go kid
    // let's go kid
    int n;
    cin >> n;
    matrix matr(n, vector<double>(n));
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++)
        {
            cin >> matr[i][j];
        }
    }
    matrix matr2 = matr;
    vector<double> eigenvals(n);
    matrix eigenvecs;

    jacobiEigenvalMethod(matr, eigenvals, eigenvecs);

    cout << "> Собственные значения (Якоби):\n";
    for (double elem : eigenvals) 
    {
        cout << elem << ' ';
    }
    cout << "\n\n> Собственные векторы (Якоби):\n";
    print(eigenvecs);


    auto [eigenval, eigenvector] = powerMethod(matr2);

    cout << "\n> Собственное значение (Степенной метод): " << eigenval << '\n';
    cout << "\n> Собственный вектор (Степенной метод): ";
    for (auto elem : eigenvector) {
        cout << elem << ' ';
    }
    cout << '\n';
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