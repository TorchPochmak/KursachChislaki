#include <bits/stdc++.h>
using namespace std;

void print(const vector<double>& vec) {
    for (double val : vec) 
    {
        cout << setw(6) << fixed << setprecision(2) << val << ' ';
    }
    cout << '\n';
}

void print(const vector<vector<double>>& vec) 
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

tuple<double, vector<double>> gauss_calc_lineq(vector<vector<double>> matrix, vector<double> coef) 
{    
    int n = matrix.size();
    double det = 1;

    //
    for (int i = 0; i < n; ++i) 
    {
        int swap_row = i;
        for (int k = i + 1; k < n; ++k) 
        {
            if (std::abs(matrix[k][i]) > std::abs(matrix[swap_row][i])) 
            {
                swap_row = k;
            }
        }

        if (swap_row != i) 
        {
            det *= -1; // ИЗМЕНИТЬ ЗНАК ПРИ ПЕРЕСТАНОВКЕ
            swap(matrix[i], matrix[swap_row]);
            swap(coef[i], coef[swap_row]);
        }

        det *= matrix[i][i];

        for (int k = i + 1; k < n; ++k) 
        {
            if (matrix[i][i] == 0) 
            {
                throw logic_error("NULL DIVISION\n");
            }
            double factor = matrix[k][i] / matrix[i][i];
            coef[k] -= factor * coef[i];
            for (int j = i; j < n; ++j) 
            {
                matrix[k][j] -= factor * matrix[i][j];
            }
        }
        std::cout << "Преобразование к верхнетреугольной, действие " << i + 1 << ": " << '\n';
        print(matrix);
        std::cout << "\n\n";
    }
    
    std::vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) 
    {
        x[i] = coef[i] / matrix[i][i];
        for (int k = i - 1; k >= 0; --k) 
        {
            coef[k] -= matrix[k][i] * x[i];
        }
    }

    return {det, x};
}

vector<vector<double>> gauss_inverse(vector<vector<double>> matrix) 
{
    int n = matrix.size();

    // расширенная
    vector<vector<double>> augmented_matrix(n, vector<double>(2 * n, 0));

    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            augmented_matrix[i][j] = matrix[i][j];
        }
        augmented_matrix[i][n + i] = 1.0;
    }
    cout << "Расширенная матрица, начало " << '\n';
    print(augmented_matrix);
    cout << '\n';

    for (int i = 0; i < n; ++i) 
    {
        int swap_row = i;
        for (int k = i + 1; k < n; ++k) 
        {
            if (std::abs(augmented_matrix[k][i]) > std::abs(augmented_matrix[swap_row][i])) 
            {
                swap_row = k;
            }
        }

        std::swap(augmented_matrix[i], augmented_matrix[swap_row]);

        if (augmented_matrix[i][i] == 0) 
        {
            throw std::logic_error("NULL DIVISION\n");
        }

        double divisor = augmented_matrix[i][i];
        for (int j = 0; j < 2 * n; ++j) 
        {
            augmented_matrix[i][j] /= divisor;
        }

        // вычитаем текущую строку из следующих строк
        for (int k = i + 1; k < n; ++k) 
        {
            double factor = augmented_matrix[k][i];
            for (int j = 0; j < 2 * n; ++j) 
            {
                augmented_matrix[k][j] -= factor * augmented_matrix[i][j];
            }
        }
        std::cout << "Расширенная матрица, прямой ход, действие " << i + 1 << '\n';
        print(augmented_matrix);
        std::cout << '\n';
    }

    // обратный ход
    for (int i = n - 1; i >= 0; --i) 
    {
        for (int k = i - 1; k >= 0; --k) 
        {
            double factor = augmented_matrix[k][i];
            for (int j = 0; j < 2 * n; ++j) 
            {
                augmented_matrix[k][j] -= factor * augmented_matrix[i][j];
            }
        }
        std::cout << "Расширенная матрица, обратный ход, действие " << (n - i) << '\n';
        print(augmented_matrix);
        std::cout << '\n';
    }

    std::vector<std::vector<double>> inverse_matrix(n, std::vector<double>(n));
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            inverse_matrix[i][j] = augmented_matrix[i][n + j];
        }
    }

    return inverse_matrix;
}

vector<vector<double>> mult_matrix(const vector<vector<double>>& a, const vector<vector<double>>& b) 
{
    int a_rows = a.size();
    int a_cols = a[0].size();
    int b_rows = b.size();
    int b_cols = b[0].size();

    if (a_cols != b_rows) 
    {
        throw invalid_argument("IMPOSSIBLE TO MULTIPLY");
    }
    vector<vector<double>> res(a_rows, vector<double>(b_cols, 0));

    for (int i = 0; i < a_rows; i++) 
    {
        for (int j = 0; j < b_cols; j++) 
        {
            for (int k = 0; k < a_cols; k++) 
            {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return res;
}

int main() 
{
    //ЧТЕНИЕ
    ios::sync_with_stdio(false);
    cin.tie(0);

    auto fx = freopen("input.txt", "r", stdin);
    auto fy = freopen("output.txt", "w", stdout);

    if(fx == NULL || fy == NULL)
    {
        cout << "ERROR: FILES NOT FOUND";
        return 2;
    }
    int n;
    cin >> n;
    vector<vector<double>> matrix(n, vector<double>(n));
    for (int i = 0; i < n; ++i) 
    {
        for (int j = 0; j < n; ++j) 
        {
            cin >> matrix[i][j];
        }
    }
    std::vector<double> coef(n);
    for (int i = 0; i < n; ++i) 
    {
        cin >> coef[i];
    }
    //ЧТЕНИЕ ЗАКОНЧЕНО

    //ПОЛУЧЕНИЕ ОПРЕДЕЛИТЕЛЯ И СТОЛБЦА Х
    auto [det, x] = std::make_tuple(0.0, std::vector<double>{});
    try 
    {
        std::tie(det, x) = gauss_calc_lineq(matrix, coef);
    } 
    catch (std::exception& e) 
    {
        std::cout << "ERROR:" << e.what();
        return 1;
    }

    std::vector<std::vector<double>> inverse;
    try 
    {
        inverse = gauss_inverse(matrix);
    } 
    catch (std::exception& e)
    {
        std::cout << "ERROR:" << e.what();
        return 1;
    }

    std::cout << "\n> Решение:\n";
    print(x);

    std::cout << "\n> A*x=COEFS:\n";

    std::vector<std::vector<double>> x_matrix(x.size(), std::vector<double>(1));

    for (int i = 0; i < x.size(); ++i) {
        x_matrix[i][0] = x[i];
    }
    try
    {
        print(mult_matrix(matrix, x_matrix));

        std::cout << "\n> Определитель:\ndet = " << det << '\n';

        std::cout << "\n> Обратная матрица:\n";
        print(inverse);

        std::cout << "\n> A*A^(-1) = E:\n";
        print(mult_matrix(matrix, inverse));
    }
    catch(const std::exception& e)
    {
        std::cout << "ERROR:" << e.what();
        return 1;
    }

}