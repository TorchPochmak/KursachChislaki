#ifndef QR_ALGORITHM_HPP
#define QR_ALGORITHM_HPP

#include <vector>
#include <cmath>
#include <complex>
#include <iostream>

// Simple matrix multiplication
std::vector<std::vector<double>> matrix_multiply(const std::vector<std::vector<double>>& A, 
                                               const std::vector<std::vector<double>>& B) {
    int n = A.size();
    std::vector<std::vector<double>> C(n, std::vector<double>(n, 0.0));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Matrix transpose
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& A) {
    int n = A.size();
    std::vector<std::vector<double>> AT(n, std::vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            AT[i][j] = A[j][i];
        }
    }
    return AT;
}

// QR decomposition using Gram-Schmidt process
void qr_decomposition(std::vector<std::vector<double>>& A, 
                     std::vector<std::vector<double>>& Q, 
                     std::vector<std::vector<double>>& R) {
    int n = A.size();
    Q = std::vector<std::vector<double>>(n, std::vector<double>(n));
    R = std::vector<std::vector<double>>(n, std::vector<double>(n, 0.0));

    // Copy A columns to Q as initial vectors
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            Q[i][j] = A[i][j];
        }
    }

    // Gram-Schmidt process
    for (int j = 0; j < n; j++) {
        for (int k = 0; k < j; k++) {
            // Calculate dot product
            double dot_product = 0.0;
            for (int i = 0; i < n; i++) {
                dot_product += Q[i][j] * Q[i][k];
            }
            R[k][j] = dot_product;

            // Subtract projection
            for (int i = 0; i < n; i++) {
                Q[i][j] -= R[k][j] * Q[i][k];
            }
        }

        // Normalize
        double norm = 0.0;
        for (int i = 0; i < n; i++) {
            norm += Q[i][j] * Q[i][j];
        }
        norm = std::sqrt(norm);

        if (norm > 1e-10) {  // Check for numerical stability
            R[j][j] = norm;
            for (int i = 0; i < n; i++) {
                Q[i][j] /= norm;
            }
        }
    }
}

// QR Algorithm to find eigenvalues
std::vector<std::complex<double>> find_eigenvalues(std::vector<std::vector<double>> A, 
                                                 int max_iter = 100, 
                                                 double tol = 1e-10) {
    int n = A.size();
    std::vector<std::complex<double>> eigenvalues(n);
    
    // Apply QR algorithm
    for (int iter = 0; iter < max_iter; iter++) {
        // Perform QR decomposition
        std::vector<std::vector<double>> Q, R;
        qr_decomposition(A, Q, R);
        
        // Check if Q and R are valid
        bool valid = true;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (!std::isfinite(Q[i][j]) || !std::isfinite(R[i][j])) {
                    valid = false;
                    break;
                }
            }
            if (!valid) break;
        }
        
        if (!valid) {
            std::cout << "Warning: Invalid values in QR decomposition at iteration " << iter << std::endl;
            break;
        }
        
        // A = RQ
        A = matrix_multiply(R, Q);
        
        // Check convergence
        bool converged = true;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (std::abs(A[i][j]) > tol) {
                    converged = false;
                    break;
                }
            }
            if (!converged) break;
        }
        
        if (converged) break;
    }
    
    // Extract eigenvalues from diagonal
    for (int i = 0; i < n; i++) {
        if (i < n - 1 && std::abs(A[i+1][i]) > tol) {
            // 2x2 block - complex conjugate pair
            double a = A[i][i];
            double b = A[i][i+1];
            double c = A[i+1][i];
            double d = A[i+1][i+1];
            
            double trace = a + d;
            double det = a*d - b*c;
            double disc = trace*trace - 4*det;
            
            if (disc < 0) {
                double real = trace/2;
                double imag = std::sqrt(-disc)/2;
                eigenvalues[i] = std::complex<double>(real, imag);
                eigenvalues[i+1] = std::complex<double>(real, -imag);
            } else {
                double sqrtDisc = std::sqrt(disc);
                eigenvalues[i] = std::complex<double>((trace + sqrtDisc)/2, 0);
                eigenvalues[i+1] = std::complex<double>((trace - sqrtDisc)/2, 0);
            }
            i++;
        } else {
            eigenvalues[i] = std::complex<double>(A[i][i], 0);
        }
    }
    
    return eigenvalues;
}

#endif // QR_ALGORITHM_HPP
