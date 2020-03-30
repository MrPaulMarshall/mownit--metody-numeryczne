#include "aghMatrix.h"
#include <iostream>
#include <vector>
#include <cmath>

void testAdditionAndMultiplication();

void testFactorizationLU();

void testFactorizationCholesky();

void testGauss();

void testJacobi();

void factorize_LU(AGHMatrix<double>& matrix, AGHMatrix<double>& L, AGHMatrix<double>& U);

void factorize_Cholesky(AGHMatrix<double>& matrix, AGHMatrix<double>& L, AGHMatrix<double>& LT);

bool solve_Gauss(AGHMatrix<double>& A, AGHMatrix<double>& X);

AGHMatrix<double> solve_Jacobi(AGHMatrix<double>& A, AGHMatrix<double>& B, const int iterations);

int main()
{
    testAdditionAndMultiplication();
    testFactorizationLU();
    testFactorizationCholesky();
    testGauss();
    testJacobi();

    return 0;
}


// functions to test my implementation

void testAdditionAndMultiplication() {
    // initialize matrices using init value
    AGHMatrix<double> mat1(5, 5, 1.2);
    AGHMatrix<double> mat2(5, 5, 2.8);

    std::cout << mat1;
    std::cout << mat2;

    // element of Chaos
    mat1(2, 3) = -100.0;

    // Uncomment when implemented
    AGHMatrix<double> matSum = mat1 + mat2;
    std::cout << matSum;

    AGHMatrix<double> matProduct = mat1 * mat2;
    std::cout << matProduct;

    std::cout << "--------------------------\n\n";
}

void testFactorizationLU() {
    // initialize matrix using specified values - LU
    std::cout << "LU starts\n\n";

    std::vector<std::vector<double>> init_LU{ { 1.0, 2.0, 3.0 },
                                            { 4.0, 5.0, 6.0 },
                                            { 7.0, 8.0, 9.0 } };
    AGHMatrix<double> A(init_LU);
    std::cout << A;

    AGHMatrix<double> L(A.get_rows(), A.get_cols(), 0);
    AGHMatrix<double> U(A.get_rows(), A.get_cols(), 0);

    factorize_LU(A, L, U);

    std::cout << L;
    std::cout << U;

    std::cout << (L * U);
    std::cout << "A == L * U: " << (A == (L * U)) << std::endl;

    std::cout << "LU ends\n\n";

    std::cout << "--------------------------\n\n";
}

void testFactorizationCholesky() {
    // initialize matrix using specified values - Choleski
    std::cout << "Cholesky starts\n\n";

    std::vector<std::vector<double>> init_Cholesky{ { 4.0, 12.0, -16.0 },
                                                { 12.0, 37.0, -43.0 },
                                                { -16.0, -43.0, 98.0 } };

    AGHMatrix<double> A(init_Cholesky);
    std::cout << A;

    AGHMatrix<double> L(A.get_rows(), A.get_cols(), 0);
    AGHMatrix<double> LT(A.get_rows(), A.get_cols(), 0);

    factorize_Cholesky(A, L, LT);

    std::cout << L;
    std::cout << LT;

    std::cout << (L * LT);
    std::cout << "A == L * LT: " << (A == (L * LT)) << std::endl;

    std::cout << "Cholesky ends\n\n";

    std::cout << "--------------------------\n\n";
}

void testGauss() {
    // eliminacja Gaussa
    std::cout << "Gauss starts\n\n";

    std::vector<std::vector<double>> init{ {1.0010, -5.030, 5.8090, 7.8320, 9.5740},
                                           {2.2660, 1.9950, 1.2120, 8.0080, 7.2190},
                                           {8.8500, 5.6810, 4.5520, 1.3020, 5.7300},
                                           {6.7750, -2.253, 2.9080, 3.9700, 6.2910} };
    AGHMatrix<double> AB(init);
    std::cout << AB;

    AGHMatrix<double> solution(AB.get_rows(), 1, 0);
    if (solve_Gauss(AB, solution) == true) {
        std::cout << AB;

        std::cout << solution;
    }

    std::vector<std::vector<double>> init_RealSol{
        {0.23312286}, {-0.00465772}, {0.59793688}, {0.74617086}
    };
    AGHMatrix<double> realSolution(init_RealSol);

    if (solution == realSolution) {
        std::cout << "Rozwiazanie jest poprawne\n\n";
    }
    else {
        std::cout << "Rozwiazanie jest bledne\n\n";
    }

    std::cout << "Gauss ends\n\n";

    std::cout << "--------------------------\n\n";
}

void testJacobi() {
    std::cout << "Jacobi starts\n\n";

    // initialization of matrices A and B
    /*
    std::vector<std::vector<double>> initA{ {1.0010, -5.030, 5.8090},
                                            {2.2660, 1.9950, 1.2120},
                                            {6.7750, -2.253, 2.9080} };

    std::vector<std::vector<double>> initB{ {9.5740},
                                            {7.2190},
                                            {6.2910} };
    */
    std::vector<std::vector<double>> initA{ {1.0010, -5.030, 5.8090, 7.8320},
                                            {2.2660, 1.9950, 1.2120, 8.0080},
                                            {8.8500, 5.6810, 4.5520, 1.3020},
                                            {6.7750, -2.253, 2.9080, 3.9700} };
    
    std::vector<std::vector<double>> initB{ {9.5740},
                                            {7.2190},
                                            {5.7300},
                                            {6.2910} };
                                            
    AGHMatrix<double> A(initA);
    std::cout << A;

    AGHMatrix<double> B(initB);
    std::cout << B;

    AGHMatrix<double> solution = solve_Jacobi(A, B, 100);
    std::cout << solution;

    std::cout << "Jacobi ends\n\n";

    std::cout << "--------------------------\n\n";
}


// ---------------------------------
// Here are implementations

void factorize_LU(AGHMatrix<double>& matrix, AGHMatrix<double>& L, AGHMatrix<double>& U) {
    if (matrix.get_rows() != matrix.get_cols()) {
        std::cout << "Matrix is not NxN\n";
        return;
    }

    unsigned n = matrix.get_rows();
    // ustawiam elementy diagonali glownej macierzy L na 1
    for (int j = 0; j < n; j++) {
        L(j, j) = 1;
    }
    
    for (int i = 0; i < n; i++) {
        // obliczam i-ty wiersz macierzy U
        for (int j = i; j < n; j++) {
            U(i, j) = matrix(i, j);
            for (int k = 0; k < i; k++) {
                U(i, j) -= L(i, k) * U(k, j);
            }
        }

        // obliczam i-tą kolumne macierzy L
        for (int j = i + 1; j < n; j++) {
            if (U(i, i) == 0) {
                std::cout << "Nie mozna dokonac rozkladu LU: 0 na diagonali glownej macierzy U\n";
                return;
            }
            L(j, i) = matrix(j, i);
            for (int k = 0; k < i; k++) {
                L(j, i) -= L(j, k) * U(k, i);
            }
            L(j, i) = L(j, i) / U(i, i);
        }
    }
}

void factorize_Cholesky(AGHMatrix<double>& matrix, AGHMatrix<double>& L, AGHMatrix<double>& LT) {
    if (matrix.isSymmetric() == false) {
        std::cout << "Matrix is not symmetric - factorization ala Cholesky is not possible\n";
        return;
    }

    unsigned n = matrix.get_rows();

    for (int i = 0; i < n; i++) {
        double elem = matrix(i, i);
        for (int k = 0; k < i; k++) {
            elem -= L(i, k) * L(i, k);
        }
        LT(i, i) = L(i, i) = sqrt(elem);
    
        if (L(i, i) == 0) {
            std::cout << "0 na diagonali glownej L: macierz A nie jest dodatnio okreslona\n";
        }
        for (int j = i + 1; j < n; j++) {
            elem = matrix(j, i);
            for (int k = 0; k < i; k++) {
                elem -= L(j, k) * L(i, k);
            }
            LT(i, j) = L(j, i) = elem / L(i, i);
        }
    }
}

bool solve_Gauss(AGHMatrix<double>& A, AGHMatrix<double>& X) {
    if (A.get_rows() + 1 != A.get_cols()) {
        std::cout << "Liczba rownan i niewiadomych nie jest zgodna\n";
        return false;
    }

    unsigned n = A.get_rows();

    if (X.get_rows() != n || X.get_cols() != 1) {
        std::cout << "Zle wymiary wektora na rozwiazanie\n";
        return false;
    }

    // Precision of comparing to 0.0
    const double EPS = 1e-7;

    // Coefficients, made only to avoid executing the same division multiple times
    double m, s;

    // Firstly, elimination of elements below main diagonal, top-bottom
    for (int i = 0; i < n - 1; i++) {
        // For each row
        for (int j = i + 1; j < n; j++) {
            // If 0 on the main diagonal
            if (abs(A(i, i)) < EPS) {
                std::cout << "Dzielenie przez 0\n";
                return false;
            }
            
            // Addition of normalized i-th row to j-th row, to make sure that below A(i, i) there are only 0s
            m = -A(j, i) / A(i, i);
            for (int k = i; k < n + 1; k++) {
                A(j, k) += m * A(i, k);
            }
        }
    }

    // Then, calculation of searched variables, bottom-up
    // For each variable ~ row
    for (int i = n - 1; i >= 0; i--) {
        // If:      A(i, i)*x(i) + Sum( A(i, j)*x(j) ) = B(i)
        // Then:    x(i) = ( B(i) - Sum(A(i, j)*x(j) ) / A(i, i)
        s = A(i, n);
        for (int j = n - 1; j > i; j--) {
            s -= A(i, j) * X(j, 0);
        }

        if (abs(A(i, i)) < EPS) {
            std::cout << "Dzielenie przez 0\n";
            return false;
        }

        X(i, 0) = s / A(i, i);
    }

    // Successfully calculated solution
    return true;
}

AGHMatrix<double> solve_Jacobi(AGHMatrix<double>& A, AGHMatrix<double>& B, const int iterations) {
    if (A.get_rows() != A.get_cols() || A.get_rows() != B.get_rows() || B.get_cols() != 1) {
        std::cout << "Zle wymiary macierzy A lub B\n";
        exit(1);
    }

    unsigned n = A.get_rows();

    // create all needed matrices - initialization with 0s
    AGHMatrix<double> L(n, n, 0);
    AGHMatrix<double> U(n, n, 0);
    AGHMatrix<double> N(n, n, 0);

    // first estimation is vector of 0s
    AGHMatrix<double> x(n, 1, 0);

    // decompose: A = L + D + U
    //      matrix D is defined as D(i, i) = A(i, i), so I don't declare it
    for (int i = 0; i < n; i++) {
        if (abs(A(i, i)) < cEPS) {
            std::cout << "0 na przekatnej glownej\n";
            exit(1);
        }

        for (int j = i + 1; j < n; j++) {
            L(j, i) = A(j, i);
            U(i, j) = A(i, j);
        }
    }

    // matrix N = D^(-1), but since D is diagonal, N is also diagonal, with elements defined as below
    for (int i = 0; i < n; i++) {
        N(i, i) = 1 / A(i, i);
    }

    // std::cout << "L:\n" << L << "N:" << N << "U:\n" << U;

    AGHMatrix<double> M = N * (L + U);
    // M = -N * (L + U), so I multiply variable M by -1
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            M(i, j) *= -1;

    // finally I apply recisive formula iterations-times
    for (int i = 0; i < iterations; i++)
        x = M * x + N * B;

    return x;
}
