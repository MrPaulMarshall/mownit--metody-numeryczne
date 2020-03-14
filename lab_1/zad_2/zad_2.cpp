#include <iostream>
#include <vector>
#include <ctime>

void prepareVector(std::vector<double>& vec, double q, int n);

double absErr(double q, const int N, double sum);
double wglErr(double q, const int N, double sum);
void normalSum(std::vector<double>& vec);
double recursiveSumInner(std::vector<double>& vec, int start, int end);
void recursiveSum(std::vector<double>& vec);
void KahanSum(std::vector<double>& vec);

double executionTime(clock_t start, clock_t end);


int main()
{
    std::vector<double> vec;
    double q = 0.3;
    const int N = 10000000;

    clock_t start, end;
    double timeRecursive, timeKahan;

    prepareVector(vec, q, N);

    start = clock();
    recursiveSum(vec);
    end = clock();
    timeRecursive = executionTime(start, end);

    start = clock();
    KahanSum(vec);
    end = clock();
    timeKahan = executionTime(start, end);

    std::cout << "Czas wykonania rekurencyjnego dodawania: " << timeRecursive << " sekund\n";
    std::cout << "Czas wykonania dodawania ala Kahan: " << timeKahan << " sekund\n";

    return 0;
}


void prepareVector(std::vector<double>& vec, double q, int n) {
    for (int i = 0; i < n; i++) {
        vec.push_back(q);
    }
}

double executionTime(clock_t start, clock_t end) {
    return (double)(end - start) / CLOCKS_PER_SEC;
}

double absErr(double q, const int N, double sum) {
    double realSum = q * N;
    return abs(sum - realSum);
}

double wglErr(double q, const int N, double sum) {
    double realSum = q * (N + 1);
    return abs(sum - realSum) / realSum;
}

void normalSum(std::vector<double>& vec) {
    double sum = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        sum += vec[i];
        if (i % 25000 == 0) {
            std::cout << "i = " << i << "\n";
            std::cout << "Bezwzgledny: " << absErr(vec[0], i, sum) << "\n";
            std::cout << "Wzgledny:    " << wglErr(vec[0], i, sum) << "\n";
        }
    }
    std::cout << "Suma po kolei: " << sum << "\n\n";
}

double recursiveSumInner(std::vector<double>& vec, int start, int end) {
    if (start > end) {
        return 0.0;
    }
    if (start == end) {
        return vec[start];
    }
    if (start + 1 == end) {
        return vec[start] + vec[end];
    }

    int mid = (start + end) / 2;
    return recursiveSumInner(vec, start, mid) + recursiveSumInner(vec, mid + 1, end);
}

void recursiveSum(std::vector<double>& vec) {
    double sum = recursiveSumInner(vec, 0, vec.size() - 1);

    std::cout.precision(std::numeric_limits<double>::digits10 + 2);
    std::cout << "Suma rekurencyjnie: " << sum << "\n";
    std::cout << "Blad bezwgledny: " << absErr(vec[0], vec.size(), sum) << "\n";
    std::cout << "Blad wzgledny: " << wglErr(vec[0], vec.size(), sum) << "\n";
    std::cout << std::endl;
}

void KahanSum(std::vector<double>& vec) {
    double sum = 0.0;
    double err = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        double y = vec[i] - err;
        double temp = sum + y;
        err = (temp - sum) - y;
        sum = temp;
    }

    std::cout.precision(std::numeric_limits<double>::digits10 + 2);
    std::cout << "Suma Kahan: " << sum << "\n";
    std::cout << "Blad bezwgledny: " << absErr(vec[0], vec.size(), sum) << "\n";
    std::cout << "Blad wzgledny: " << wglErr(vec[0], vec.size(), sum) << "\n";
    std::cout << std::endl;
}
