#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>

void prepareVector(std::vector<float>& vec, float q, const int n);

float absErr(float q, const int N, float sum);
float wglErr(float q, const int N, float sum);
void normalSum(std::vector<float>& vec);
float recursiveSumInner(std::vector<float>& vec, int start, int end);
void recursiveSum(std::vector<float>& vec);

float executionTime(clock_t start, clock_t end);

int main() {
    std::vector<float> vec;
    const int N = 10000000;
    float q = 0.3f;

    clock_t start, end;
    float timeNormal, timeRecursive;

    prepareVector(vec, q, N);
    std::cout.precision(std::numeric_limits<float>::digits10 + 2);

    start = clock();
    normalSum(vec);
    end = clock();
    timeNormal = executionTime(start, end);

    start = clock();
    recursiveSum(vec);
    end = clock();
    timeRecursive = executionTime(start, end);

    std::cout << "Czas wykonania zwyklego dodawania: " << timeNormal << " sekund\n";
    std::cout << "Czas wykonania rekurencyjnego dodawania: " << timeRecursive << " sekund\n";

    return 0;
}


void prepareVector(std::vector<float>& vec, float q, const int n) {
    for (int i = 0; i < n; i++) {
        vec.push_back(q);
    }
}

float executionTime(clock_t start, clock_t end) {
    return (float)(end - start) / CLOCKS_PER_SEC;
}

float absErr(float q, const int N, float sum) {
    float realSum = q * N;
    return abs(sum - realSum);
}

float wglErr(float q, const int N, float sum) {
    float realSum = (N + 1) * q;
    return abs(sum - realSum) / realSum;
}

void normalSum(std::vector<float>& vec) {
    float sum = 0.0f;

    std::ofstream zapis("dane.txt");

    for (int i = 0; i < vec.size(); i++) {
        sum += vec[i];
        if ((i+1) % 25000 == 0) {
            // std::cout << "i = " << i << "\n";
            // std::cout << "Bezwzgledny: " << absErr(vec[0], i, sum) << "\n";
            // std::cout << "Wzgledny:    " << wglErr(vec[0], i, sum) << "\n";
        
            // zapis do pliku
            zapis << (i+1) << " " << wglErr(vec[0], i, sum) << "\n";
        }
    }

    zapis.close();

    std::cout << std::endl;
    std::cout << "Suma po kolei: " << sum << "\n";
    std::cout << "Bezwzgledny: " << absErr(vec[0], vec.size(), sum) << "\n";
    std::cout << "Wzgledny:    " << wglErr(vec[0], vec.size(), sum) << "\n\n";
}


float recursiveSumInner(std::vector<float>& vec, int start, int end) {
    if (start > end) {
        return 0.0f;
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

void recursiveSum(std::vector<float>& vec) {
    float sum = recursiveSumInner(vec, 0, vec.size() - 1);
    std::cout << "Suma rekurencyjnie: " << sum << "\n";
    std::cout << "Blad bezwgledny: " << absErr(vec[0], vec.size(), sum) << "\n";
    std::cout << "Blad wzgledny: " << wglErr(vec[0], vec.size(), sum) << "\n\n";
}
