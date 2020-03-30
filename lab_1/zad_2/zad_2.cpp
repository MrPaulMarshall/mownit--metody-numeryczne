#include <iostream>
#include <vector>
#include <ctime>

void prepareVector(std::vector<float>& vec, float q, int n);

float absErr(float q, const int N, float sum);
float wglErr(float q, const int N, float sum);
void KahanSum(std::vector<float>& vec);

float executionTime(clock_t start, clock_t end);


int main()
{
    std::vector<float> vec;
    float q = 0.3f;
    const int N = 10000000;

    clock_t start, end;
    float timeKahan;

    prepareVector(vec, q, N);

    start = clock();
    KahanSum(vec);
    end = clock();
    timeKahan = executionTime(start, end);

    std::cout << "Czas wykonania dodawania ala Kahan: " << timeKahan << " sekund\n";

    return 0;
}


void prepareVector(std::vector<float>& vec, float q, int n) {
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
    float realSum = q * (N + 1);
    return abs(sum - realSum) / realSum;
}

void KahanSum(std::vector<float>& vec) {
    float sum = 0.0f;
    float err = 0.0f;
    for (int i = 0; i < vec.size(); i++) {
        float y = vec[i] - err;
        float temp = sum + y;
        err = (temp - sum) - y;
        sum = temp;
    }

    std::cout.precision(std::numeric_limits<float>::digits10 + 2);
    std::cout << "Suma Kahan: " << sum << "\n";
    std::cout << "Blad Kahan: " << err << "\n";
    std::cout << "Blad wzgledny Kahan: " << abs(err) / (vec[0] * vec.size()) << "\n";
    std::cout << "Blad bezwgledny: " << absErr(vec[0], vec.size(), sum) << "\n";
    std::cout << "Blad wzgledny: " << wglErr(vec[0], vec.size(), sum) << "\n";
    std::cout << std::endl;
}
