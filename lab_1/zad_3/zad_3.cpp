#include <iostream>
#include <vector>

void prepareVector(std::vector<float>& vec, const int n);
void prepareVector(std::vector<double>& vec, const int n);

float sumForwFloat(const int n);
float sumBackFloat(const int n);
double sumForwDouble(const int n);
double sumBackDouble(const int n);

float sumKahanFloat(const int n);
double sumKahanDouble(const int n);

int main()
{
    std::vector<int> vec;
    vec.push_back(50);
    vec.push_back(100);
    vec.push_back(200);
    vec.push_back(500);
    vec.push_back(800);

    for (int i = 0; i < vec.size(); i++) {
        std::cout.precision(std::numeric_limits<double>::digits10 + 2);

        std::cout << "sumForwFloat(" << vec[i] << ") = " << sumForwFloat(vec[i]) << std::endl;
        std::cout << "sumBackFloat(" << vec[i] << ") = " << sumBackFloat(vec[i]) << std::endl;
        std::cout << "sumKahanFloat(" << vec[i] << ") = " << sumKahanFloat(vec[i]) << std::endl;
        std::cout << std::endl;

        std::cout << "sumForwDouble(" << vec[i] << ") = " << sumForwDouble(vec[i]) << std::endl;
        std::cout << "sumBackDouble(" << vec[i] << ") = " << sumBackDouble(vec[i]) << std::endl;
        std::cout << "sumKahanDouble(" << vec[i] << ") = " << sumKahanDouble(vec[i]) << std::endl;
        std::cout << std::endl << std::endl;
    }
    
    return 0;
}


void prepareVector(std::vector<float>& vec, const int n) {
    float aK = 1.0f;
    // jezeli liczymy dla aK = 1/2^(k+1)
    // to nie dodajemy elementu vec[0] = 0.5

    for (int k = 0; k < n + 1; k++) {
        aK /= 2;
        vec.push_back(aK);
    }
}

void prepareVector(std::vector<double>& vec, const int n) {
    double aK = 1.0;
    // jezeli liczymy dla aK = 1/2^(k+1)
    // to nie dodajemy elementu vec[0] = 0.5

    for (int k = 0; k < n + 1; k++) {
        aK /= 2;
        vec.push_back(aK);
    }
}

float sumForwFloat(const int n) {
    float sum = 0.0f;
    std::vector<float> vec;
    prepareVector(vec, n);
    for (int k = 1; k < n + 1; k++) {
        sum += vec[k];
    }

    return sum;
}

float sumBackFloat(const int n) {
    float sum = 0.0f;
    std::vector<float> vec;
    prepareVector(vec, n);
    for (int k = n; k > 0; k--) {
        sum += vec[k];
    }

    return sum;
}

double sumForwDouble(const int n) {
    double sum = 0.0f;
    std::vector<double> vec;
    prepareVector(vec, n);
    for (int k = 1; k < n + 1; k++) {
        sum += vec[k];
    }

    return sum;
}

double sumBackDouble(const int n) {
    double sum = 0.0;
    std::vector<double> vec;
    prepareVector(vec, n);
    for (int k = n; k > 0; k--) {
        sum += vec[k];
    }

    return sum;
}

float sumKahanFloat(const int n) {
    float sum = 0.0f;
    float err = 0.0f;
    std::vector<float> vec;
    prepareVector(vec, n);
    for (int k = 1; k < n + 1; k++) {
        float y = vec[k] - err;
        float temp = sum + y;
        err = (temp - y) - sum;
        sum = temp;
    }

    return sum;
}

double sumKahanDouble(const int n) {
    double sum = 0.0;
    double err = 0.0;
    std::vector<double> vec;
    prepareVector(vec, n);
    for (int k = 0; k < n + 1; k++) {
        double y = vec[k] - err;
        double temp = sum + y;
        err = (temp - y) - sum;
        sum = temp;
    }

    return sum;
}

/*
float sumForwFloat(const int n) {
    float sum = 0.0f;
    float aK = 0.5f;
    for (int k = 2; k <= n + 1; k++) {
        aK /= 2;
        sum += aK;
    }
    return sum;
}

float sumBackFloat(const int n) {
    float sum = 0.0f;
    float aK = 0.5f;
    for (int k = 2; k <= n + 1; k++) {
        aK /= 2;
    }
    for (int k = n + 1; k > 1; k--) {
        sum += aK;
        aK *= 2;
    }
    return sum;
}

double sumForwDouble(const int n) {
    double sum = 0.0;
    double aK = 0.5;
    for (int k = 2; k <= n + 1; k++) {
        aK /= 2;
        sum += aK;
    }
    return sum;
}

double sumBackDouble(const int n) {
    double sum = 0.0;
    double aK = 0.5;
    for (int k = 2; k <= n + 1; k++) {
        aK /= 2;
    }
    for (int k = n + 1; k > 1; k--) {
        sum += aK;
        aK *= 2;
    }
    return sum;
}
*/
