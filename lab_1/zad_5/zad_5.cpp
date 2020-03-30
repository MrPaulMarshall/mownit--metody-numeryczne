#include <iostream>
#include <cstdlib>

// bede analizowal wielomian 3 stopnia

struct Coefficients {
    double a;
    double b;
    double c;
    double d;
};

double f(Coefficients c, double x) {
    return c.a * x * x * x + c.b * x * x + c.c * x + c.d;
}

double df(Coefficients c, double x) {
    return 3 * c.a * x * x + 2 * c.b * x + c.c;
}

double d2f(Coefficients c, double x) {
    return 6 * c.a * x + 2 * c.b;
}

const double EPS = 0.0001;

double nextStep(Coefficients c, double x) {
    return x - f(c, x) / df(c, x);
}

void findZero(Coefficients c, double x) {
    std::cout << "For x0 = " << x << std::endl;
    int steps = 0;
    while (abs(f(c, x)) > EPS) {
        x = nextStep(c, x);
        steps++;
    }
    std::cout << "x = " << x << "; found after " << steps << " steps" << std::endl << std::endl;
}

int main()
{
    Coefficients c;
    c.a = 1.0;
    c.b = 0.0;
    c.c = -12.0;
    c.d = 17.0;

    findZero(c, -2.01);
    findZero(c, -1.99);
}
