#include <iostream>
#include <cmath>
#include <cstdlib>

double f1(double x) {
	return cos(x) * cosh(x) - 1;
}
double f2(double x) {
	return 1 / x - tan(x);
}
double f3(double x) {
	return pow(2, -x) + exp(x) + 2 * cos(x) - 6;
}

struct Raport {
	double x0;
	int number_of_iterations;
};

Raport bisection(double (*f)(double), double a, double b, const double EPS);

int main() {
	Raport r1 = bisection(f1, 1.5 * M_PI, 2 * M_PI, 1e-7);
	std::cout << "Miejsce zerowe: " << r1.x0 << std::endl;
	std::cout << "Po " << r1.number_of_iterations << " iteracjach" << std::endl;

	return 0;
}

Raport bisection(double (*f)(double), double a, double b, const double EPS) {
	Raport raport;
	int iter = 0;

	if (abs((*f)(a)) < EPS) {
		raport.x0 = a;
		raport.number_of_iterations = 0;
	}
	if (abs((*f)(b)) < EPS) {
		raport.x0 = b;
		raport.number_of_iterations = 0;
	}

	if ((*f)(a) * (*f)(b) > 0) {
		raport.x0 = 0.0;
		raport.number_of_iterations = -1;
		return raport;
	}

	while (abs(b - a) < EPS) {
		iter++;
		double x = (a + b) / 2;

		if (abs((*f)(x)) < EPS) {
			break;
		}
		else if ((*f)(a) * (*f)(x) < 0) {
			b = x;
		}
		else {
			a = x;
		}
	}
	raport.x0 = (a + b) / 2;
	raport.number_of_iterations = iter;

	return raport;
}
