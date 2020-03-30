#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>

struct Raport {
	double x0;
	int number_of_iterations;
};

double f1(double x);
double f2(double x);
double f3(double x);

Raport bisection(double (*f)(double), double a, double b, const double EPS);
Raport newton(double (*f)(double), double seed, const double EPS, const int iter_limit);
Raport secants(double (*f)(double), double a, double b, const double EPS, const int iter_limit);
void printRaport(Raport raport, std::string name, double valueAtZero, const double precision);

int main() {
	std::vector<double> EPS = { 1e-7, 1e-15, 1e-33 };
	std::vector<double (*)(double)> functions = { f1, f2, f3 };


	// zad 2 - metoda bisekcji
	for (int i = 0; i < EPS.size(); i++) {
		std::cout << "Metoda: bisekcja" << std::endl << std::endl;

		Raport r1 = bisection(functions[0], 1.5 * M_PI, 2 * M_PI, EPS[i]);
		printRaport(r1, "f1", functions[0](r1.x0), EPS[i]);

		Raport r2 = bisection(functions[1], 0.0, 0.5 * M_PI, EPS[i]);
		printRaport(r2, "f2", functions[1](r2.x0), EPS[i]);

		Raport r3 = bisection(functions[2], 1.0, 3.0, EPS[i]);
		printRaport(r3, "f3", functions[2](r3.x0), EPS[i]);

		std::cout << std::endl << "-------------------------" << std::endl << std::endl;
	}

	// zad 3 - metoda Newtona
	const int max_iter = 10000000;
	for (int i = 0; i < EPS.size(); i++) {
		std::cout << "Metoda: Newton" << std::endl << std::endl;

		Raport r1 = newton(functions[0], 1.5 * M_PI, EPS[i], max_iter);
		printRaport(r1, "f1", functions[0](r1.x0), EPS[i]);

		Raport r2 = newton(functions[1], 0.1, EPS[i], max_iter);
		printRaport(r2, "f2", functions[1](r2.x0), EPS[i]);

		Raport r3 = newton(functions[2], 1.0, EPS[i], max_iter);
		printRaport(r3, "f3", functions[2](r3.x0), EPS[i]);

		std::cout << std::endl << "-------------------------" << std::endl << std::endl;
	}

	// zad 4 - metoda siecznych
	for (int i = 0; i < EPS.size(); i++) {
		std::cout << "Metoda: siecznych" << std::endl << std::endl;

		Raport r1 = secants(functions[0], 1.5 * M_PI, 2 * M_PI, EPS[i], max_iter);
		printRaport(r1, "f1", functions[0](r1.x0), EPS[i]);

		Raport r2 = secants(functions[1], 0.0, 0.5 * M_PI, EPS[i], max_iter);
		printRaport(r2, "f2", functions[1](r2.x0), EPS[i]);

		Raport r3 = secants(functions[2], 1.0, 3.0, EPS[i], max_iter);
		printRaport(r3, "f3", functions[2](r3.x0), EPS[i]);

		std::cout << std::endl << "-------------------------" << std::endl << std::endl;
	}

	return 0;
}


// ------------

// Funkcje testowe
double f1(double x) {
	return cos(x) * cosh(x) - 1;
}
double f2(double x) {
	return 1 / x - tan(x);
}
double f3(double x) {
	return pow(2, -x) + exp(x) + 2 * cos(x) - 6;
}

// -----------

void printRaport(Raport raport, std::string name, double valueAtZero, const double precision) {
	std::cout.precision(std::numeric_limits<double>::digits10 + 2);

	std::cout << "Miejsce zerowe funkcji " << name << ": " << raport.x0 << std::endl;
	std::cout << "Wyznaczone dla prezycji: " << precision << std::endl;
	std::cout << "Po " << raport.number_of_iterations << " iteracjach" << std::endl;
	std::cout << name << "(" << raport.x0 << ") = " << valueAtZero << std::endl << std::endl;
}

Raport bisection(double (*f)(double), double a, double b, const double EPS) {
	Raport raport;
	int iter = 0;

	// std::cout << "f(" << a << ") = " << (*f)(a) << std::endl;
	// std::cout << "f(" << b << ") = " << (*f)(b) << std::endl;

	if (abs((*f)(a)) < EPS) {
		raport.x0 = a;
		raport.number_of_iterations = 0;
		return raport;
	}
	if (abs((*f)(b)) < EPS) {
		raport.x0 = b;
		raport.number_of_iterations = 0;
		return raport;
	}

	if ((*f)(a) * (*f)(b) > 0) {
		raport.x0 = 0.0;
		raport.number_of_iterations = -1;
		return raport;
	}

	while (abs(b - a) > EPS) {
		iter++;
		double x = (a + b) / 2;
		// std::cout << "f(" << x << ") = " << (*f)(x) << std::endl;

		if (abs((*f)(x)) < EPS) {
			break;
		}
		else if ((*f)(a) * (*f)(x) < 0) {
			b = x;
		}
		else {
			a = x;
		}

		if (iter > 10000000) {
			std::cout << "After 10 milion iterations, still not found; x = " << x << std::endl;
			raport.x0 = x;
			raport.number_of_iterations = -1;
			return raport;
		}
	}
	raport.x0 = (a + b) / 2;
	raport.number_of_iterations = iter;

	return raport;
}

double diff(double (*f)(double), double x, const double h) {
	return ((*f)(x + h) - (*f)(x - h)) / (2*h);
}

Raport newton(double (*f)(double), double seed, const double EPS, const int iter_limit) {
	Raport raport;

	const double h = 1e-7;
	double x0 = seed;
	double x1 = x0 - (*f)(x0) / diff(f, x0, h);
	int iter = 0;

	do {
		x0 = x1;
		x1 = x0 - (*f)(x0) / diff(f, x0, h);
		iter++;
	} while (iter < iter_limit && abs((*f)(x1)) > EPS && abs(x1 - x0) > EPS);

	raport.x0 = x1;
	raport.number_of_iterations = iter;

	return raport;
}

Raport secants(double (*f)(double), double a, double b, const double EPS, const int iter_limit) {
	Raport raport;

	int iter = 0;
	double x0 = a;
	double x1 = b;
	double x2;

	do {
		iter++;
		x2 = ((*f)(x1) * x0 - (*f)(x0) * x1) / ((*f)(x1) - (*f)(x0));
		x0 = x1;
		x1 = x2;
	} while (iter < iter_limit && abs((*f)(x2)) > EPS && abs(x1 - x0) > EPS);

	raport.x0 = x2;
	raport.number_of_iterations = iter;

	return raport;
}
