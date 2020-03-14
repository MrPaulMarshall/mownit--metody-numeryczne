#include <iostream>

double machineEpsilon(double nextEPS) {
    double EPS = 1.0;

    while ((1 + nextEPS) != 1) {
        EPS = nextEPS;
        nextEPS /= 2;
    }

    return EPS;
}

int main()
{
    std::cout << "Machine epsilon = " << machineEpsilon(0.5) << std::endl;
}
