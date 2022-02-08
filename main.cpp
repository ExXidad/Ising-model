#include <iostream>
#include "Solver.h"
#include <cmath>

int main()
{
    system("mkdir res");
    system("rm res/*");

    Solver solver(50, 1, 0.01, 0.0001);
    solver.temperatureSweep(-4, 4, 2000, pow(10., 6), true);
//    solver.temperatureSweep(-4, -3.5, 10, pow(10., 5), false;
    solver.solve(5*pow(10,5));
    return 0;
}
