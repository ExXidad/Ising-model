//
// Created by Ivan Kalesnikau on 07.02.2022.
//

#ifndef MCMC_SOLVER_H
#define MCMC_SOLVER_H

#include <random>
#include <iostream>
#include <ios>
#include <fstream>
#include <string>

class Solver
{
private:
    int n;
    double J, H, T;
    int **spinGrid;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_int_distribution<> uid;
    std::uniform_real_distribution<> urd;

public:


private:
    void initialSeed();

    int randomSpin();

    int periodicBC(const int &j);

    double energy(const int &j, const int &i);

    double crossEnergy(const int &j, const int &i);

    void flip(const int &j, const int &i);


public:
    Solver(const int &n, const double &J, const double &H, const double &T);

    ~Solver();

    void solve(const int &N);

    void save(const std::string &name);

    void temperatureSweep(const double &TMin, const double &TMax, const int &nop, const int &N, const bool &exp);
};


#endif //MCMC_SOLVER_H
