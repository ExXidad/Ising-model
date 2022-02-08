//
// Created by Ivan Kalesnikau on 07.02.2022.
//

#include "Solver.h"

Solver::Solver(const int &n, const double &J, const double &H, const double &T)
{
    this->n = n;
    this->J = J;
    this->H = H;
    this->T = T;

    spinGrid = new int *[n];
    for (int i = 0; i < n; ++i)
        spinGrid[i] = new int[n];

    gen = std::mt19937(rd());
    uid = std::uniform_int_distribution<>(0, 1);
    urd = std::uniform_real_distribution<>(0, 1);
}

Solver::~Solver()
{
    for (int i = 0; i < n; ++i)
        delete spinGrid[i];
    delete spinGrid;
}

int Solver::randomSpin()
{
    return uid(gen) == 1 ? 1 : -1;
}

void Solver::initialSeed()
{
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            spinGrid[j][i] = randomSpin();
        }
    }
}

void Solver::solve(const int &N)
{
    initialSeed();

    int saveEveryXIterations = std::max(N / 100, 1);

    std::fstream file;
    file.open("energy", std::ios::out);

    for (int iter = 0; iter < N + saveEveryXIterations; ++iter)
    {
        int j = urd(gen) * n, i = urd(gen) * n;

        double eOld = crossEnergy(j, i);

        if (iter % saveEveryXIterations == 0)
        {
            save(std::to_string(iter));

            double meanEnergyPerSpin = 0;
            double meanSpin = 0;
            for (int k = 0; k < n; ++k)
                for (int l = 0; l < n; ++l)
                {
                    meanEnergyPerSpin += energy(k, l);
                    meanSpin += spinGrid[k][l];
                }

            meanEnergyPerSpin /= n * n;
            meanSpin /= n * n;
            file << iter << "\t" << meanEnergyPerSpin << "\t" << meanSpin << "\t" << meanSpin / H << std::endl;
        }


        flip(j, i);

        double eNew = crossEnergy(j, i);

        double acceptProbability = std::min(1., std::exp((eOld - eNew) / T));

        urd(gen) <= acceptProbability ? void() : flip(j, i);
    }

    file.close();

}

int Solver::periodicBC(const int &j)
{
    return std::abs(j % n);
}

double Solver::energy(const int &j, const int &i)
{
    double e = -H * spinGrid[j][i];
    for (int k = -1; k <= 1; ++k)
    {
        for (int l = -1; l <= 1; ++l)
        {
            if (k xor l) e += -J / 2 * spinGrid[j][i] * spinGrid[periodicBC(j + k)][periodicBC(i + l)];
        }
    }
    return e;
}

void Solver::save(const std::string &name)
{
    std::fstream file;
    file.open("res/" + name, std::ios::out);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            file << spinGrid[j][i] << "\t";
        }
        file << std::endl;

    }
    file.close();

}

double Solver::crossEnergy(const int &j, const int &i)
{
    double e = 0;
    for (int k = -1; k <= 1; ++k)
    {
        for (int l = -1; l <= 1; ++l)
        {
            if (k * l == 0) e += energy(periodicBC(j + k), periodicBC(i + l));
        }
    }
    return e;
}

void Solver::flip(const int &j, const int &i)
{
    spinGrid[j][i] = spinGrid[j][i] == 1 ? -1 : 1;
}

void
Solver::temperatureSweep(const double &TMin, const double &TMax, const int &nop, const int &N, const bool &exp)
{
    double dT = (TMax - TMin) / double(nop);
    double oldT = this->T;

    int separator = std::max(N / 5, 1);

    std::fstream file;
    file.open("temperatureSweep", std::ios::out);

    int counter = 0;
    double tmp = 0;

    for (double T = TMin; T <= TMax + 0.5 * dT; T += dT)
    {
        initialSeed();
        this->T = exp ? pow(10., T) : T;
        double meps = 0, ms = 0, mms = 0;

        for (int iter = 0; iter <= N; ++iter)
        {
            int j = urd(gen) * n, i = urd(gen) * n;

            double eOld = crossEnergy(j, i);

            if (iter >= N - separator)
            {
                double meanEnergyPerSpin = 0;
                double meanSpin = 0;
                for (int k = 0; k < n; ++k)
                    for (int l = 0; l < n; ++l)
                    {
                        meanEnergyPerSpin += energy(k, l);
                        meanSpin += spinGrid[k][l];
                    }

                meanEnergyPerSpin /= n * n;
                meanSpin /= n * n;
                meps += meanEnergyPerSpin;
                ms += meanSpin;
                mms += meanSpin / H;
                ++counter;
                tmp = meanEnergyPerSpin;
            }

            flip(j, i);

            double eNew = crossEnergy(j, i);

            double acceptProbability = std::min(1., std::exp((eOld - eNew) / this->T));

            urd(gen) <= acceptProbability ? void() : flip(j, i);
        }

        file << (exp ? pow(10., T) : T) << "\t" << meps / double(separator + 1) << "\t" << ms / double(separator + 1) << "\t"
             << mms / H / double(separator + 1) << std::endl;
        std::cout << TMin << ":" << T << ":" << TMax << " through " << dT << std::endl;
        counter = 0;


    }
    file.close();

    this->T = oldT;
}
