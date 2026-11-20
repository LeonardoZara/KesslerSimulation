#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <cstdint>
#include <algorithm>
#include <filesystem>
#include <Eigen/Dense>
#include "KesslerSimulation.hpp"

// Funzioni definite in analisi.cpp
double lambdaAnnuale(const Eigen::Matrix4d& A, const Eigen::Vector4d& X, bool& valutabile);
void AnalisiMonteCarlo(const std::vector<std::vector<std::vector<double>>>& lambda,
                       const std::vector<std::vector<std::vector<char>>>& valutabile,
                       int M, int T, int numGusci,
                       const std::vector<double>& lanci, std::uint64_t seedBase);

int main() {
    std::filesystem::create_directories("output");
    //configurazione run:
    // M simulazioni Monte Carlo indipendenti, ciascuna di T anni, tutte dalla stessa configurazione iniziale e con gli stessi lanci annui costanti.
    const int    M_SIMULAZIONI = 50;
    const int    T_ANNI        = 20;
    const double DT            = 1.0; //anni

    // Seed di base. 0 = scelto a caso e stampato a terminale; Diverso da 0 = run riproducibile. La simulazione r usa il seme derivato da (SEED_BASE, r).
    const std::uint64_t SEED_BASE = 0;

    //4 gusci: 0 = Atmosfera, 1 = 300-400, 2 = 400-500, 3 = 500-600.
    const int numGusci = 4;

    // Lanci annui di satelliti attivi per guscio. Valori di riferimento: gamma * P0 = {9, 85, 4}.
    const std::vector<double> LANCI_P = {0.0, 9.0, 85.0, 4.0};

    // Sezione d'urto e massa dei satelliti: usate sia per la popolazione iniziale sia per i satelliti lanciati.
    // [MODIFICA] prima i lanciati avevano 130 m^2 e 800 kg (nel .cpp) mentre quelli iniziali 37.5 m^2 e 500 kg.
    const double SAT_SEZIONE = 37.5;    // m^2
    const double SAT_MASSA   = 500.0;   // kg

    // Fattore di scala della popolazione iniziale.
    const double scala = 1.0;

    KesslerSimulation sim(numGusci);

    auto n = [scala](int quantità){ return static_cast<int>(quantità * scala);};

    // CONFIGURAZIONE ORBITALE INIZIALE (modificabile a mano)
    // Parametri: Guscio, Quantità, Tipo, Sezione d'Urto (m^2), Massa (kg)
    // Satelliti: SAT_SEZIONE e SAT_MASSA. Frammenti: i valori comuni definiti in SpaceObject.hpp

    // 300-400 km
    sim.addObjectsToShell(1, n(945),   ObjectType::ACTIVE_PAYLOAD, SAT_SEZIONE, SAT_MASSA);
    sim.addObjectsToShell(1, n(140),   ObjectType::DEAD_PAYLOAD,   SAT_SEZIONE, SAT_MASSA);
    sim.addObjectsToShell(1, n(3000),  ObjectType::FRAGMENT, fragmentSigma, fragmentMass);

    // 400-500 km
    sim.addObjectsToShell(2, n(8520),  ObjectType::ACTIVE_PAYLOAD, SAT_SEZIONE, SAT_MASSA);
    sim.addObjectsToShell(2, n(1265),  ObjectType::DEAD_PAYLOAD,   SAT_SEZIONE, SAT_MASSA);
    sim.addObjectsToShell(2, n(12000), ObjectType::FRAGMENT, fragmentSigma, fragmentMass);

    // 500-600 km
    sim.addObjectsToShell(3, n(408),   ObjectType::ACTIVE_PAYLOAD, SAT_SEZIONE, SAT_MASSA);
    sim.addObjectsToShell(3, n(60),    ObjectType::DEAD_PAYLOAD,   SAT_SEZIONE, SAT_MASSA);
    sim.addObjectsToShell(3, n(85000), ObjectType::FRAGMENT, fragmentSigma, fragmentMass);

    // Lanci della run e copia immutabile della configurazione iniziale
    sim.setLaunches(LANCI_P, SAT_SEZIONE, SAT_MASSA);
    sim.saveInitialState();

    std::uint64_t seedBase = SEED_BASE;
    if (seedBase == 0) {
        std::random_device rd;
        seedBase = (static_cast<std::uint64_t>(rd()) << 32) | rd();
        if (seedBase == 0) seedBase = 1;
    }

    // Risultati: lambda[r][s][t] = lambda_max della matrice empirica del guscio s nell'anno t+1 della simulazione r;
    // valutabile[r][s][t] = 1 se in quell'anno la popolazione di frammenti era >= soglia (vedi analisi.cpp).
    std::vector<std::vector<std::vector<double>>> lambda(M_SIMULAZIONI,
        std::vector<std::vector<double>>(numGusci, std::vector<double>(T_ANNI, 0.0)));
    std::vector<std::vector<std::vector<char>>> valutabile(M_SIMULAZIONI,
        std::vector<std::vector<char>>(numGusci, std::vector<char>(T_ANNI, 0)));

    std::cout << "Avvio Monte Carlo: M = " << M_SIMULAZIONI << " simulazioni indipendenti da T = " << T_ANNI
              << " anni, seed base = " << seedBase << '\n';
    auto t0 = std::chrono::steady_clock::now();

    const int passoMessaggi = std::max(1, M_SIMULAZIONI / 10);

    for (int r = 0; r < M_SIMULAZIONI; ++r) {

        // Ripristino: gusci, tracker, tempo e seed di questa simulazione.
        sim.restoreInitialState(seedBase, static_cast<std::uint64_t>(r));

        for (int anno = 1; anno <= T_ANNI; ++anno) {

            sim.step(DT);

            // lambda_max dell'anno per ogni guscio, calcolato in memoria
            for (int s = 1; s < numGusci; ++s) {   // il guscio 0 (atmosfera) non ha significato
                bool val = false;
                lambda[r][s][anno - 1] = lambdaAnnuale(sim.getTracker().computeEmpiricalMatrix(s),
                                                       sim.getTracker().getInitialPopulation(s), val);
                valutabile[r][s][anno - 1] = val ? 1 : 0;
            }
        }

        if (r == 0 || (r + 1) % passoMessaggi == 0 || r + 1 == M_SIMULAZIONI) {
            double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::cout << "Completata simulazione " << (r + 1) << " / " << M_SIMULAZIONI
                      << "  (" << sec << " s)" << '\n';
        }
    }

    std::cout << "\nSimulazioni completate." << '\n';
    AnalisiMonteCarlo(lambda, valutabile, M_SIMULAZIONI, T_ANNI, numGusci, LANCI_P, seedBase);

    return 0;
}