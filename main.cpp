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
bool annoInstabile(const Eigen::Matrix4d& A, const Eigen::Vector4d& X);
void AnalisiMonteCarlo(const std::vector<std::vector<int>>& anniInstabili,
                       int M, int T, int numGusci,
                       const std::vector<double>& lanci, std::uint64_t seedBase, int kAnniMinimi);

int main() {
    std::filesystem::create_directories("output");
    
    // M simulazioni Monte Carlo indipendenti, ciascuna di T anni, tutte dalla stessa configurazione iniziale e con gli stessi lanci annui (costanti negli anni).
    const int M_SIMULAZIONI = 20;
    const int T_ANNI = 10;
    const double DT = 1.0;   // 1 anno, non cambiare

    // Una simulazione e' "instabile" se ha almeno K_ANNI_MINIMI anni con lambda_max > 1 (anno instabile).
    const int K_ANNI_MINIMI = 1;

    // Seed di base. 0 = scelto a caso e stampato a video (cosi' la run si puo' rifare);
    const std::uint64_t SEED_BASE = 0;

    const int numGusci = 4;

    const std::vector<double> LANCI_P = {0.0, 304.0, 2744.0, 131.0};

    const double SAT_SEZIONE = 37.5;    // m^2
    const double SAT_MASSA   = 500.0;   // kg

    // Fattore di scala della popolazione iniziale
    const double scala = 1.0;

    KesslerSimulation sim(numGusci);

    auto n = [scala](int quantità){ return static_cast<int>(quantità * scala);};

    // Configurazione orbitale iniziale:
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

    sim.setLaunches(LANCI_P, SAT_SEZIONE, SAT_MASSA);
    sim.saveInitialState();

    std::uint64_t seedBase = SEED_BASE;
    if (seedBase == 0) {
        std::random_device rd;
        seedBase = (static_cast<std::uint64_t>(rd()) << 32) | rd();
        if (seedBase == 0) seedBase = 1;
    }

    // Risultato di ogni simulazione: anniInstabili[s][r] = quanti anni (su T) il guscio s della simulazione r e' stato instabile
    std::vector<std::vector<int>> anniInstabili(numGusci, std::vector<int>(M_SIMULAZIONI, 0));

    std::cout << "Avvio Monte Carlo: M = " << M_SIMULAZIONI << " simulazioni indipendenti da T = " << T_ANNI
              << " anni, seed base = " << seedBase << '\n';
    auto t0 = std::chrono::steady_clock::now();

    const int passoMessaggi = std::max(1, M_SIMULAZIONI / 10);

    for (int r = 0; r < M_SIMULAZIONI; ++r) {

        // Ripristino: gusci, tracker, tempo e RNG (seed di questa simulazione)
        sim.restoreInitialState(seedBase, static_cast<std::uint64_t>(r));

        for (int anno = 1; anno <= T_ANNI; ++anno) {

            sim.step(DT);

            // Conta gli anni instabili di ogni guscio (lambda_max > 1)
            for (int s = 1; s < numGusci; ++s) {
                if (annoInstabile(sim.getTracker().computeEmpiricalMatrix(s),
                                  sim.getTracker().getInitialPopulation(s))) {
                    anniInstabili[s][r]++;
                }
            }
        }

        if (r == 0 || (r + 1) % passoMessaggi == 0 || r + 1 == M_SIMULAZIONI) {
            double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            std::cout << "Completata simulazione " << (r + 1) << " / " << M_SIMULAZIONI
                      << "  (" << sec << " s)" << '\n';
        }
    }

    std::cout << "\nSimulazioni completate." << '\n';
    AnalisiMonteCarlo(anniInstabili, M_SIMULAZIONI, T_ANNI, numGusci, LANCI_P, seedBase, K_ANNI_MINIMI);

    return 0;
}