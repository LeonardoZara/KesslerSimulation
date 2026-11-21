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

    // =====================================================================
    //                 CONFIGURAZIONE DELLA RUN (modificare qui)
    // =====================================================================
    // M simulazioni Monte Carlo INDIPENDENTI, ciascuna di T anni, tutte dalla stessa
    // configurazione iniziale e con gli stessi lanci annui (costanti negli anni).
    const int    M_SIMULAZIONI = 200;
    const int    T_ANNI        = 50;
    const double DT            = 1.0;   // passo = 1 anno. Non cambiarlo: l'analisi e i lanci assumono dt = 1

    // Una simulazione e' "instabile" se ha almeno K_ANNI_MINIMI anni con lambda_max > 1 (anno instabile).
    // Scelta fissata a 1: basta un anno instabile. Tenere lo stesso valore per tutte le run che si confrontano.
    const int    K_ANNI_MINIMI = 1;

    // Seed di base. 0 = scelto a caso e STAMPATO a video (cosi' la run si puo' rifare);
    // un valore diverso da 0 = run riproducibile. La simulazione r usa il seme derivato da (SEED_BASE, r).
    const std::uint64_t SEED_BASE = 0;

    //4 gusci: 0 = Atmosfera, 1 = 300-400, 2 = 400-500, 3 = 500-600.
    const int numGusci = 4;

    // Lanci annui di satelliti attivi per guscio (indice 0 = atmosfera, ignorato).
    // Interi: il codice arrotonda a intero (lround) i lanci di ogni anno.
    // Valori di riferimento: gamma * P0 = {9, 85, 4}. Cambiarli a mano tra una run e l'altra.
    const std::vector<double> LANCI_P = {0.0, 9.0, 85.0, 4.0};

    // Sezione d'urto e massa dei satelliti: usate sia per la popolazione iniziale sia per i satelliti lanciati.
    // [MODIFICA] prima i lanciati avevano 130 m^2 e 800 kg (nel .cpp) mentre quelli iniziali 37.5 m^2 e 500 kg.
    const double SAT_SEZIONE = 37.5;    // m^2
    const double SAT_MASSA   = 500.0;   // kg

    // Fattore di scala della popolazione iniziale. 1.0 = dati del PDF.
    const double scala = 1.0;
    // =====================================================================

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

    // Lanci della run e copia immutabile della configurazione iniziale (da qui in poi non si aggiunge piu' nulla a mano)
    sim.setLaunches(LANCI_P, SAT_SEZIONE, SAT_MASSA);
    sim.saveInitialState();

    // Seed di base della run
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

        // Ripristino COMPLETO: gusci, tracker, tempo e RNG (seed di questa simulazione)
        sim.restoreInitialState(seedBase, static_cast<std::uint64_t>(r));

        for (int anno = 1; anno <= T_ANNI; ++anno) {

            sim.step(DT);

            // Conta gli anni instabili di ogni guscio (lambda_max > 1), in memoria
            for (int s = 1; s < numGusci; ++s) {   // il guscio 0 (atmosfera) non ha significato
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