#pragma once

#include <vector>
#include <random>
#include <cstdint>
#include "OrbitalShell.hpp"
#include "TransitionTracker.hpp"

class KesslerSimulation {
private:
    int numShells;
    double currentTime;
    std::vector<OrbitalShell> shells;
    TransitionTracker tracker;
    std::mt19937 rng;
    std::vector<OrbitalShell> initialShells;

    // satelliti attivi lanciati per anno in ciascun guscio, e loro sezione d'urto / massa.
    std::vector<double> launchRates;
    double launchCrossSection;
    double launchMass;

public:
    // Costruttore
    KesslerSimulation(int numShells)
        : numShells(numShells),
          currentTime(0.0),
          tracker(numShells),
          rng(std::random_device{}()), // Inizializza il seed in modo puramente casuale dall'hardware (poi sovrascritto da restoreInitialState)
          launchRates(numShells, 0.0),
          launchCrossSection(37.5),
          launchMass(500.0)
    {
        // Creiamo FISICAMENTE i gusci in memoria
        shells.push_back(OrbitalShell(0, 0.0 , 300.0));
        for (int i = 1; i < numShells; ++i) {
            // Valori di esempio: Indice, Quota Base (km), Spessore (km), Volume (km^3)
            // Per la tesi potrai tarare il volume esatto in base alla quota
            shells.push_back(OrbitalShell(i, 300.0 + ((i - 1) * 100.0), 100.0));
        }
    }

    // Permette di caricare i dati (es. da un file CSV o simulando un lancio)
    void addObjectsToShell(int shellIndex, int count, ObjectType type, double crossSection, double mass);

    // [MODIFICA] Imposta i lanci della run: ratesPerYear ha un elemento per guscio (l'indice 0,
    // l'atmosfera, viene ignorato); crossSection (m^2) e mass (kg) valgono per i satelliti lanciati.
    void setLaunches(const std::vector<double>& ratesPerYear, double crossSection, double mass);

    // [MODIFICA] Salva la configurazione iniziale corrente (da chiamare UNA volta, dopo aver
    // aggiunto tutti gli oggetti con addObjectsToShell).
    void saveInitialState();

    // [MODIFICA] Riporta TUTTO allo stato iniziale (gusci, tracker, tempo) e re-inizializza l'RNG
    // in modo deterministico da (baseSeed, runIndex): ogni simulazione ha il suo flusso casuale.
    void restoreInitialState(std::uint64_t baseSeed, std::uint64_t runIndex);

    // Esegue UN singolo step temporale (es. dt = 1.0 anni)
    void step(double dt);

    // Restituisce il tracker per permettere al main di stampare le matrici
    const TransitionTracker& getTracker() const {
        return tracker;
    }

    // Restituisce l'anno attuale della simulazione
    double getCurrentTime() const {
        return currentTime;
    }

private:

    void processDecay(double dt);
    void processFailures(double dt);   // guasti P -> N (gamma e kappa)
    void processCollisions(double dt);
    void processLaunches(double dt);   // nuovi satelliti (forzante esterna g_k)
    void generateFragments(int shellIndex, const SpaceObject& obj1, const SpaceObject& obj2);
};