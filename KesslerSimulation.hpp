#pragma once

#include <vector>
#include <random>
#include <string>
#include "OrbitalShell.hpp"
#include "TransitionTracker.hpp"

class KesslerSimulation {
private:
    int numShells;
    double currentTime; // Tiene traccia degli anni passati
    
    // Le due componenti principali del nostro universo
    std::vector<OrbitalShell> shells;
    TransitionTracker tracker;

    std::mt19937 rng;

public:
    // Costruttore
    KesslerSimulation(int numShells) 
        : numShells(numShells), 
          currentTime(0.0), 
          tracker(numShells), 
          rng(std::random_device{}()) // Inizializza il seed in modo puramente casuale dall'hardware
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