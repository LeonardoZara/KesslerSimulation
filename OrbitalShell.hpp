#pragma once

#include <vector>
#include <array>
#include "SpaceObject.hpp"

class OrbitalShell {
public:
    int index;              // Indice topologico (0 = Atmosfera, 1 = LEO Bassa, ecc.)
    double baseAltitude;    // Quota di base del guscio in km (es. 500.0)
    double thickness;       // Spessore del guscio in km (es. 50.0, per il guscio 500-550 km)
    double volume;          // Volume totale in km^3 (V: serve per la densità spaziale)
    
    // Il cuore della struttura: un vettore dinamico che contiene tutti gli oggetti 
    // fisicamente presenti in questo guscio.
    std::vector<SpaceObject> objects;


    // Costruttore della classe. Inizializza i parametri spaziali del guscio.
    OrbitalShell(int idx, double alt, double thick, double vol) 
        : index(idx), baseAltitude(alt), thickness(thick), volume(vol) {
            // Un trucco di ottimizzazione C++: pre-allochiamo memoria per evitare 
            // riallocazioni continue quando aggiungiamo migliaia di detriti.
            // (10.000 è un numero arbitrario per dare spazio iniziale)
            objects.reserve(10000); 
        }

    // --- FUNZIONI FISICHE ---

    // Calcola la densità spaziale N/V in questo preciso istante.
    // Nota: oggetti.size() restituisce il numero N di oggetti attualmente nel vettore.
    double getDensity() const {
        if (volume == 0.0) return 0.0; // Prevenzione divisione per zero
        return objects.size() / volume;
    }

    // Pulisce il guscio svuotando il vettore (utile per i reset, se necessari)
    void clear() {
        objects.clear();
    }

    // Restituisce il numero di oggetti nel guscio senza esporre il vettore
    int getPopulationCount() const {
        return objects.size();
    }

    int getPopulationTypes(const std::vector<SpaceObject>& objects, int n){
        std::array<int, 4> counting ={0, 0, 0, 0};
        for(const auto& o : objects){
            counting[static_cast<size_t>(o.type)]++;
        }
        return counting[n];
    }
};