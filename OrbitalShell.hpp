#pragma once

#include <vector>
#include <array>
#include "SpaceObject.hpp"

class OrbitalShell {
public:
    int index;              // numero guscio
    double baseAltitude;    // Quota di base del guscio in km
    double thickness;       // Spessore del guscio in km 
    double volume;          // Volume totale in km^3
    
    std::vector<SpaceObject> objects;

    OrbitalShell(int idx, double alt, double thick) 
        : index(idx), baseAltitude(alt), thickness(thick) {
            double r_earth = 6371.0; // Raggio terrestre in km
            double r_out = r_earth + baseAltitude + thickness;
            double r_in = r_earth + baseAltitude;
            volume = (4.0 * 3.14159 / 3.0) * (std::pow(r_out, 3) - std::pow(r_in, 3));

            // Un trucco di ottimizzazione C++: pre-allochiamo memoria per evitare 
            // riallocazioni continue quando aggiungiamo migliaia di detriti.
            objects.reserve(110000); 
        }

    double getDensity() const {
        if (volume == 0.0) return 0.0;
        return objects.size() / volume;
    }

    void clear() {
        objects.clear();
    }

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