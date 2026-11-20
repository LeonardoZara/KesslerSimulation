#pragma once

enum class ObjectType {
    ACTIVE_PAYLOAD=0, 
    DEAD_PAYLOAD=1,   
    ROCKET_BODY=2,    
    FRAGMENT=3
};

constexpr double fragmentSigma = 0.00785; // m^2 (sotto 0.01: non tracciabile, non schivabile)
constexpr double fragmentMass  = 1.0;   // kg

struct SpaceObject {
    ObjectType type;
    double crossSection;    // Area in m^2 (sigma: per calcolare la probabilità di urto)
    double mass;            // Massa in kg (per calcolare il numero di frammenti post-urto)
    int initialShellIndex;  // Salva l'indice del guscio in cui l'oggetto si trovava all'inizio dell'anno simulato.
    
    ObjectType initialType;
    
    bool isManeuverable() const {
        return type == ObjectType::ACTIVE_PAYLOAD;
    }

    bool isTrackable() const {
        return crossSection > 0.01; 
    }
};