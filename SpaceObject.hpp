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
    double crossSection;
    double mass;
    int initialShellIndex;
    
    ObjectType initialType;
    
    bool isManeuverable() const {
        return type == ObjectType::ACTIVE_PAYLOAD;
    }

    bool isTrackable() const {
        return crossSection > 0.01; 
    }
};