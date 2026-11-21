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
    std::vector<double> launchRates;
    double launchCrossSection;
    double launchMass;

public:
    KesslerSimulation(int numShells)
        : numShells(numShells),
          currentTime(0.0),
          tracker(numShells),
          rng(std::random_device{}()),
          launchRates(numShells, 0.0),
          launchCrossSection(37.5),
          launchMass(500.0)
    {
        shells.push_back(OrbitalShell(0, 0.0 , 300.0));
        for (int i = 1; i < numShells; ++i) {
            // Indice, Quota Base (km), Spessore (km), Volume (km^3)
            shells.push_back(OrbitalShell(i, 300.0 + ((i - 1) * 100.0), 100.0));
        }
    }

    void addObjectsToShell(int shellIndex, int count, ObjectType type, double crossSection, double mass);
    void setLaunches(const std::vector<double>& ratesPerYear, double crossSection, double mass);
    void saveInitialState();
    void restoreInitialState(std::uint64_t baseSeed, std::uint64_t runIndex);
    void step(double dt);
    const TransitionTracker& getTracker() const {
        return tracker;
    }
    double getCurrentTime() const {
        return currentTime;
    }

private:
    void processDecay(double dt);
    void processFailures(double dt);
    void processCollisions(double dt);
    void processLaunches(double dt);
    void generateFragments(int shellIndex, const SpaceObject& obj1, const SpaceObject& obj2);
};