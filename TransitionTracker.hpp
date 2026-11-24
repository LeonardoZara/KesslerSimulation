#pragma once

#include <Eigen/Dense>
#include <vector>
#include <limits>
#include "OrbitalShell.hpp"
 
class TransitionTracker {
private:
    int numShells;
 
    std::vector<Eigen::Matrix4d> flowMatrices;
    std::vector<Eigen::Vector4d> initialPop; // X_k: popolazione iniziale di ogni tipo, per ogni guscio, all'inizio del tick
 
public:
    TransitionTracker(int shells)
        : numShells(shells),
          flowMatrices(shells, Eigen::Matrix4d::Zero()),
          initialPop(shells, Eigen::Vector4d::Zero()) {}
 
    void reset(const std::vector<OrbitalShell>& systemShells) {
        for (int s = 0; s < numShells; ++s) {
            flowMatrices[s].setZero();
 
            // Fotografa la popolazione iniziale X_k del guscio s
            initialPop[s].setZero();
            for (const auto& obj : systemShells[s].objects) {
                initialPop[s](static_cast<int>(obj.type)) += 1.0;
            }
        }
    }
 
    void recordFlow(int shell, ObjectType fromType, ObjectType toType, double count = 1.0) {
        // Aggiunge 'count' alla casella (destinazione, origine) della matrice del guscio
        flowMatrices[shell](static_cast<int>(toType), static_cast<int>(fromType)) += count;
    }
 
    Eigen::Matrix4d computeEmpiricalMatrix(int shell) const {
        Eigen::Matrix4d A_emp = Eigen::Matrix4d::Zero();
 
        for (int j = 0; j < 4; ++j) {
            if (initialPop[shell](j) > 0.0) {
                A_emp.col(j) = flowMatrices[shell].col(j) / initialPop[shell](j);
            } else {
                // Se un tipo era completamente assente dal guscio all'inizio del turno, il coefficiente per unita' di x_j non e' definito: segnato con nan
                A_emp.col(j).setConstant(std::numeric_limits<double>::quiet_NaN());
            }
        }
        return A_emp;
    }
 
    Eigen::Vector4d getInitialPopulation(int shell) const {
        return initialPop[shell];
    }
};