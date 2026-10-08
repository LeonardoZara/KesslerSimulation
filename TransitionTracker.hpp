#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "OrbitalShell.hpp"

class TransitionTracker {
private:
    int numShells;
    
    // F: La matrice dei flussi assoluti (chi va dove)
    // Le righe (i) sono la destinazione, le colonne (j) sono l'origine.
    Eigen::MatrixXd flowMatrix;

    // X_k: Il vettore di stato all'inizio del tick
    Eigen::VectorXd initialPop;

    Eigen::VectorXd initialShellPop;
public:
    //Matrici di transizione dei ogni guscio guscio
    std::vector<Eigen::Matrix4d> shellMatrixes;
    int i=4;
    // Costruttore: inizializza le matrici a zero con la grandezza giusta
    TransitionTracker(int shells) : numShells(shells), shellMatrixes(shells, Eigen::MatrixXd::Zero()) {
        flowMatrix = Eigen::MatrixXd::Zero(shells, shells);
        initialPop = Eigen::VectorXd::Zero(shells);
    }

    // --- FASE 1: Snapshot Iniziale ---
    // Chiamato all'inizio di ogni anno simulato
    void reset(const std::vector<OrbitalShell>& systemShells) {
        // Azzera la matrice dei flussi
        flowMatrix.setZero();
        
        // Fotografa la popolazione iniziale X_k
        for (int i = 0; i < numShells; ++i) {
            initialPop(i) = systemShells[i].getPopulationCount();
        }

    }

    // --- FASE 2 & 3: Registrazione Eventi ---
    // Il motore fisico chiama questa funzione ogni volta che succede qualcosa
    // count = 1 di default (es. un satellite decade). Ma se c'è un urto 
    // e si creano 200 frammenti, chiamerai recordFlow(guscio, guscio, 200).
    void recordFlow(int fromShell, int toShell, double count = 1.0) {
        // Aggiunge 'count' alla casella corrispondente
        flowMatrix(toShell, fromShell) += count;
    }

    // --- FASE 4: La Magia dell'Algebra Lineare ---
    // Calcola e restituisce A_empirica a fine turno
    Eigen::MatrixXd computeEmpiricalMatrix() const {
        Eigen::MatrixXd A_emp = Eigen::MatrixXd::Zero(numShells, numShells);
        
        for (int j = 0; j < numShells; ++j) {
            // Prevenzione divisione per zero: se un guscio era completamente 
            // vuoto all'inizio del turno, la sua colonna in A resta a zero.
            if (initialPop(j) > 0.0) {
                // Estrae la colonna j dalla matrice dei flussi e la divide 
                // per la popolazione iniziale di quel guscio
                A_emp.col(j) = flowMatrix.col(j) / initialPop(j);
            }
        }
        
        return A_emp;
    }

    // (Opzionale) Restituisce X_k, utile per stampare i dati nel main
    Eigen::VectorXd getInitialPopulation() const {
        return initialPop;
    }
};