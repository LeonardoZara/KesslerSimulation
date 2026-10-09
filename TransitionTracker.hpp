#pragma once
 
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <limits>
#include "OrbitalShell.hpp"
 
class TransitionTracker {
private:
    int numShells;
 
    // F: una matrice 4x4 dei flussi assoluti PER OGNI GUSCIO (3 matrici separate).
    // Le righe (i) sono il tipo di destinazione, le colonne (j) il tipo di origine.
    // L'ordine dei tipi e' quello dell'enum ObjectType:
    //   0 = ACTIVE_PAYLOAD (P), 1 = DEAD_PAYLOAD (N), 2 = ROCKET_BODY (U), 3 = FRAGMENT (F)
    std::vector<Eigen::Matrix4d> flowMatrices;
 
    // X_k: popolazione iniziale di ogni tipo, per ogni guscio, all'inizio del tick
    std::vector<Eigen::Vector4d> initialPop;
 
public:
    // Costruttore: inizializza le matrici a zero con la grandezza giusta
    TransitionTracker(int shells)
        : numShells(shells),
          flowMatrices(shells, Eigen::Matrix4d::Zero()),
          initialPop(shells, Eigen::Vector4d::Zero()) {}
 
    // --- FASE 1: Snapshot Iniziale ---
    // Chiamato all'inizio di ogni anno simulato
    void reset(const std::vector<OrbitalShell>& systemShells) {
        for (int s = 0; s < numShells; ++s) {
            // Azzera la matrice dei flussi del guscio s
            flowMatrices[s].setZero();
 
            // Fotografa la popolazione iniziale X_k del guscio s, tipo per tipo
            initialPop[s].setZero();
            for (const auto& obj : systemShells[s].objects) {
                initialPop[s](static_cast<int>(obj.type)) += 1.0;
            }
        }
    }
 
    // --- FASE 2 & 3: Registrazione Eventi ---
    // Il motore fisico chiama questa funzione ogni volta che succede qualcosa
    // DENTRO il guscio 'shell'. Esempi:
    //   - un oggetto sopravvive:           recordFlow(s, tipo, tipo)
    //   - un frammento nasce da un P:      recordFlow(s, ACTIVE_PAYLOAD, FRAGMENT)
    void recordFlow(int shell, ObjectType fromType, ObjectType toType, double count = 1.0) {
        // Aggiunge 'count' alla casella (destinazione, origine) della matrice del guscio
        flowMatrices[shell](static_cast<int>(toType), static_cast<int>(fromType)) += count;
    }
 
    // --- FASE 4: La Magia dell'Algebra Lineare ---
    // Calcola e restituisce A_empirica (4x4) del guscio richiesto, a fine turno
    Eigen::Matrix4d computeEmpiricalMatrix(int shell) const {
        Eigen::Matrix4d A_emp = Eigen::Matrix4d::Zero();
 
        for (int j = 0; j < 4; ++j) {
            if (initialPop[shell](j) > 0.0) {
                // Estrae la colonna j dalla matrice dei flussi e la divide
                // per la popolazione iniziale di quel tipo in quel guscio
                A_emp.col(j) = flowMatrices[shell].col(j) / initialPop[shell](j);
            } else {
                // Se un tipo era completamente assente dal guscio all'inizio del turno,
                // il coefficiente "per unita' di x_j" non e' definito: lo segniamo con nan
                // (uno 0 qui sembrerebbe "tutti morti" e falserebbe la media nel tempo).
                A_emp.col(j).setConstant(std::numeric_limits<double>::quiet_NaN());
            }
        }
 
        return A_emp;
    }
 
    // (Opzionale) Restituisce X_k del guscio richiesto, utile per stampare i dati nel main
    Eigen::Vector4d getInitialPopulation(int shell) const {
        return initialPop[shell];
    }
};
 