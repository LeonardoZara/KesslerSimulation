#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>

#include "KesslerSimulation.hpp"

// Funzione helper per salvare la matrice Eigen in CSV
void salvaMatriceCSV(const std::string& nomeFile, const Eigen::MatrixXd& matrice) {
    std::ofstream file(nomeFile);
    if (file.is_open()) {
        const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
        file << matrice.format(CSVFormat);
        file.close();
    } else {
        std::cerr << "ERRORE: Impossibile aprire il file in scrittura: " << nomeFile << "\n"
                  << "Assicurati di aver creato la cartella 'output/' nella directory principale!" << '\n';
    }
}

void AnalisiSpettrale();

int main() {
    // 1. INIZIALIZZAZIONE DELLO SPAZIO
    // Creiamo 4 gusci: 0 = Atmosfera, 1 = 300-400, 2 = 400-500, 3 = 500-600.
    int numGusci = 4;
    KesslerSimulation sim(numGusci);

    // 2. POPOLAMENTO INIZIALE (Qui puoi sbizzarrirti per la tesi)
    // Parametri: Guscio, Quantità, Tipo, Sezione d'Urto (m^2), Massa (kg)
    //DEVI INSERIRE MASSA E SEZIONE D'URTO CORRETTI
    //RICORDA DI METTERE UNA MATRICE PER GUSCIO
    
    sim.addObjectsToShell(1, 945, ObjectType::ACTIVE_PAYLOAD, 10.0, 500.0);
    sim.addObjectsToShell(1, 140, ObjectType::DEAD_PAYLOAD, 10.0, 500.0);
    sim.addObjectsToShell(1, 3000, ObjectType::FRAGMENT, 10.0, 500.0);


    sim.addObjectsToShell(2, 8520, ObjectType::ACTIVE_PAYLOAD, 15.0, 260.0);
    sim.addObjectsToShell(2, 1625, ObjectType::DEAD_PAYLOAD, 10.0, 500.0);
    sim.addObjectsToShell(2, 12000, ObjectType::FRAGMENT, 10.0, 500.0);

    sim.addObjectsToShell(3, 408, ObjectType::ACTIVE_PAYLOAD, 15.0, 260.0);
    sim.addObjectsToShell(3, 60, ObjectType::DEAD_PAYLOAD, 10.0, 500.0);
    sim.addObjectsToShell(3, 85000, ObjectType::FRAGMENT, 10.0, 500.0);

    // 3. ESECUZIONE DELLA SIMULAZIONE
    int anniDaSimulare = 50;
    std::cout << "Avvio simulazione Kessler per " << anniDaSimulare << " anni..." << '\n';

    for (int anno = 1; anno <= anniDaSimulare; ++anno) {
        
        // Fai avanzare il tempo di 1 anno
        sim.step(1.0);
        
        // Estrai la matrice dei flussi empirici calcolata in questo anno
        Eigen::MatrixXd A_empirica = sim.getTracker().computeEmpiricalMatrix();

        // Salva la matrice nel file CSV
        std::string nomeFile = "output/matrice_empirica_anno_" + std::to_string(anno) + ".csv";
        salvaMatriceCSV(nomeFile, A_empirica);
        
        // Feedback a terminale (stampa un aggiornamento ogni 10 anni)
        if (anno % 10 == 0 || anno == 1) {
            std::cout << "Simulato e salvato anno " << anno << " / " << anniDaSimulare << '\n';
        }
    }

    std::cout << "\nSimulazione completata con successo!" << '\n';
    std::cout << "\nAvvio dell'analisi dati e calcolo degli autovalori." << '\n';
    AnalisiSpettrale(); 

    return 0;
}