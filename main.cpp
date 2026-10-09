#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
 
#include "KesslerSimulation.hpp"
 
// Funzione helper per salvare la matrice Eigen in CSV
void salvaMatriceCSV(const std::string& nomeFile, const Eigen::MatrixXd& matrice) {
    std::ofstream file(nomeFile);
    if (file.is_open()) {
        const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
        file << matrice.format(CSVFormat);
        file.close();
    } else {
        std::cerr << "ERRORE: Impossibile aprire il file in scrittura: " << nomeFile << "\n"
                  << "Assicurati di aver creato la cartella 'output/' nella directory principale!" << '\n';
    }
}
 
void AnalisiSpettrale(int anni);
 
int main() {
    //4 gusci: 0 = Atmosfera, 1 = 300-400, 2 = 400-500, 3 = 500-600.
    int numGusci = 4;
    KesslerSimulation sim(numGusci);

    // Parametri: Guscio, Quantità, Tipo, Sezione d'Urto (m^2), Massa (kg)
    //DEVI INSERIRE MASSA E SEZIONE D'URTO CORRETTI
    //RICORDA DI METTERE UNA MATRICE PER GUSCIO
    sim.addObjectsToShell(1, 945, ObjectType::ACTIVE_PAYLOAD, 10.0, 500.0);
    sim.addObjectsToShell(1, 140, ObjectType::DEAD_PAYLOAD, 10.0, 500.0);
    sim.addObjectsToShell(1, 3000, ObjectType::FRAGMENT, 10.0, 1.0);

    sim.addObjectsToShell(2, 8520, ObjectType::ACTIVE_PAYLOAD, 15.0, 260.0);
    sim.addObjectsToShell(2, 1625, ObjectType::DEAD_PAYLOAD, 10.0, 500.0);
    sim.addObjectsToShell(2, 12000, ObjectType::FRAGMENT, 10.0, 1.0);

    sim.addObjectsToShell(3, 408, ObjectType::ACTIVE_PAYLOAD, 15.0, 260.0);
    sim.addObjectsToShell(3, 60, ObjectType::DEAD_PAYLOAD, 10.0, 500.0);
    sim.addObjectsToShell(3, 85000, ObjectType::FRAGMENT, 10.0, 1.0);

    // 3. ESECUZIONE DELLA SIMULAZIONE
    int anniDaSimulare = 50;
    std::cout << "Avvio simulazione Kessler per " << anniDaSimulare << " anni..." << '\n';
 
    for (int anno = 1; anno <= anniDaSimulare; ++anno) {
        
        // Fai avanzare il tempo di 1 anno
        sim.step(1.0);
        
        // Per ogni guscio: estrai la matrice 4x4 (tipo di oggetto -> tipo di oggetto)
        // calcolata in questo anno e salvala nel suo file CSV
        for (int s = 0; s < numGusci; ++s) {
            Eigen::Matrix4d A_empirica = sim.getTracker().computeEmpiricalMatrix(s);
 
            std::string nomeFile = "output/matrice_guscio_" + std::to_string(s)
                                 + "_anno_" + std::to_string(anno) + ".csv";
            salvaMatriceCSV(nomeFile, A_empirica);
        }
        
        // Feedback a terminale (stampa un aggiornamento ogni 10 anni)
        if (anno % 10 == 0 || anno == 1) {
            std::cout << "Simulato e salvato anno " << anno << " / " << anniDaSimulare << '\n';
        }
    }
 
    std::cout << "\nSimulazione completata." << '\n';
    std::cout << "\nAvvio dell'analisi dati e calcolo degli autovalori" << '\n';
    AnalisiSpettrale(anniDaSimulare);
 
    return 0; // Fine definitiva di tutto il programma
}