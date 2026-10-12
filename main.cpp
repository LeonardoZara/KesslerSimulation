#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
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
    std::filesystem::create_directories("output");
 
    //4 gusci: 0 = Atmosfera, 1 = 300-400, 2 = 400-500, 3 = 500-600.
    int numGusci = 4;
 
    KesslerSimulation sim(numGusci);
 
    // Fattore di scala della popolazione iniziale. 1.0 = dati del PDF.
    // Cambiandolo (0.2, 0.5, 2.0 ...) e rilanciando si vede dove lambda_max attraversa 1.
    double scala = 1.0;
    auto n = [scala](int quantità){ return static_cast<int>(quantità * scala);};
 
    // Parametri: Guscio, Quantità, Tipo, Sezione d'Urto (m^2), Massa (kg)
    // Satelliti: 130 m^2 e 800 kg. Frammenti: i valori comuni definiti in SpaceObject.hpp

    // 300-400 km
    sim.addObjectsToShell(1, n(945),   ObjectType::ACTIVE_PAYLOAD, 130.0, 800.0);
    sim.addObjectsToShell(1, n(140),   ObjectType::DEAD_PAYLOAD,   130.0, 800.0);
    sim.addObjectsToShell(1, n(3000),  ObjectType::FRAGMENT, fragmentSigma, fragmentMass);
 
    // 400-500 km
    // [MODIFICA] N0 = 1265 come nel PDF (nel tuo file c'era 1625, cifre scambiate).
    sim.addObjectsToShell(2, n(8520),  ObjectType::ACTIVE_PAYLOAD, 130.0, 800.0);
    sim.addObjectsToShell(2, n(1265),  ObjectType::DEAD_PAYLOAD,   130.0, 800.0);
    sim.addObjectsToShell(2, n(12000), ObjectType::FRAGMENT, fragmentSigma, fragmentMass);
 
    // 500-600 km
    sim.addObjectsToShell(3, n(408),   ObjectType::ACTIVE_PAYLOAD, 130.0, 800.0);
    sim.addObjectsToShell(3, n(60),    ObjectType::DEAD_PAYLOAD,   130.0, 800.0);
    sim.addObjectsToShell(3, n(85000), ObjectType::FRAGMENT, fragmentSigma, fragmentMass);
 
    int anniDaSimulare = 50;
    std::cout << "Avvio simulazione Kessler per " << anniDaSimulare << " anni..." << '\n';
 
    for (int anno = 1; anno <= anniDaSimulare; ++anno) {
 
        sim.step(1.0);
 
        // Per ogni guscio: estrae la matrice 4x4 (tipo di oggetto -> tipo di oggetto) calcolata in questo anno e salvala nel suo file CSV
        for (int s = 0; s < numGusci; ++s) {
            Eigen::Matrix4d A_empirica = sim.getTracker().computeEmpiricalMatrix(s);
 
            std::string nomeFile = "output/matrice_guscio_" + std::to_string(s)
                                 + "_anno_" + std::to_string(anno) + ".csv";
            salvaMatriceCSV(nomeFile, A_empirica);
 
            // Salva anche la popolazione iniziale X_k (P, N, U, F) di quell'anno.
            // L'analisi pesa ogni anno con X_k. Un anno in cui restano 10 frammenti
            // da' un coefficiente (flusso/10) enorme e rumoroso; con la media semplice un solo
            // anno cosi' bastava a far risultare "Kessler" un guscio dove i frammenti si estinguono.
            std::string nomePop = "output/popolazione_guscio_" + std::to_string(s)
                                + "_anno_" + std::to_string(anno) + ".csv";
            salvaMatriceCSV(nomePop, sim.getTracker().getInitialPopulation(s));
        }
 
        // Feedback a terminale
        if (anno % 10 == 0 || anno == 1) {
            std::cout << "Simulato e salvato anno " << anno << " / " << anniDaSimulare << '\n';
        }
    }
 
    std::cout << "\nSimulazione completata." << '\n';
    std::cout << "\nAvvio dell'analisi dati e calcolo degli autovalori" << '\n';
    AnalisiSpettrale(anniDaSimulare);
 
    return 0;
}