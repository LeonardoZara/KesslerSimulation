#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace fs = std::filesystem;

// Funzione helper per leggere un CSV e trasformarlo in MatrixXd
Eigen::MatrixXd leggiMatriceCSV(const std::string& nomeFile, int numGusci) {
    Eigen::MatrixXd matrice = Eigen::MatrixXd::Zero(numGusci, numGusci);
    std::ifstream file(nomeFile);
    
    if (!file.is_open()) {
        std::cerr << "Errore nella lettura del file: " << nomeFile << '\n';
        return matrice;
    }

    std::string linea, cella;
    int riga = 0;
    
    // Leggi il file riga per riga
    while (std::getline(file, linea) && riga < numGusci) {
        std::stringstream rigaStream(linea);
        int col = 0;
        
        // Separa i valori usando la virgola
        while (std::getline(rigaStream, cella, ',') && col < numGusci) {
            matrice(riga, col) = std::stod(cella); // stod converte string in double
            col++;
        }
        riga++;
    }
    
    file.close();
    return matrice;
}

void AnalisiSpettrale() {
    int numGusci = 3; // Deve coincidere con i gusci della tua simulazione
    std::string cartellaOutput = "output/";

    // 1. Prepariamo una matrice vuota per accumulare la somma
    Eigen::MatrixXd matriceSomma = Eigen::MatrixXd::Zero(numGusci, numGusci);
    int numeroFileLetti = 0;

    std::cout << "Inizio lettura dei file CSV in " << cartellaOutput << "." << '\n';

    // 2. Esploriamo la cartella!
    for (const auto& entry : fs::directory_iterator(cartellaOutput)) {
        
        // Controlliamo che sia un file normale e che finisca per .csv
        if (entry.is_regular_file() && entry.path().extension() == ".csv") {
            
            // Estraiamo il percorso completo (es. "output/matrice_empirica_anno_1.csv")
            std::string percorsoFile = entry.path().string();
            
            // Leggiamo la matrice
            Eigen::MatrixXd A = leggiMatriceCSV(percorsoFile, numGusci);
            
            // La sommiamo al totale
            matriceSomma += A;
            numeroFileLetti++;
        }
    }

    // 3. Calcolo e Salvataggio della Media Finale
    if (numeroFileLetti > 0) {
        Eigen::MatrixXd matriceMedia = matriceSomma / numeroFileLetti;
        
        std::cout << "Letti con successo " << numeroFileLetti << " file." << '\n';
        std::cout << "\n Matrice media: " << '\n';
        std::cout << matriceMedia << '\n';
        
        // (Opzionale) Puoi salvare questa matrice media in un file a parte
        // salvaMatriceCSV("output/MATRICE_FINALE_MEDIA.csv", matriceMedia);

        // 1. Inizializziamo il risolutore di Eigen per matrici generiche asimmetriche
        Eigen::EigenSolver<Eigen::MatrixXd> solver(matriceMedia);

        // 2. Eigen restituisce un vettore di numeri complessi (std::complex<double>).
        // Questo perché matematicamente gli autovalori possono avere parte immaginaria.
        Eigen::VectorXcd autovalori = solver.eigenvalues();

        // 3. Troviamo l'autovalore con la parte reale più grande (il dominante)
        double lambda_max = -1000.0;
            
        std::cout << "Lista degli autovalori trovati:" << '\n';
        for (int i = 0; i < autovalori.size(); i++) {
            double parteReale = autovalori(i).real();
            double parteImmaginaria = autovalori(i).imag();
                
            std::cout << "  lambda_" << i << " = " << parteReale;
            // Stampiamo la parte immaginaria solo se esiste (per eleganza)
            if (std::abs(parteImmaginaria) > 1e-9) {
                std::cout << (parteImmaginaria > 0 ? " + " : " - ") 
                            << std::abs(parteImmaginaria) << "i";
            }
            std::cout << '\n';

            // Aggiorniamo il massimo
            if (parteReale > lambda_max) {
                lambda_max = parteReale;
            }
        }

        std::cout << "\nAutovalore dominante (lambda_max): " << lambda_max << '\n';

        // 4. Il Verdetto Fisico
        std::cout << "\nVerdetto della simulazione: ";
        if (lambda_max > 1.0) {
            std::cout << "sindrome di Kessler innescata!" << '\n';
            std::cout << "La popolazione dei detriti crescera' in modo esponenziale." << '\n';
        } else if (lambda_max < 1.0) {
            std::cout << "ambiente spaziale stabile!" << '\n';
            std::cout << "I detriti decadranno fisiologicamente nel tempo." << '\n';
        } else {
            std::cout << "equilibrio critico." << '\n';
        }

    } else {
        std::cout << "Nessun file CSV trovato nella cartella." << '\n';
    }
}