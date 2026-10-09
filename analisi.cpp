#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
 
// Funzione helper per leggere un CSV e trasformarlo in MatrixXd (dim x dim)
Eigen::MatrixXd leggiMatriceCSV(const std::string& nomeFile, int dim) {
    Eigen::MatrixXd matrice = Eigen::MatrixXd::Zero(dim, dim);
    std::ifstream file(nomeFile);
 
    if (!file.is_open()) {
        std::cerr << "Errore nella lettura del file: " << nomeFile << '\n';
        return matrice;
    }
 
    std::string linea, cella;
    int riga = 0;
 
    // Leggi il file riga per riga
    while (std::getline(file, linea) && riga < dim) {
        std::stringstream rigaStream(linea);
        int col = 0;
 
        // Separa i valori usando la virgola
        while (std::getline(rigaStream, cella, ',') && col < dim) {
            matrice(riga, col) = std::stod(cella); // stod converte string in double
            col++;
        }
        riga++;
    }
 
    file.close();
    return matrice;
}
 
void AnalisiSpettrale(int anni) {
    int numGusci = 4; // Deve coincidere con i gusci della tua simulazione
    int numTipi = 4;  // P, N, U, F (ordine dell'enum ObjectType): matrici 4x4
    std::string cartellaOutput = "output/";
    std::string nomiGusci[4] = {"Atmosfera", "300-400 km", "400-500 km", "500-600 km"};
 
    if (anni <= 0) {
        std::cout << "Nessun anno da analizzare." << '\n';
        return;
    }
 
    // 1. Tre matrici "somma" SEPARATE, una per guscio (niente matrice a blocchi).
    // Per ogni colonna (tipo di oggetto) contiamo anche in quanti anni era definita,
    // perche' la media va fatta solo sugli anni in cui quel tipo era presente.
    std::vector<Eigen::MatrixXd> matriciSomma(numGusci, Eigen::MatrixXd::Zero(numTipi, numTipi));
    std::vector<Eigen::VectorXd> anniValidi(numGusci, Eigen::VectorXd::Zero(numTipi));
 
    std::cout << "Inizio lettura dei file CSV in " << cartellaOutput << "..." << '\n';
 
    // 2. Per ogni guscio, leggiamo i CSV di tutti gli anni e li sommiamo
    for (int s = 0; s < numGusci; ++s) {
        for (int anno = 1; anno <= anni; ++anno) {
            // es. "output/matrice_guscio_2_anno_17.csv"
            std::string percorsoFile = cartellaOutput + "matrice_guscio_" + std::to_string(s)
                                     + "_anno_" + std::to_string(anno) + ".csv";
            Eigen::MatrixXd A = leggiMatriceCSV(percorsoFile, numTipi);
 
            for (int j = 0; j < numTipi; ++j) {
                if (std::isnan(A(0, j))) continue; // tipo assente a inizio anno: colonna non definita
                matriciSomma[s].col(j) += A.col(j);
                anniValidi[s](j) += 1.0;
            }
        }
    }
    std::cout << "Letti " << anni << " file per ciascuno dei " << numGusci << " gusci." << '\n';
 
    // 3. Per ogni guscio: media temporale, autovalori, verdetto
    std::vector<bool> kessler(numGusci, false);
    std::vector<std::string> esito(numGusci, "non applicabile");
 
    for (int s = 0; s < numGusci; ++s) {
        // Media temporale colonna per colonna (un tipo mai presente resta a 0)
        Eigen::MatrixXd matriceMedia = Eigen::MatrixXd::Zero(numTipi, numTipi);
        for (int j = 0; j < numTipi; ++j) {
            if (anniValidi[s](j) > 0.0) {
                matriceMedia.col(j) = matriciSomma[s].col(j) / anniValidi[s](j);
            }
        }
 
        std::cout << "\nGUSCIO " << s << " (" << nomiGusci[s] << ")" << '\n';
        std::cout << "MATRICE MEDIA (righe/colonne = P, N, U, F)" << '\n';
        std::cout << matriceMedia << '\n';
 
        // Il risolutore di Eigen per matrici generiche asimmetriche
        Eigen::EigenSolver<Eigen::MatrixXd> solver(matriceMedia);
 
        // Vettore di numeri complessi (std::complex<double>): gli autovalori
        // possono avere parte immaginaria.
        Eigen::VectorXcd autovalori = solver.eigenvalues();
 
        // L'autovalore dominante e' quello con la parte reale piu' grande
        double lambda_max = -1000.0;
 
        std::cout << "\nAutovalori:" << '\n';
        for (int i = 0; i < autovalori.size(); i++) {
            double parteReale = autovalori(i).real();
            double parteImmaginaria = autovalori(i).imag();
 
            std::cout << "  lambda_" << i << " = " << parteReale;
            if (std::abs(parteImmaginaria) > 1e-9) {
                std::cout << (parteImmaginaria > 0 ? " + " : " - ")
                          << std::abs(parteImmaginaria) << "i";
            }
            std::cout << '\n';
 
            if (parteReale > lambda_max) {
                lambda_max = parteReale;
            }
        }
        std::cout << "Autovalore dominante (lambda_max): " << lambda_max << '\n';
 
        // 4. Il Verdetto Fisico, guscio per guscio
        if (s == 0) {
            // Il guscio 0 e' l'atmosfera: stato assorbente, non puo' innescare Kessler
            std::cout << "Atmosfera (stato assorbente): verdetto non applicabile" << '\n';
        } else if (lambda_max > 1.0) {
            kessler[s] = true;
            esito[s] = "KESSLER (instabile)";
            std::cout << "Sindrome di Kessler innescata in questo guscio" << '\n';
        } else if (lambda_max < 1.0) {
            esito[s] = "stabile";
            std::cout << "Guscio stabile" << '\n';
        } else {
            esito[s] = "equilibrio critico";
            std::cout << "Equilibrio critico" << '\n';
        }
    }
 
    // 5. Riepilogo finale
    std::cout << "\nVerdetto del modello:" << '\n';
    bool almenoUno = false;
    for (int s = 1; s < numGusci; ++s) {
        std::cout << "Guscio " << s << " (" << nomiGusci[s] << "): " << esito[s] << '\n';
        if (kessler[s]) almenoUno = true;
    }
    if (almenoUno) {
        std::cout << "La popolazione di detriti cresce in modo esponenziale nei gusci instabili." << '\n';
    } else {
        std::cout << "Nessun guscio in regime di Kessler." << '\n';
    }
}