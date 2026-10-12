#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
 
// Funzione helper per leggere un CSV e trasformarlo in MatrixXd (dim x dim)
bool leggiMatriceCSV(const std::string& nomeFile, int righe, int colonne, Eigen::MatrixXd& matrice) {
    matrice = Eigen::MatrixXd::Zero(righe, colonne);
    std::ifstream file(nomeFile);
 
    if (!file.is_open()) {
        std::cerr << "Errore nella lettura del file: " << nomeFile << '\n';
        return false;
    }
 
    std::string linea, cella;
    int riga = 0;
 
    while (std::getline(file, linea) && riga < righe) {
        std::stringstream rigaStream(linea);
        int col = 0;
 
        while (std::getline(rigaStream, cella, ',') && col < colonne) {
            matrice(riga, col) = std::stod(cella); // stod converte string in double
            col++;
        }
        riga++;
    }
 
    file.close();
    return true;
}
 
// Parte reale piu' grande tra gli autovalori. Le colonne non definite (nan: tipo assente a inizio anno) valgono 0. Serve per calcolare lambda_max anno per anno.
double lambdaMax(Eigen::MatrixXd A) {
    for (int i = 0; i < A.rows(); ++i)
        for (int j = 0; j < A.cols(); ++j)
            if (std::isnan(A(i, j))) A(i, j) = 0.0;
 
    Eigen::EigenSolver<Eigen::MatrixXd> solver(A);
    Eigen::VectorXcd autovalori = solver.eigenvalues();
 
    double massimo = -1000.0;
    for (int i = 0; i < autovalori.size(); i++) {
        if (autovalori(i).real() > massimo) massimo = autovalori(i).real();
    }
    return massimo;
}
 
void AnalisiSpettrale(int anni) {
    int numGusci = 4; // Deve coincidere con i gusci della simulazione
    int numTipi = 4;  // P, N, U, F: matrici 4x4
    std::string cartellaOutput = "output/";
    std::string nomiGusci[4] = {"Atmosfera", "300-400 km", "400-500 km", "500-600 km"};
 
    // Tolleranza sul confronto con 1.
    const double TOLLERANZA = 1e-4;
 
    if (anni <= 0) {
        std::cout << "Nessun anno da analizzare." << '\n';
        return;
    }
 
    // 1. Tre matrici "somma" SEPARATE, una per guscio (niente matrice a blocchi).
    // Media pesata con la popolazione iniziale invece della media semplice.
    // Ogni colonna j di ogni anno entra con peso X_k(j), cioé il numero di oggetti di quel tipo a inizio anno. Equivale a (somma dei flussi) / (somma delle popolazioni).
    // a_ij = flusso / X_j; se X_j è piccolo il rapporto è enorme e casuale (10 frammenti e un urto: a_44 = 28).
    // Un tipo assente (X_j = 0) ha peso 0, quindi sparisce da solo.
    std::vector<Eigen::MatrixXd> matriciSomma(numGusci, Eigen::MatrixXd::Zero(numTipi, numTipi));
    std::vector<Eigen::VectorXd> anniValidi(numGusci, Eigen::VectorXd::Zero(numTipi)); // somma dei pesi X_k(j)
 
    // Sotto questa popolazione una colonna é troppo rumorosa per dire se l'anno é instabile (usata solo per il conteggio degli anni con lambda_max > 1). 
    // Un singolo urto catastrofico genera ~250 frammenti per genitore, quindi con X oggetti sposta a_44 di 250/X: con 1000 oggetti al massimo di 0.25, con 100 di 2.5 (basta un urto per superare 1).
    const double POPOLAZIONE_MINIMA = 1000.0;
 
    // lambda_max di ogni singolo anno, per ogni guscio.
    // Dopo una cascata i bersagli sono finiti, le matrici degli ultimi anni tornano stabili e la media su 50 anni puo' scendere sotto 1 anche se la sindrome é scoppiata nei primi anni.
    std::vector<int> anniInstabili(numGusci, 0);   // anni con lambda_max > 1
    std::vector<int> primoAnnoInstabile(numGusci, 0);
    std::vector<double> lambdaPicco(numGusci, 0.0);
 
    std::cout << "Inizio lettura dei file CSV in " << cartellaOutput << "..." << '\n';
 
    for (int s = 0; s < numGusci; ++s) {
        for (int anno = 1; anno <= anni; ++anno) {
            std::string percorsoFile = cartellaOutput + "matrice_guscio_" + std::to_string(s)
                                     + "_anno_" + std::to_string(anno) + ".csv";
            Eigen::MatrixXd A;
            if (!leggiMatriceCSV(percorsoFile, numTipi, numTipi, A)) continue;
 
            std::string percorsoPop = cartellaOutput + "popolazione_guscio_" + std::to_string(s)
                                    + "_anno_" + std::to_string(anno) + ".csv";
            Eigen::MatrixXd X;
            if (!leggiMatriceCSV(percorsoPop, numTipi, 1, X)) continue;
 
            for (int j = 0; j < numTipi; ++j) {
                if (std::isnan(A(0, j)) || X(j, 0) <= 0.0) continue; // tipo assente: colonna non definita
                matriciSomma[s].col(j) += X(j, 0) * A.col(j); // peso = X_k(j)
                anniValidi[s](j) += X(j, 0);
            }
 
            // lambda_max di quest'anno (per il guscio 0 non ha significato: salta).
            // Le colonne con meno di popolazione minima oggetti valgono 0 (troppo rumorose).
            if (s > 0) {
                for (int j = 0; j < numTipi; ++j) {
                    if (X(j, 0) < POPOLAZIONE_MINIMA) A.col(j).setZero();
                }
                double lam = lambdaMax(A);
                if (lam > lambdaPicco[s]) lambdaPicco[s] = lam;
                if (lam > 1.0 + TOLLERANZA) {
                    anniInstabili[s]++;
                    if (primoAnnoInstabile[s] == 0) primoAnnoInstabile[s] = anno;
                }
            }
        }
    }
    std::cout << "Letti " << anni << " file per ciascuno dei " << numGusci << " gusci." << '\n';
 
    std::vector<bool> kessler(numGusci, false);
    std::vector<std::string> esito(numGusci, "non applicabile");
 
    for (int s = 0; s < numGusci; ++s) {
        // Media temporale pesata colonna per colonna (un tipo mai presente resta a 0)
        Eigen::MatrixXd matriceMedia = Eigen::MatrixXd::Zero(numTipi, numTipi);
        for (int j = 0; j < numTipi; ++j) {
            if (anniValidi[s](j) > 0.0) {
                matriceMedia.col(j) = matriciSomma[s].col(j) / anniValidi[s](j);
            }
        }
 
        std::cout << "\nGUSCIO " << s << " (" << nomiGusci[s] << ")" << '\n';
        std::cout << "MATRICE MEDIA (righe/colonne = P, N, U, F)" << '\n';
        std::cout << matriceMedia << '\n';
 
        Eigen::EigenSolver<Eigen::MatrixXd> solver(matriceMedia);
 
        // Vettore di numeri complessi (std::complex<double>): gli autovalori possono avere parte immaginaria.
        Eigen::VectorXcd autovalori = solver.eigenvalues();
 
        // L'autovalore dominante é quello con la parte reale piu' grande
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
 
        if (s == 0) {
            std::cout << "Atmosfera (stato assorbente): verdetto non applicabile" << '\n';
            continue;
        }
 
        std::cout << "lambda_max anno per anno: picco " << lambdaPicco[s] << ", maggiore di 1 in "
                  << anniInstabili[s] << " anni su " << anni;
        if (anniInstabili[s] > 0) std::cout << " (primo: anno " << primoAnnoInstabile[s] << ")";
        std::cout << '\n';
 
        if (lambda_max > 1.0 + TOLLERANZA) {
            kessler[s] = true;
            esito[s] = "KESSLER (instabile)";
            std::cout << "Sindrome di Kessler innescata in questo guscio" << '\n';
        } else if (anniInstabili[s] > 0) {
            // La media é < 1 ma in alcuni anni lambda_max era > 1: la cascata é partita e poi si é fermata (bersagli esauriti o frammenti decaduti).
            kessler[s] = true;
            esito[s] = "KESSLER TRANSITORIA (instabile nei primi anni, poi esaurita)";
            std::cout << "Media stabile, ma cascata innescata in " << anniInstabili[s] << " anni" << '\n';
        } else if (lambda_max < 1.0 - TOLLERANZA) {
            esito[s] = "stabile";
            std::cout << "Guscio stabile" << '\n';
        } else {
            esito[s] = "equilibrio critico";
            std::cout << "Equilibrio critico" << '\n';
        }
    }
 
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