#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <filesystem>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

static const double TOLLERANZA = 1e-4;

// Sotto questa popolazione una colonna é troppo rumorosa per dire se l'anno é instabile.
static const double POPOLAZIONE_MINIMA = 1000.0;

// Parte reale piu' grande tra gli autovalori. Le colonne non definite (nan: tipo assente a inizio anno) valgono 0.
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

// Un anno di una simulazione in UN guscio e' instabile se lambda_max della matrice empirica A di quell'anno e' > 1.
bool annoInstabile(const Eigen::Matrix4d& A4, const Eigen::Vector4d& X) {
    Eigen::MatrixXd A = A4;
    for (int j = 0; j < 4; ++j) {
        if (X(j) < POPOLAZIONE_MINIMA) A.col(j).setZero();
    }
    return lambdaMax(A) > 1.0 + TOLLERANZA;
}

// Riepilogo di una run (M simulazioni indipendenti con gli stessi lanci).
// anniInstabili[s][r] = numero di anni (su T) in cui il guscio s della simulazione r e' risultato instabile.
// Due numeri per guscio:
//  1. % di simulazioni instabili: una simulazione e' instabile se ha almeno kAnniMinimi anni instabili.
//     Margine di errore (regola pratica): +/- 100/sqrt(M) punti percentuali (prudente: e' il caso peggiore, p = 50%).
//  2. % di anni instabili: (somma degli anni instabili di tutte le simulazioni) / (M*T).
//     Margine: 1.96 * errore standard calcolato sulle M frazioni (una per simulazione), NON sugli M*T anni,
//     perche' gli anni di una stessa simulazione non sono indipendenti.
void AnalisiMonteCarlo(const std::vector<std::vector<int>>& anniInstabili,
                       int M, int T, int numGusci,
                       const std::vector<double>& lanci, std::uint64_t seedBase, int kAnniMinimi) {
    std::string nomiGusci[4] = {"Atmosfera", "300-400 km", "400-500 km", "500-600 km"}; // Deve coincidere con i gusci della simulazione

    if (M <= 0 || T <= 0) {
        std::cout << "Nessuna simulazione da analizzare." << '\n';
        return;
    }

    const double margineSimulazioni = 100.0 / std::sqrt(static_cast<double>(M));

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "\nRisultati:" << '\n';
    std::cout << "M = " << M << " simulazioni indipendenti, T = " << T << " anni, seed base = " << seedBase << '\n';
    std::cout << "Lanci annui di satelliti per guscio:";
    for (int s = 1; s < numGusci; ++s) std::cout << "  g" << s << " = " << lanci[s];
    std::cout << '\n';
    std::cout << "Anno instabile: lambda_max > 1 (anni con meno di " << POPOLAZIONE_MINIMA << " frammenti contano come stabili)" << '\n';
    std::cout << "Simulazione instabile: almeno " << kAnniMinimi << " anno/i instabile/i" << '\n';

    // Parte comune di ogni riga del CSV: lanci, M, T, k, seed
    std::ostringstream prefisso;
    prefisso << std::setprecision(10);
    for (int s = 1; s < numGusci; ++s) prefisso << lanci[s] << ",";
    prefisso << M << "," << T << "," << kAnniMinimi << "," << seedBase << ",";

    std::vector<std::string> righeCSV;

    for (int s = 1; s < numGusci; ++s) {
        int nInstabili = 0;
        std::vector<double> frazioni(M, 0.0);   // frazione di anni instabili di ogni simulazione
        for (int r = 0; r < M; ++r) {
            if (anniInstabili[s][r] >= kAnniMinimi) nInstabili++;
            frazioni[r] = static_cast<double>(anniInstabili[s][r]) / T;
        }

        double percSimulazioni = 100.0 * nInstabili / M;

        double media = 0.0;
        for (double f : frazioni) media += f;
        media /= M;
        double errStd = 0.0;
        if (M > 1) {
            double ss = 0.0;
            for (double f : frazioni) ss += (f - media) * (f - media);
            errStd = std::sqrt(ss / (M - 1)) / std::sqrt(static_cast<double>(M));
        }
        double percAnni = 100.0 * media;
        double margineAnni = 100.0 * 1.96 * errStd;

        std::cout << "\nGUSCIO " << s << " (" << nomiGusci[s] << ")" << '\n';
        std::cout << "  simulazioni instabili : " << nInstabili << " su " << M << " = " << percSimulazioni
                  << "%   (margine circa +/- " << margineSimulazioni << " punti)" << '\n';
        std::cout << "  anni instabili        : " << percAnni << "% degli anni simulati   (+/- " << margineAnni << ")" << '\n';

        std::ostringstream riga;
        riga << prefisso.str() << s << "," << std::fixed << std::setprecision(4)
             << nInstabili << "," << percSimulazioni << "," << margineSimulazioni << ","
             << percAnni << "," << margineAnni;
        righeCSV.push_back(riga.str());
    }

    // Accoda il riepilogo a output/risultati_run.csv (una riga per guscio per ogni esecuzione):
    // serve a fare i grafici "lanci vs % simulazioni instabili" e "lanci vs % anni instabili" mettendo insieme piu' run.
    const std::string percorso = "output/risultati_run.csv";
    bool scriviIntestazione = !std::filesystem::exists(percorso) || std::filesystem::file_size(percorso) == 0;
    std::ofstream out(percorso, std::ios::app);
    if (!out.is_open()) {
        std::cerr << "ERRORE: impossibile scrivere " << percorso << '\n';
        return;
    }
    if (scriviIntestazione) {
        for (int s = 1; s < numGusci; ++s) out << "lanci_g" << s << ",";
        out << "M,T,k_anni_minimi,seed_base,guscio,n_sim_instabili,perc_sim_instabili,margine_perc_sim,"
               "perc_anni_instabili,margine_perc_anni\n";
    }
    for (const auto& riga : righeCSV) out << riga << "\n";
    std::cout << "\nRiepilogo in " << percorso << '\n';
}