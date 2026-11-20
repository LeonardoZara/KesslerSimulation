#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <filesystem>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Tolleranza sul confronto con 1.
static const double TOLLERANZA = 1e-4;

// Sotto questa popolazione una colonna é troppo rumorosa per dire se l'anno é instabile.
// Un singolo urto catastrofico genera ~250 frammenti per genitore, quindi con X oggetti sposta a_44 di 250/X: con 1000 oggetti al massimo di 0.25, con 100 di 2.5 (basta un urto per superare 1).
static const double POPOLAZIONE_MINIMA = 1000.0;

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

// lambda_max della matrice empirica di UN guscio in UN anno di UNA simulazione.
// X = popolazione iniziale dell'anno (P, N, U, F). Le colonne con X < POPOLAZIONE_MINIMA valgono 0 (stessa regola di prima).
// valutabile = true se la popolazione di frammenti F era >= POPOLAZIONE_MINIMA. Con la struttura attuale la matrice e'
// triangolare bassa e P, N, U hanno a_ii <= 1: un anno puo' avere lambda_max > 1 solo attraverso a_44 (frammenti), quindi
// un anno con F sotto soglia non e' valutabile (non puo' risultare instabile) e non entra nelle statistiche di lambda_max.
double lambdaAnnuale(const Eigen::Matrix4d& A4, const Eigen::Vector4d& X, bool& valutabile) {
    Eigen::MatrixXd A = A4;
    for (int j = 0; j < 4; ++j) {
        if (X(j) < POPOLAZIONE_MINIMA) A.col(j).setZero();
    }
    valutabile = (X(3) >= POPOLAZIONE_MINIMA);
    return lambdaMax(A);
}

// ---------------------------------------------------------------------------------------------
// Statistica di una proporzione k/M (M run indipendenti): errore standard e intervallo di Wilson al 95%
// (preferito a media +/- 1.96*SE perche' funziona anche per p vicino a 0 o 1, e per k = 0 o k = M).
struct StatBinomiale { double p, se, lo, hi; };

static StatBinomiale statBinomiale(int k, int M) {
    const double z = 1.959964;
    StatBinomiale r;
    r.p  = static_cast<double>(k) / M;
    r.se = std::sqrt(r.p * (1.0 - r.p) / M);
    double z2 = z * z;
    double denom  = 1.0 + z2 / M;
    double centro = (r.p + z2 / (2.0 * M)) / denom;
    double semi   = z * std::sqrt(r.p * (1.0 - r.p) / M + z2 / (4.0 * M * static_cast<double>(M))) / denom;
    r.lo = std::max(0.0, centro - semi);
    r.hi = std::min(1.0, centro + semi);
    return r;
}

// Statistiche descrittive di un campione (una cifra per ogni run).
struct StatDescrittive { int n; double media, sd, se, minimo, mediana, massimo; };

static StatDescrittive statDescrittive(std::vector<double> v) {
    StatDescrittive s{0, NAN, NAN, NAN, NAN, NAN, NAN};
    s.n = static_cast<int>(v.size());
    if (s.n == 0) return s;
    std::sort(v.begin(), v.end());
    double somma = 0.0;
    for (double x : v) somma += x;
    s.media = somma / s.n;
    double ss = 0.0;
    for (double x : v) ss += (x - s.media) * (x - s.media);
    s.sd = (s.n > 1) ? std::sqrt(ss / (s.n - 1)) : 0.0;
    s.se = s.sd / std::sqrt(static_cast<double>(s.n));
    s.minimo  = v.front();
    s.massimo = v.back();
    s.mediana = (s.n % 2 == 1) ? v[s.n / 2] : 0.5 * (v[s.n / 2 - 1] + v[s.n / 2]);
    return s;
}

// Analisi finale. lambda[r][s][t] e valutabile[r][s][t]: simulazione r (0..M-1), guscio s, anno t+1 (t = 0..T-1).
// L'unita' statistica e' la SIMULAZIONE (M campioni indipendenti), non l'anno: gli anni di una stessa
// simulazione sono correlati tra loro, quindi M*T non e' il numero di campioni indipendenti.
void AnalisiMonteCarlo(const std::vector<std::vector<std::vector<double>>>& lambda,
                       const std::vector<std::vector<std::vector<char>>>& valutabile,
                       int M, int T, int numGusci,
                       const std::vector<double>& lanci, std::uint64_t seedBase) {
    std::string nomiGusci[4] = {"Atmosfera", "300-400 km", "400-500 km", "500-600 km"}; // Deve coincidere con i gusci della simulazione

    if (M <= 0 || T <= 0) {
        std::cout << "Nessuna simulazione da analizzare." << '\n';
        return;
    }

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n================ ANALISI MONTE CARLO ================" << '\n';
    std::cout << "M = " << M << " simulazioni indipendenti, T = " << T << " anni, seed base = " << seedBase << '\n';
    std::cout << "Lanci annui di satelliti per guscio:";
    for (int s = 1; s < numGusci; ++s) std::cout << "  g" << s << " = " << lanci[s];
    std::cout << '\n';
    std::cout << "Criterio: lambda_max(A_t) > 1 + " << TOLLERANZA << " in un anno t con almeno " << POPOLAZIONE_MINIMA
              << " frammenti a inizio anno (anni con meno frammenti: non valutabili)" << '\n';

    std::vector<char> instabileQualsiasi(M, 0);   // per run: instabile in almeno un guscio e un anno
    std::vector<char> finaleQualsiasi(M, 0);      // per run: instabile in almeno un guscio all'anno T

    // Riga di riepilogo per ogni guscio, da accodare al CSV
    std::vector<std::string> righeCSV;
    std::ostringstream prefissoLanci;
    prefissoLanci << std::setprecision(10);
    for (int s = 1; s < numGusci; ++s) prefissoLanci << lanci[s] << ",";
    prefissoLanci << M << "," << T << "," << seedBase << ",";

    for (int s = 1; s < numGusci; ++s) {
        int kAlmenoUna = 0, kFinale = 0;
        std::vector<double> picchi, medie, primiAnni;
        double sommaAnniVal = 0.0;

        for (int r = 0; r < M; ++r) {
            int nVal = 0;
            double picco = -1.0e300, somma = 0.0;
            int primo = 0;
            for (int t = 0; t < T; ++t) {
                if (!valutabile[r][s][t]) continue;
                double lam = lambda[r][s][t];
                nVal++;
                somma += lam;
                if (lam > picco) picco = lam;
                if (lam > 1.0 + TOLLERANZA && primo == 0) primo = t + 1;
            }
            sommaAnniVal += nVal;

            bool inst = (primo > 0);
            bool fin  = valutabile[r][s][T - 1] && (lambda[r][s][T - 1] > 1.0 + TOLLERANZA);
            if (inst) { kAlmenoUna++; primiAnni.push_back(primo); instabileQualsiasi[r] = 1; }
            if (fin)  { kFinale++;    finaleQualsiasi[r] = 1; }
            if (nVal > 0) { picchi.push_back(picco); medie.push_back(somma / nVal); }
        }

        StatBinomiale bAny = statBinomiale(kAlmenoUna, M);
        StatBinomiale bFin = statBinomiale(kFinale, M);
        StatDescrittive dPicco = statDescrittive(picchi);
        StatDescrittive dMedia = statDescrittive(medie);
        StatDescrittive dPrimo = statDescrittive(primiAnni);
        double anniValMedi = sommaAnniVal / M;

        std::cout << "\nGUSCIO " << s << " (" << nomiGusci[s] << ")" << '\n';
        std::cout << "  lambda_max > 1 in almeno un anno : " << kAlmenoUna << "/" << M << " = " << bAny.p
                  << "   SE = " << bAny.se << "   IC95% Wilson [" << bAny.lo << ", " << bAny.hi << "]" << '\n';
        std::cout << "  lambda_max > 1 all'anno finale T : " << kFinale << "/" << M << " = " << bFin.p
                  << "   SE = " << bFin.se << "   IC95% Wilson [" << bFin.lo << ", " << bFin.hi << "]" << '\n';
        if (kAlmenoUna > 0) {
            std::cout << "  primo anno con lambda_max > 1 (mediana tra i " << kAlmenoUna << " run instabili): " << dPrimo.mediana << '\n';
        }
        std::cout << "  anni valutabili per simulazione (media): " << anniValMedi << " su " << T << '\n';
        std::cout << "  lambda_max di picco (max sugli anni valutabili, " << dPicco.n << " run): media = " << dPicco.media
                  << " +/- " << dPicco.se << " (SE), SD = " << dPicco.sd << ", min = " << dPicco.minimo
                  << ", mediana = " << dPicco.mediana << ", max = " << dPicco.massimo << '\n';
        std::cout << "  lambda_max medio (media sugli anni valutabili, " << dMedia.n << " run): media = " << dMedia.media
                  << " +/- " << dMedia.se << " (SE), SD = " << dMedia.sd << ", min = " << dMedia.minimo
                  << ", mediana = " << dMedia.mediana << ", max = " << dMedia.massimo << '\n';

        std::ostringstream riga;
        riga << prefissoLanci.str() << s << "," << std::setprecision(6)
             << kAlmenoUna << "," << bAny.p << "," << bAny.se << "," << bAny.lo << "," << bAny.hi << ","
             << kFinale << "," << bFin.p << "," << bFin.lo << "," << bFin.hi << ","
             << dPicco.n << "," << dPicco.media << "," << dPicco.se << "," << dPicco.sd << ","
             << dPicco.minimo << "," << dPicco.mediana << "," << dPicco.massimo << ","
             << dMedia.media << "," << anniValMedi;
        righeCSV.push_back(riga.str());
    }

    // Riga "qualsiasi guscio": evento = in quella simulazione almeno un guscio e' instabile
    int kAnyShell = 0, kFinShell = 0;
    for (int r = 0; r < M; ++r) { kAnyShell += instabileQualsiasi[r]; kFinShell += finaleQualsiasi[r]; }
    StatBinomiale bAnyShell = statBinomiale(kAnyShell, M);
    StatBinomiale bFinShell = statBinomiale(kFinShell, M);
    std::cout << "\nQUALSIASI GUSCIO (almeno un guscio instabile nella stessa simulazione)" << '\n';
    std::cout << "  in almeno un anno : " << kAnyShell << "/" << M << " = " << bAnyShell.p
              << "   SE = " << bAnyShell.se << "   IC95% Wilson [" << bAnyShell.lo << ", " << bAnyShell.hi << "]" << '\n';
    std::cout << "  all'anno finale T : " << kFinShell << "/" << M << " = " << bFinShell.p
              << "   SE = " << bFinShell.se << "   IC95% Wilson [" << bFinShell.lo << ", " << bFinShell.hi << "]" << '\n';
    {
        std::ostringstream riga;
        riga << prefissoLanci.str() << "qualsiasi," << std::setprecision(6)
             << kAnyShell << "," << bAnyShell.p << "," << bAnyShell.se << "," << bAnyShell.lo << "," << bAnyShell.hi << ","
             << kFinShell << "," << bFinShell.p << "," << bFinShell.lo << "," << bFinShell.hi << ","
             << "0,nan,nan,nan,nan,nan,nan,nan,nan";
        righeCSV.push_back(riga.str());
    }

    // Accoda il riepilogo a output/riepilogo_run.csv (una riga per guscio per ogni esecuzione):
    // serve a costruire il grafico "lanci vs frazione di simulazioni instabili" mettendo insieme piu' run.
    const std::string percorso = "output/riepilogo_run.csv";
    bool scriviIntestazione = !std::filesystem::exists(percorso) || std::filesystem::file_size(percorso) == 0;
    std::ofstream out(percorso, std::ios::app);
    if (!out.is_open()) {
        std::cerr << "ERRORE: impossibile scrivere " << percorso << '\n';
        return;
    }
    if (scriviIntestazione) {
        for (int s = 1; s < numGusci; ++s) out << "lanci_g" << s << ",";
        out << "M,T,seed_base,guscio,k_almeno_un_anno,p_almeno_un_anno,se,ic95_basso,ic95_alto,"
               "k_anno_finale,p_anno_finale,ic95_basso_finale,ic95_alto_finale,"
               "n_run_con_anni_valutabili,lambda_picco_media,lambda_picco_se,lambda_picco_sd,"
               "lambda_picco_min,lambda_picco_mediana,lambda_picco_max,lambda_medio_media,anni_valutabili_medi\n";
    }
    for (const auto& riga : righeCSV) out << riga << "\n";
    std::cout << "\nRiepilogo accodato a " << percorso << '\n';
}
