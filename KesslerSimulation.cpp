#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include "KesslerSimulation.hpp"


// Indici dei tipi: 0 = P (attivo), 1 = N (morto), 2 = U (razzo), 3 = F (frammento)
// Indici dei gusci: 0 = atmosfera, 1 = 300-400, 2 = 400-500, 3 = 500-600 km
namespace {
    const double secondsInYear = 3.15576e7;
    const double v_rel = 10.0;    // km/s, velocita' relativa media in LEO
 
    const double S_CAM = 0.999;

    const double GAMMA = 0.01; // Guasti dei satelliti attivi: gamma = guasto a fine missione
    const double KAPPA = 5.3; // kappa = perdita di manovrabilita' per piccoli urti con i frammenti.
 
    const double SOGLIA_J_PER_G = 40.0; // Soglia di frammentazione catastrofica (PDF, dopo Eq. 5): 40 J per grammo di bersaglio.
 
    const double C_FRAMMENTI[4] = {500.0, 500.0, 3000.0, 1.0}; // Frammenti generati dalla distruzione di un oggetto di ciascun tipo (PDF §5/§6.3):
 
    // Tasso di decadimento verso il guscio sotto, epsilon[tipo][guscio], per anno.
    //  - P: 0, i satelliti attivi mantengono la quota con i propulsori.
    //  - F: PDF §6.3 (0.90 / 0.75 / 0.60).
    //  - N e U: Eq. 28 con A/m = 0.075 e rho(h) = 3.7e-12 * exp(-(h-400)/60) (3.7e-12 é la media dei due valori del PDF a 400 km; la scala di 60 km é una mia stima).
    //    Risultato 10.4 / 1.98 / 0.38; i primi due superano 1 (escono dal guscio in meno di un anno) e li limito a 1.
    const double EPSILON[4][4] = {
        {0.0, 0.0,  0.0,  0.0 },   // P
        {0.0, 1.0,  1.0,  0.38},   // N
        {0.0, 1.0,  1.0,  0.38},   // U
        {0.0, 0.90, 0.75, 0.60}    // F
    };
 
    // Lanci di nuovi satelliti attivi per anno e per guscio (forzante esterna g_k del PDF).
    // Ho scelto gamma * P0 (9, 85, 4): rimpiazzano i satelliti che si guastano, cosi' la
    // flotta resta circa costante. Metti 0 per non lanciare nulla.
    const double LANCI_P[4] = {0.0, 9.0, 85.0, 4.0};
    const double payloadSigma = 130.0;  // m^2, come in main
    const double payloadMass   = 800.0;  // kg,  come in main
 
    // Tetto di sicurezza al numero di coppie candidate per guscio e per anno
    // (poisson_distribution<int> va in overflow oltre ~2e9).
    const double maxCandidates = 1.0e8;
 
    // Sezione d'urto geometrica di una coppia: pi*(r1+r2)^2 con r = sqrt(A/pi)
    // si semplifica in (sqrt(A1) + sqrt(A2))^2. Aree in m^2, risultato in m^2.
    double sigmaCoppia(double areaA, double areaB) {
        double somma = std::sqrt(areaA) + std::sqrt(areaB);
        return somma * somma;
    }
}





void KesslerSimulation::addObjectsToShell(int shellIndex, int count, ObjectType type, double crossSection, double mass) {
    if (shellIndex < 0 || shellIndex >= numShells) return;
    
    for (int i = 0; i < count; ++i) {
        SpaceObject obj;
        obj.type = type;
        obj.crossSection = crossSection;
        obj.mass = mass;
        obj.initialShellIndex = shellIndex; // Verrà comunque sovrascritto dalla Fase 1
        obj.initialType = type;
        shells[shellIndex].objects.push_back(obj);
    }
}
 
void KesslerSimulation::step(double dt) {
    // FASE 1: Snapshot Iniziale
    // Aggiorniamo l'indice iniziale di tutti gli oggetti e diciamo al tracker di prepararsi
    for (int i = 0; i < numShells; ++i) {
        for (auto& obj : shells[i].objects) {
            obj.initialShellIndex = i;
            obj.initialType = obj.type;
        }
    }
    tracker.reset(shells);
 
    processCollisions(dt);
    processFailures(dt);
    processDecay(dt);
    processLaunches(dt);

    // Costruzione delle Matrici dei Flussi (una per guscio)
    // Contiamo i SOPRAVVISSUTI di ogni guscio: gli oggetti che a inizio anno erano
    // in questo stesso guscio e ci sono ancora (diagonale della matrice: a_ii).
    // Restano fuori:
    //  - gli oggetti arrivati da un altro guscio per decadimento (sono un input esterno, g_k);
    //  - i frammenti nati quest'anno (initialShellIndex = -1): sono gia' stati
    //    registrati in generateFragments, e non vanno contati due volte.
    for (int s = 0; s < numShells; ++s) {
        for (const auto& obj : shells[s].objects) {
            if (obj.initialShellIndex == s) {
                tracker.recordFlow(s, obj.initialType, obj.type, 1.0);
            }
        }
    }
 
    shells[0].clear();
    currentTime += dt;
}

void KesslerSimulation::processDecay(double dt) {
    // Distribuzione uniforme per il lancio del dado [0.0, 1.0)
    std::uniform_real_distribution<double> dist(0.0, 1.0);
 
    // Iteriamo dal guscio più basso (1) a quello più alto.
    // Il guscio 0 è l'atmosfera (stato assorbente), non calcoliamo il decadimento da lì.
    for (int i = 1; i < numShells; ++i) {
        std::vector<SpaceObject> survivors;
        survivors.reserve(shells[i].objects.size());
 
        for (const auto& obj : shells[i].objects) {

            // Gli oggetti nati quest'anno (frammenti, lanci: initialShellIndex = -1)
            // non decadono nell'anno di nascita: nel PDF a_41..a_44 contano i frammenti
            // generati, senza sottrarre il loro decadimento, che avviene dall'anno dopo.
            if (obj.initialShellIndex != i) {
                survivors.push_back(obj);
                continue;
            }

            // Probabilita' di decadimento dalla tabella EPSILON[tipo][guscio].
            // COSA: prima 0.05 per i frammenti e 0.001 per tutto il resto, uguale in ogni guscio.
            // PERCHé: il PDF usa un epsilon diverso per guscio e per tipo (a bassa quota
            // l'attrito é molto piu' forte); epsilon é la frazione che scende ogni anno,
            // come in a_ii = 1 - epsilon - ..., quindi p = epsilon*dt (al massimo 1).
            // Si usa il tipo INIZIALE: un satellite che si guasta quest'anno decade come P
            // (epsilon_P) e come N (epsilon_N) solo dall'anno dopo, come impone la matrice.
            double p_decay = std::min(1.0, EPSILON[static_cast<int>(obj.initialType)][i] * dt);

            if (dist(rng) < p_decay) {
                // Decade! Va nel guscio sottostante (i - 1)
                shells[i - 1].objects.push_back(obj);
            } else {
                // Sopravvive! Resta nel guscio attuale
                survivors.push_back(obj);
            }
        }
        // Aggiornamento del guscio con solo con i sopravvissuti
        shells[i].objects = std::move(survivors);
    }
}


void KesslerSimulation::processFailures(double dt) {
    // Guasti dei satelliti attivi, termini gamma e kappa di a_11 e a_31.
    // Ogni P diventa N (resta dov'è, cambia solo tipo) con probabilità p = gamma*dt + kappa * (urti attesi P-F per satellite).
    // Gli urti attesi P-F per satellite sono (F/V)*sigma(P,F)*v*dt; kappa li moltiplica perché i piccoli urti non distruttivi sono molti di piu' di quelli esplicitamente simulati.
    std::uniform_real_distribution<double> dist(0.0, 1.0);
 
    for (int i = 1; i < numShells; ++i) {
        // Frammenti presenti a inizio anno (X_k): quelli nati negli urti di quest'anno hanno initialShellIndex = -1 e non contano.
        double nF = 0.0;
        for (const auto& obj : shells[i].objects) {
            if (obj.type == ObjectType::FRAGMENT && obj.initialShellIndex == i) nF += 1.0;
        }
        double V = shells[i].volume;
 
        for (auto& obj : shells[i].objects) {
            if (obj.type != ObjectType::ACTIVE_PAYLOAD) continue;
 
            double sigma_km2 = sigmaCoppia(obj.crossSection, fragmentSigma) * 1.0e-6;
            double urti_PF = (nF / V) * sigma_km2 * v_rel * dt * secondsInYear;
            double p_guasto = std::min(1.0, GAMMA * dt + KAPPA * urti_PF);
 
            if (dist(rng) < p_guasto) {
                obj.type = ObjectType::DEAD_PAYLOAD;
            }
        }
    }
}
 
void KesslerSimulation::processLaunches(double dt) {
    // Aggiungiamo LANCI_P[guscio]*dt satelliti con initialShellIndex = -1.
    SpaceObject nuovo;
    nuovo.type = ObjectType::ACTIVE_PAYLOAD;
    nuovo.crossSection = payloadSigma;
    nuovo.mass = payloadMass;
    nuovo.initialShellIndex = -1;
    nuovo.initialType = ObjectType::ACTIVE_PAYLOAD;
 
    for (int i = 1; i < numShells; ++i) {
        int quanti = static_cast<int>(std::lround(LANCI_P[i] * dt));
        if (quanti > 0) {
            shells[i].objects.insert(shells[i].objects.end(), quanti, nuovo);
        }
    }
}

void KesslerSimulation::processCollisions(double dt) {
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
 
    for (int i = 1; i < numShells; ++i) { // Ignoriamo l'atmosfera (guscio 0)
        int N = static_cast<int>(shells[i].objects.size());
        if (N < 2) continue; // Servono almeno 2 oggetti per scontrarsi!
 
        double V = shells[i].volume;
 
        // Calcolo del Valore Atteso (Media teorica degli urti)
        // Vengono prese le coppie candidate con il sigma MASSIMO del guscio (due oggetti piu' grandi), poi ogni coppia viene accettata con probabilita'
        // sigma_coppia/sigma_max. Il risultato é esattamente un processo di Poisson con tasso proporzionale a sigma(i,j), senza dover sommare su tutte le coppie.
        double A_max = 0.0;
        for (const auto& obj : shells[i].objects) {
            if (obj.crossSection > A_max) A_max = obj.crossSection;
        }
        double sigma_max_m2 = sigmaCoppia(A_max, A_max);
        double sigma_max_km2 = sigma_max_m2 * 1.0e-6;   // m^2 -> km^2
 
        // Formula della teoria cinetica dei gas per collisioni intra-guscio
        double candidati_attesi = 0.5 * (static_cast<double>(N) * (N - 1) / V) * sigma_max_km2 * v_rel * dt * secondsInYear;
 
        // Tetto di sicurezza: poisson_distribution<int> va in overflow oltre ~2e9.
        candidati_attesi = std::min(candidati_attesi, maxCandidates);
        if (candidati_attesi <= 0.0) continue;
 
        // Estrazione del numero vero di coppie candidate da una Poissoniana
        std::poisson_distribution<int> poisson_dist(candidati_attesi);
        int candidati = poisson_dist(rng);
 
        if (candidati == 0) continue;
 
        std::vector<char> distrutto(N, 0);
        std::uniform_int_distribution<int> object_picker(0, N - 1);
 
        for (int c = 0; c < candidati; ++c) {
            // 2 oggetti a caso (solo tra i N presenti a inizio anno: i frammenti appena nati sono in coda e non vengono mai pescati)
            int idx1 = object_picker(rng);
            int idx2 = object_picker(rng);
            while (idx1 == idx2) idx2 = object_picker(rng); // Assicurati che siano diversi
 
            if (distrutto[idx1] || distrutto[idx2]) continue;
 
            const SpaceObject obj1 = shells[i].objects[idx1];
            const SpaceObject obj2 = shells[i].objects[idx2];
 
            // Viene accettata la coppia con probabilita' sigma_coppia / sigma_max
            double sigma_coppia_m2 = sigmaCoppia(obj1.crossSection, obj2.crossSection);
            if (uniform_dist(rng) * sigma_max_m2 >= sigma_coppia_m2) continue;
 
            bool collisionOccurs = true;
 
            if (obj1.isManeuverable() && obj2.isTrackable()) {
                if (uniform_dist(rng) < S_CAM) collisionOccurs = false;
            } else if (obj2.isManeuverable() && obj1.isTrackable()) {
                if (uniform_dist(rng) < S_CAM) collisionOccurs = false;
            }
            if (!collisionOccurs) continue;

            // L'urto è catastrofico solo se l'energia del proiettile (il piu' leggero) per grammo di bersaglio (il piu' pesante) supera 40 J/g.
            // E/m_bersaglio = 0.5 * (m_leggero/m_pesante) * v^2, con v in m/s; dividendo per 1000 si passa da J/kg a J/g.
            double m_leggero = std::min(obj1.mass, obj2.mass);
            double m_pesante = std::max(obj1.mass, obj2.mass);
            double v_m_s = v_rel * 1000.0;
            double energia_J_per_g = 0.5 * (m_leggero / m_pesante) * v_m_s * v_m_s / 1000.0;
            if (energia_J_per_g < SOGLIA_J_PER_G) continue;
            distrutto[idx1] = 1;
            distrutto[idx2] = 1;
            generateFragments(i, obj1, obj2);
        }
 
        // Rimozione oggetti distrutti dal guscio
        std::vector<SpaceObject> rimasti;
        rimasti.reserve(shells[i].objects.size());
        for (int k = 0; k < static_cast<int>(shells[i].objects.size()); ++k) {
            if (k >= N || !distrutto[k]) rimasti.push_back(shells[i].objects[k]);
        }
        shells[i].objects = std::move(rimasti);
    }
}
 
void KesslerSimulation::generateFragments(int shellIndex, const SpaceObject& obj1, const SpaceObject& obj2) {
    // Numero di frammenti = c_tipo1 + c_tipo2
    int numNewFragments = static_cast<int>(C_FRAMMENTI[static_cast<int>(obj1.type)] + C_FRAMMENTI[static_cast<int>(obj2.type)]);
 
    SpaceObject frag;
    frag.type = ObjectType::FRAGMENT;
    frag.crossSection = fragmentSigma;
    frag.mass = fragmentMass;
    // -1 = "nato quest'anno": step() non lo contera' tra i sopravvissuti
    frag.initialShellIndex = -1;
    frag.initialType = ObjectType::FRAGMENT;
 
    // I frammenti nascono subito in coda al guscio (non vengono pescati in questo anno, vedi object_picker); inseriti tutti insieme con una sola chiamata.
    shells[shellIndex].objects.insert(shells[shellIndex].objects.end(), numNewFragments, frag);
 
    // Registrazione nel tracker: meta' dei frammenti a ciascun genitore.
    // il termine collisionale é bilineare e va diviso tra le due colonne.
    // Riga = FRAGMENT (destinazione), colonna = tipo del genitore (origine).
    tracker.recordFlow(shellIndex, obj1.type, ObjectType::FRAGMENT, 0.5 * numNewFragments);
    tracker.recordFlow(shellIndex, obj2.type, ObjectType::FRAGMENT, 0.5 * numNewFragments);
}