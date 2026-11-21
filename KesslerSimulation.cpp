#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include "KesslerSimulation.hpp"

// Indici dei tipi: 0 = P (attivo), 1 = N (morto), 2 = U (razzo), 3 = F (frammento)
// Indici dei gusci: 0 = atmosfera, 1 = 300-400, 2 = 400-500, 3 = 500-600 km
namespace {
    const double secondsInYear = 3.15576e7;
    const double v_rel = 9.73;    // km/s, velocita' relativa media in LEO
    const double S_CAM = 0.999;
    const double GAMMA = 0.01; // Guasti dei satelliti attivi: gamma = guasto a fine missione
    const double KAPPA = 5.3; // kappa = perdita di manovrabilita' per piccoli urti con i frammenti.
 
    const double SOGLIA_J_PER_G = 40.0; // Soglia di frammentazione catastrofica: 40 J per grammo di bersaglio.
 
    const double C_FRAMMENTI[4] = {500.0, 500.0, 3000.0, 1.0}; // Frammenti generati dalla distruzione di un oggetto di ciascun tipo
 
    const double EPSILON[4][4] = {
        {0.0, 0.0,  0.0,  0.0 },   // P
        {0.0, 1.0,  1.0,  0.27},   // N
        {0.0, 0.61, 0.1,  0.02},   // U
        {0.0, 0.85, 0.14, 0.03}    // F
    };

    // Tetto di sicurezza al numero di coppie candidate per guscio e per anno.
    const double maxCandidates = 1.0e8;
 
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
        obj.initialShellIndex = shellIndex;
        obj.initialType = type;
        shells[shellIndex].objects.push_back(obj);
    }
}
 
void KesslerSimulation::setLaunches(const std::vector<double>& ratesPerYear, double crossSection, double mass) {
    if (static_cast<int>(ratesPerYear.size()) != numShells) {
        throw std::invalid_argument("setLaunches: servono esattamente numShells valori (uno per guscio, atmosfera inclusa)");
    }
    launchRates = ratesPerYear;
    launchRates[0] = 0.0; // nessun lancio "nell'atmosfera"
    launchCrossSection = crossSection;
    launchMass = mass;
}

void KesslerSimulation::saveInitialState() {
    initialShells = shells;
}

void KesslerSimulation::restoreInitialState(std::uint64_t baseSeed, std::uint64_t runIndex) {
    if (initialShells.size() != shells.size()) {
        throw std::logic_error("restoreInitialState: chiamare prima saveInitialState()");
    }
    shells = initialShells;
    tracker = TransitionTracker(numShells);
    currentTime = 0.0;
    // RNG: seed deterministico (baseSeed, runIndex).
    // Run diverse danno flussi casuali diversi; stesse (baseSeed, runIndex) allora stessa traiettoria.
    std::seed_seq seq{ static_cast<std::uint32_t>(baseSeed & 0xFFFFFFFFu),
                       static_cast<std::uint32_t>(baseSeed >> 32),
                       static_cast<std::uint32_t>(runIndex & 0xFFFFFFFFu),
                       static_cast<std::uint32_t>(runIndex >> 32) };
    rng.seed(seq);
}

void KesslerSimulation::step(double dt) {
    //Snapshot iniziale
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

    // Costruzione delle matrici dei flussi (una per guscio)
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
    // Distribuzione uniforme per il lancio del random[0.0, 1.0)
    std::uniform_real_distribution<double> dist(0.0, 1.0);
 
    // Il guscio 0 è l'atmosfera, non serve calcolare il decadimento da lì.
    for (int i = 1; i < numShells; ++i) {
        std::vector<SpaceObject> survivors;
        survivors.reserve(shells[i].objects.size());
 
        for (const auto& obj : shells[i].objects) {

            // Gli oggetti nati quest'anno (frammenti, lanci: initialShellIndex = -1) non decadono nell'anno di nascita
            if (obj.initialShellIndex != i) {
                survivors.push_back(obj);
                continue;
            }

            // Probabilita' di decadimento dalla tabella EPSILON[tipo][guscio].
            double p_decay = std::min(1.0, EPSILON[static_cast<int>(obj.initialType)][i] * dt);

            if (dist(rng) < p_decay) {
                shells[i - 1].objects.push_back(obj);
            } else {
                // Resta nel guscio attuale
                survivors.push_back(obj);
            }
        }
        // Aggiornamento del guscio con solo con i sopravvissuti
        shells[i].objects = std::move(survivors);
    }
}


void KesslerSimulation::processFailures(double dt) {
    // Guasti dei satelliti attivi, termini gamma e kappa di a_11 e a_31.
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
    SpaceObject nuovo;
    nuovo.type = ObjectType::ACTIVE_PAYLOAD;
    nuovo.crossSection = launchCrossSection;
    nuovo.mass = launchMass;
    nuovo.initialShellIndex = -1;
    nuovo.initialType = ObjectType::ACTIVE_PAYLOAD;

    for (int i = 1; i < numShells; ++i) {
        int quanti = static_cast<int>(std::lround(launchRates[i] * dt));
        if (quanti > 0) {
            shells[i].objects.insert(shells[i].objects.end(), quanti, nuovo);
        }
    }
}

void KesslerSimulation::processCollisions(double dt) {
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
 
    for (int i = 1; i < numShells; ++i) { 
        int N = static_cast<int>(shells[i].objects.size());
        if (N < 2) continue;
 
        double V = shells[i].volume;
 
        // Calcolo del Valore Atteso
        // Vengono prese le coppie candidate con il sigma massimo del guscio (due oggetti piu' grandi), poi ogni coppia viene accettata con probabilità
        // sigma_coppia/sigma_max. Il risultato é esattamente un processo di Poisson con tasso proporzionale a sigma(i,j), senza dover sommare su tutte le coppie.
        // Le coppie frammento-frammento (F-F) non vengono campionate, perchè un urto tra frammenti non aggiunge frammenti.
        std::vector<int> grandi; // indici degli oggetti che non sono frammenti (P, N, U)
        std::vector<int> frammenti;
        double A_max = 0.0;
        for (int k = 0; k < N; ++k) {
            const SpaceObject& obj = shells[i].objects[k];
            if (obj.type == ObjectType::FRAGMENT) frammenti.push_back(k); else grandi.push_back(k);
            if (obj.crossSection > A_max) A_max = obj.crossSection;
        }
        const double nG = static_cast<double>(grandi.size());
        const double nF = static_cast<double>(frammenti.size());
        const double coppieGG = 0.5 * nG * (nG - 1.0);   
        const double coppieGF = nG * nF;                
        const double coppie   = coppieGG + coppieGF;
        if (coppie <= 0.0) continue;

        double sigma_max_m2 = sigmaCoppia(A_max, A_max);
        double sigma_max_km2 = sigma_max_m2 * 1.0e-6;

        double candidati_attesi = (coppie / V) * sigma_max_km2 * v_rel * dt * secondsInYear;

        // Tetto di sicurezza: poisson_distribution<int> va in overflow oltre 2e9.
        // Se il tetto scatta gli urti sono SOTTOSTIMATI e i risultati non sono affidabili: avviso a terminale.
        if (candidati_attesi > maxCandidates) {
            static bool avvisato = false;
            if (!avvisato) {
                std::cerr << "ATTENZIONE: tetto maxCandidates raggiunto nel guscio " << i
                          << ": gli urti sono sottostimati, risultati non affidabili." << '\n';
                avvisato = true;
            }
            candidati_attesi = maxCandidates;
        }
        if (candidati_attesi <= 0.0) continue;
 
        // Estrazione del numero vero di coppie candidate da una poissoniana
        std::poisson_distribution<int> poisson_dist(candidati_attesi);
        int candidati = poisson_dist(rng);
 
        if (candidati == 0) continue;
 
        std::vector<char> distrutto(N, 0);
        std::uniform_int_distribution<int> pickG(0, std::max(0, static_cast<int>(grandi.size()) - 1));
        std::uniform_int_distribution<int> pickF(0, std::max(0, static_cast<int>(frammenti.size()) - 1));

        for (int c = 0; c < candidati; ++c) {
            int idx1, idx2;
            if (coppieGF <= 0.0 || uniform_dist(rng) * coppie < coppieGG) {
                idx1 = grandi[pickG(rng)];                   
                idx2 = grandi[pickG(rng)];
                while (idx1 == idx2) idx2 = grandi[pickG(rng)];
            } else {
                idx1 = grandi[pickG(rng)];                     
                idx2 = frammenti[pickF(rng)];
            }

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
            // E/m_bersaglio = 0.5 * (m_leggero/m_pesante) * v^2, con v in m/s;
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
    int numNewFragments = static_cast<int>(C_FRAMMENTI[static_cast<int>(obj1.type)] + C_FRAMMENTI[static_cast<int>(obj2.type)]);
 
    SpaceObject frag;
    frag.type = ObjectType::FRAGMENT;
    frag.crossSection = fragmentSigma;
    frag.mass = fragmentMass;
    frag.initialShellIndex = -1; // -1 = "nato quest'anno": step() non lo contera' tra i sopravvissuti
    frag.initialType = ObjectType::FRAGMENT;
 
    shells[shellIndex].objects.insert(shells[shellIndex].objects.end(), numNewFragments, frag);
 
    // Registrazione nel tracker: meta' dei frammenti a ciascun genitore.
    // Riga = FRAGMENT (destinazione), colonna = tipo del genitore (origine)
    tracker.recordFlow(shellIndex, obj1.type, ObjectType::FRAGMENT, 0.5 * numNewFragments);
    tracker.recordFlow(shellIndex, obj2.type, ObjectType::FRAGMENT, 0.5 * numNewFragments);
}