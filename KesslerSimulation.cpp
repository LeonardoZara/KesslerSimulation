#include <iostream>
#include <cmath>

#include "KesslerSimulation.hpp"

 
void KesslerSimulation::addObjectsToShell(int shellIndex, int count, ObjectType type, double crossSection, double mass) {
    if (shellIndex < 0 || shellIndex >= numShells) return;
    
    for (int i = 0; i < count; ++i) {
        SpaceObject obj;
        obj.type = type;
        obj.crossSection = crossSection;
        obj.mass = mass;
        obj.initialShellIndex = shellIndex; // Verrà comunque sovrascritto dalla Fase 1
        
        shells[shellIndex].objects.push_back(obj);
    }
}
 
void KesslerSimulation::step(double dt) {
    // FASE 1: Snapshot Iniziale
    // Aggiorniamo l'indice iniziale di tutti gli oggetti e diciamo al tracker di prepararsi
    for (int i = 0; i < numShells; ++i) {
        for (auto& obj : shells[i].objects) {
            obj.initialShellIndex = i;
        }
    }
    tracker.reset(shells);
 
    // FASE 2: Fisica del Decadimento
    processDecay(dt);
 
    // FASE 3: Fisica delle Collisioni
    processCollisions(dt);
 
    // FASE 4: Costruzione delle Matrici dei Flussi (una per guscio)
    // Contiamo i SOPRAVVISSUTI di ogni guscio: gli oggetti che a inizio anno erano
    // in questo stesso guscio e ci sono ancora (diagonale della matrice: a_ii).
    // Restano fuori:
    //  - gli oggetti arrivati da un altro guscio per decadimento (sono un input esterno, g_k);
    //  - i frammenti nati quest'anno (initialShellIndex = -1): sono gia' stati
    //    registrati in generateFragments, e non vanno contati due volte.
    for (int s = 0; s < numShells; ++s) {
        for (const auto& obj : shells[s].objects) {
            if (obj.initialShellIndex == s) {
                tracker.recordFlow(s, obj.type, obj.type, 1.0);
            }
        }
    }
 
    currentTime += dt;
}
 
// --- IMPLEMENTAZIONE DELLA FISICA PRIVATA ---
 
void KesslerSimulation::processDecay(double dt) {
    // Distribuzione uniforme per il lancio del dado [0.0, 1.0)
    std::uniform_real_distribution<double> dist(0.0, 1.0);
 
    // Iteriamo dal guscio più basso (1) a quello più alto.
    // Il guscio 0 è l'atmosfera (stato assorbente), non calcoliamo il decadimento da lì.
    for (int i = 1; i < numShells; ++i) {
        std::vector<SpaceObject> survivors;
        survivors.reserve(shells[i].objects.size());
 
        for (const auto& obj : shells[i].objects) {
            // Calcolo Probabilità di decadimento (Semplificata per l'esempio)
            // Nella tesi citi epsilon_F = 0.05. I satelliti attivi/morti decadono molto meno.
            double p_decay = 0.0;
            if (obj.type == ObjectType::FRAGMENT) {
                p_decay = 0.05 * dt; 
            } else {
                p_decay = 0.001 * dt; // Satelliti grandi e pesanti risentono meno dell'attrito a quote alte
            }
 
            // Lancio del dado!
            if (dist(rng) < p_decay) {
                // Decade! Lo spingiamo nel guscio sottostante (i - 1)
                shells[i - 1].objects.push_back(obj);
            } else {
                // Sopravvive! Resta nel guscio attuale
                survivors.push_back(obj);
            }
        }
        // Aggiorniamo il guscio solo con i sopravvissuti
        shells[i].objects = std::move(survivors);
    }
}
 
void KesslerSimulation::processCollisions(double dt) {
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double v_rel = 10.0; // km/s (Come dimostrato nella tua tesi)
 
    for (int i = 1; i < numShells; ++i) { // Ignoriamo l'atmosfera (guscio 0)
        double N = shells[i].objects.size();
        if (N < 2) continue; // Servono almeno 2 oggetti per scontrarsi!
 
        double V = shells[i].volume;
        double density = N / V;
 
        // 1. Calcolo del Valore Atteso (Media teorica degli urti)
        // Sigma_media semplificata (in un modello reale faresti la media esatta)
        double sigma_avg = 0.01; // km^2 (10 m^2)
        
        // Formula della teoria cinetica dei gas per collisioni intra-guscio
        double expected_collisions = 0.5 * (N * (N - 1) / V) * sigma_avg * v_rel * dt;
 
        // 2. Stocasticità: Estraiamo il numero VERO di urti da una Poissoniana
        std::poisson_distribution<int> poisson_dist(expected_collisions);
        int actual_collisions = poisson_dist(rng);
 
        if (actual_collisions == 0) continue; // Che fortuna, nessun urto quest'anno!
 
        std::vector<int> to_destroy; // Indici degli oggetti distrutti
        std::uniform_int_distribution<int> object_picker(0, N - 1);
 
        for (int c = 0; c < actual_collisions; ++c) {
            // Peschiamo 2 oggetti a caso
            int idx1 = object_picker(rng);
            int idx2 = object_picker(rng);
            while (idx1 == idx2) idx2 = object_picker(rng); // Assicurati che siano diversi
 
            const SpaceObject obj1 = shells[i].objects[idx1];
            const SpaceObject obj2 = shells[i].objects[idx2];
 
            // 3. Meccanica di Evitamento (s_CAM)
            bool collisionOccurs = true;
            double s_CAM = 0.999; // Probabilità di schivata riuscita
 
            if (obj1.isManeuverable() && obj2.isTrackable()) {
                if (uniform_dist(rng) < s_CAM) collisionOccurs = false;
            } else if (obj2.isManeuverable() && obj1.isTrackable()) {
                if (uniform_dist(rng) < s_CAM) collisionOccurs = false;
            }
 
            if (collisionOccurs) {
                // Li segniamo per la distruzione
                to_destroy.push_back(idx1);
                to_destroy.push_back(idx2);
                
                // 4. Generazione dei frammenti (Cascata di Kessler)
                generateFragments(i, obj1, obj2);
            }
        }
 
        // 5. Pulizia: Rimuoviamo gli oggetti distrutti dal guscio
        // (Usiamo un trucco per non sballare gli indici mentre eliminiamo)
        if (!to_destroy.empty()) {
            std::sort(to_destroy.begin(), to_destroy.end(), std::greater<int>());
            to_destroy.erase(std::unique(to_destroy.begin(), to_destroy.end()), to_destroy.end());
            
            for (int idx : to_destroy) {
                shells[i].objects.erase(shells[i].objects.begin() + idx);
            }
        }
    }
}
 
void KesslerSimulation::generateFragments(int shellIndex, const SpaceObject& obj1, const SpaceObject& obj2) {
    // Formula base empirica: Num frammenti proporzionale alla massa totale
    // C_p e C_u citati nella tua tesi (500 per satelliti, 3000 per upper stages)
    double totalMass = obj1.mass + obj2.mass;
    int numNewFragments = static_cast<int>(totalMass * 0.5); // Es. 1 kg = 0.5 frammenti macroscopici
 
    for (int i = 0; i < numNewFragments; ++i) {
        SpaceObject frag;
        frag.type = ObjectType::FRAGMENT;
        frag.crossSection = 0.005; // 50 cm^2
        frag.mass = 1.0;           // 1 kg
        
        // FONDAMENTALE PER LA MATRICE: I frammenti "ereditano" il TIPO dai genitori.
        // Scegliamo a caso di chi è "figlio" questo frammento: la sua colonna (origine)
        // nella matrice del guscio sara' il tipo di quel genitore (a_41, a_42, a_43, a_44).
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        const SpaceObject& genitore = (dist(rng) < 0.5) ? obj1 : obj2;
 
        // -1 = "nato quest'anno": la FASE 4 di step() non lo contera' tra i sopravvissuti
        frag.initialShellIndex = -1;
 
        // Nota: non li aggiungiamo subito a 'shells[shellIndex].objects' per non farli 
        // scontrare nel turno in cui sono appena nati (e non invalidare l'iteratore).
        // Nel codice di produzione reale li aggiungeresti a un vettore temporaneo 
        // e li uniresti alla fine di processCollisions. Per semplicità qui li mettiamo subito in coda.
        shells[shellIndex].objects.push_back(frag);
        
        // E lo registriamo nel tracker, nella matrice del guscio dove e' avvenuto l'urto:
        // riga = FRAGMENT (destinazione), colonna = tipo del genitore (origine)
        tracker.recordFlow(shellIndex, genitore.type, ObjectType::FRAGMENT, 1.0);
    }
}