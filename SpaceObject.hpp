#pragma once

// Definiamo i tipi di oggetti spaziali. Usare 'enum class' è molto più sicuro
// in C++ moderno rispetto ai vecchi enum, perché impedisce conversioni accidentali in numeri.
enum class ObjectType {
    ACTIVE_PAYLOAD=0, // Satelliti operativi con propulsione (s_CAM)
    DEAD_PAYLOAD=1,   // Satelliti guasti o a fine vita (passivi)
    ROCKET_BODY=2,    // Secondi stadi dei razzi (massicci, passivi)
    FRAGMENT=3        // Detriti da frammentazione (piccoli, passivi)
};

// La struct base. Deve essere il più leggera possibile in memoria.
struct SpaceObject {
    ObjectType type;
    
    // Parametri fisici
    double crossSection;    // Area in m^2 (sigma: serve per calcolare la probabilità di urto)
    double mass;            // Massa in kg (serve per calcolare il numero di frammenti post-urto)
    
    // Parametro topologico (fondamentale per estrarre la Matrice Empirica)
    int initialShellIndex;  // Salva l'indice del guscio in cui l'oggetto si trovava
                            // all'inizio dell'anno simulato.

    // --- FUNZIONI LOGICHE (Helper) ---
    
    // Verifica se l'oggetto può tentare una manovra di elusione (s_CAM)
    bool isManeuverable() const {
        return type == ObjectType::ACTIVE_PAYLOAD;
    }

    // Verifica se l'oggetto è abbastanza grande da essere tracciato dai radar a terra.
    // Se non è tracciabile, un Active Payload non può schivarlo!
    // (Soglia tipica: 10 cm di diametro, approssimabile a ~0.01 m^2 di sezione d'urto)
    bool isTrackable() const {
        return crossSection > 0.01; 
    }
};