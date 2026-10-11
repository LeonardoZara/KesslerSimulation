# KesslerSimulation
Kessler sindrome stochastic model


# Simulazione Stocastica della Sindrome di Kessler

Questo repository contiene l'implementazione in **C++** di un modello stocastico per la simulazione della **Sindrome di Kessler**. Il progetto modella l'ambiente dell'Orbita Terrestre Bassa (LEO) utilizzando un approccio macroscopico basato sull'algebra delle catene di Markov e sulla teoria delle reti, al fine di prevedere matematicamente l'innesco dell'instabilità collisionale spaziale.

Il modello e la teoria alla base di questa simulazione sono descritti in dettaglio nella tesina: *"Modello stocastico della Sindrome di Kessler"* (Leonardo Zara).

## 🛰️ Il Modello Fisico-Matematico

Lo spazio delle fasi (l'orbita LEO) è suddiviso in gusci orbitali concentrici. All'interno di ogni guscio, la popolazione di oggetti spaziali non è omogenea, ma suddivisa in quattro nodi (classi fisiche):

1. **$P$ (Payloads):** Satelliti attivi o manovrabili.
2. **$U$ (Upper Stages):** Stadi superiori dei vettori abbandonati.
3. **$N$ (Non-Maneuverable):** Satelliti inattivi o guasti.
4. **$F$ (Fragments):** Detriti spaziali (il pozzo entropico del sistema).

L'evoluzione temporale del sistema è governata dall'equazione di stato matriciale linearizzata lungo la traiettoria:

$$X_{k+1} \approx A_k X_k + g_k(t)$$

Dove:
* $X_k$ è il vettore colonna $[P_k, U_k, N_k, F_k]^T$ al tempo $t_k$.
* $A_k$ è la matrice di transizione di stato che codifica la topologia della rete, includendo il decadimento atmosferico (sulla diagonale) e la generazione non-lineare di frammenti dovuta alle collisioni (fuori diagonale).
* $g_k(t)$ è il vettore delle forzanti esterne (nuovi lanci).

### ⚠️ Condizione di Criticità e Innesco
La stabilità del sistema spaziale in un dato guscio è determinata analiticamente dalle proprietà spettrali della matrice $A_k$. Il punto di innesco della Sindrome di Kessler equivale a una transizione di fase topologica, che si verifica quando il raggio spettrale $\lambda_{max}$ della matrice supera l'unità:

$$\lambda_{max}(A_k) > 1$$

## 💻 La Simulazione C++

Il programma contenuto in questo repository esegue una **simulazione Monte Carlo** per campionare stocasticamente gli eventi di collisione all'interno di un guscio orbitale specifico. 

**Flusso principale del software:**
1. Inizializzazione delle popolazioni $P_0, U_0, N_0, F_0$ e dei parametri fisici (tassi di guasto $\gamma$, decadimento atmosferico $\epsilon$, moltiplicatori di massa $c$).
2. Stima delle frequenze di collisione empiriche.
3. Popolamento dinamico dei coefficienti della matrice non-stocastica $A_k$.
4. Calcolo degli autovalori della matrice per estrarre il raggio spettrale $\lambda_{max}$.

### Caso di Studio: Costellazione Starlink
Il modello è stato testato utilizzando i dati reali della costellazione Starlink di SpaceX (quote 300-600 km), per valutare come l'elevata densità di satelliti attivi influenzi i termini collisionali non lineari, contrapposta al rapido decadimento atmosferico a quelle quote.

---

## 🚀 Guida all'Uso

### Prerequisiti e Dipendenze
Per compilare ed eseguire il progetto, è necessario avere installato:
* Un compilatore C++ (es. `g++` o `clang++`) con supporto a C++[INSERISCI VERSIONE, es. 17 o 20].
* `CMake` (versione [INSERISCI VERSIONE] o superiore) - *[Lascia questa riga solo se usi CMake]*
* [DA COMPILARE: Inserisci eventuali librerie esterne utilizzate, ad esempio Eigen per il calcolo delle matrici e degli autovalori, e indica come installarle].

### Compilazione
```bash
# Clona il repository
git clone https://github.com/LeonardoZara/KesslerSimulation.git
cd KesslerSimulation

# [DA COMPILARE: Inserisci qui i comandi esatti per compilare il tuo progetto]
# Esempio se usi g++ direttamente:
# g++ -o kessler_sim main.cpp -O3

# Esempio se usi CMake:
# mkdir build && cd build
# cmake ..
# make
```

### Esecuzione
```bash
# [DA COMPILARE: Inserisci il comando per eseguire il programma]
# Esempio:
# ./kessler_sim --config starlink_params.txt
```

---

## 📊 Output e Risultati Attesi

Eseguendo la simulazione con i parametri di default (scenario sub-critico e super-critico), il programma restituirà in output:
* Le matrici di transizione $A_{subcritica}$ e $A_{supercritica}$.
* Gli autovalori calcolati per ogni step della simulazione.
* Il tracciamento del superamento della soglia $\lambda_{max} = 1$.

*[DA COMPILARE: Se il tuo programma genera file .csv, grafici o log particolari, specificalo qui. Ad esempio: "Il programma genera un file `results.csv` contenente l'andamento di lambda rispetto all'aumento dei frammenti F, utile per generare plot."]*

## 📌 Limitazioni del Modello
Come discusso nella tesina, la simulazione attuale presenta alcune limitazioni di cui tenere conto:
* Assenza di un database esatto per i frammenti $F$ di piccole dimensioni (trattati in modo parametrico).
* Fluttuazioni del decadimento atmosferico dovute al ciclo solare undecennale non attualmente integrate dinamicamente.
* Mancato scambio di massa inter-shell (il modello analizza un guscio per volta).

## Autore
**Leonardo Zara** - *Collegio Superiore* - Micro macro computational models.
