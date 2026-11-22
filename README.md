# Simulazione Stocastica della Sindrome di Kessler

Questo repository contiene l'implementazione in **C++** di un modello stocastico per la simulazione della **Sindrome di Kessler**. Il progetto modella l'ambiente dell'orbita terrestre bassa (LEO) utilizzando un approccio macroscopico al fine di prevedere l'innesco dell'instabilità collisionale spaziale.

##  Il Modello Fisico-Matematico

L'orbita è suddivisa in gusci orbitali concentrici. All'interno di ogni guscio, la popolazione di oggetti spaziali è suddivisa in quattro nodi:

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

###  Condizione di Criticità e Innesco
La stabilità del sistema spaziale in un dato guscio è determinata analiticamente dalle proprietà spettrali della matrice $A_k$. Il punto di innesco della Sindrome di Kessler si verifica quando il raggio spettrale $\lambda_{max}$ della matrice supera l'unità:

$$\lambda_{max}(A_k) > 1$$

##  La Simulazione C++

Il programma contenuto in questo repository esegue una simulazione Monte Carlo per campionare stocasticamente gli eventi di collisione all'interno di un guscio orbitale specifico. 

Il modello è stato testato utilizzando i dati reali della costellazione Starlink di SpaceX (quote 300-600 km).


##  Guida all'Uso

Per compilare ed eseguire il progetto, è necessario avere installato la libreria Eigen.
