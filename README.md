# MacroMisuraG

# Analisi Dati: Esperimento di Cavendish (Misura di $G$)

Questo repository contiene la macro ROOT in C++ utilizzata per l'analisi dei dati dell'esperimento di Cavendish per la determinazione della costante di gravitazione universale $G$.

## Obiettivo dell'Analisi

Il codice effettua un fit non lineare sui dati della posizione dello spot luminoso (in cm) in funzione del tempo (in secondi). Il modello implementato tiene conto di un'oscillazione sinusoidale smorzata sovrapposta a una deriva lineare dell'equilibrio:
$$f(x) = p_0 \cdot \sin(p_1 \cdot x + p_2) \cdot e^{-\frac{x}{p_3}} + p_4 + p_5 \cdot x$$

### Significato dei Parametri:
* `p0`: Ampiezza iniziale dell'oscillazione.
* `p1`: Pulsazione ($\omega = \frac{2\pi}{T}$), utile per ricavare il periodo d'oscillazione $T$.
* `p2`: Angolo di fase iniziale.
* `p3`: Tempo di smorzamento caratteristico ($\tau$).
* `p4`: Posizione di equilibrio iniziale del manubrio a $t = 0$.
* `p5`: Coefficiente angolare della retta di deriva.

---

## Struttura dei Dati Richiesti

Per far girare la macro, è necessario avere nella stessa cartella un file di testo chiamato **`dati.txt`**. Il file deve essere strutturato in 4 colonne separate da spazi o tabulazioni:

```text
# Tempo(s)  Posizione(cm)
0           25.3        
10          25.8          
20          26.4
