# Multigrid Verfahren für Differentialgleichungen 
### Aufruf
./mgsolve lvl iter
wobei:
lvl = Anzahl der Level, um die das Gitter verfeinert werden soll (Level 1 wäre das gröbste Gitter)
iter = Anzahl der Iterationen 

### Ausgabedatei
Es entsteht eine Ausgabe Datei solution.txt, in der die jeweiligen Lösungswerte für jede Stelle (x,y) des Gitters stehen. Das Gitter kann mit der splot Funktion von gnuplot visualisiert werden.

### Ausgabe auf Kommando-Zeile
Auf der Kommandozeile wird zunächst nochmal die Anzahl der Level und Iterationen geprintet und die dazu berechneten Werte zu number_of_unknowns und stepsize.
Number_of_unknowns beschreibt hierbei die Anzahl aller noch unbekannten Werte im Gitter, stepsize ist der Abstand der zwischen den Werten im Gitter liegt.
Nach jedem V-Cycle wird die diskrete L2-Norm des Residuums und der geschätzte Konvergenzfaktor ausgegeben und zum Schluss wird die Laufzeit der Berechnung und der totale Fehler der Annäherung geprintet

### Grid-Klasse
Zur Berechnung wurde eine eigene Klasse grid erstellt, die jeweils ein Gitter erstellen kann. Je nachdem, ob die Randwerte auch 0 sein sollen oder nicht, wird das Gitter dementsprechend initialisiert.
Die Klasse stellt Methoden bereit, um auf den Gittern die Restriktion und Interpolation zu berechnen.