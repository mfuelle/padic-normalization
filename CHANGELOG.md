# Änderungen vom 13. November 2026

- Allgemeine normierte Erweiterungspolynome über ZZ zugelassen, sofern sie über Q_p irreduzibel sind. Die Prüfung erfolgt exakt über die Primidealzerlegung im Zahlkörper. Dies umfasst insbesondere Erweiterungen von Q_2, deren Polynom weder Eisenstein noch modulo 2 irreduzibel ist.
- Vollständige Ganzheitsbasis und Multiplikationstabelle ergänzen bei Bedarf die nichtmaximale Ordnung ZZ[a]. Die Kurvennormalisierung beginnt dann über dem vollen Ganzheitsring.
- Neue Methoden `coefficient_ring`, `ring_over_integers`, `normalization_map` und `coefficient_map`; `integral_basis=True` fordert den neuen Ausgabeweg auch für die bisherigen Erweiterungstypen an.
- Der p-adische allgemeine Ausgabeweg verwendet eine Gesamtpräsentation über Z_p mit expliziter O_K-Struktur. Er benötigt keine native Sage-DVR-Klasse für das allgemeine Erweiterungspolynom. Einschränkungen von Sage bei Quotientenarithmetik sind in der README dokumentiert.
- Der Idealquotient wird durch sequenzielle Syzygienrechnungen berechnet. Vorher werden redundante Erzeuger modulo dem Relationsideal entfernt. Dies vermeidet die großen Blockmatrizen, die bei den neuen Testmodellen Singular-Abstürze verursachten.
- Die bisherigen Eingabebeschränkungen für die Kurvengleichung bleiben bestehen: ganzzahlige Koeffizienten in a,x,y, Normiertheit in einer Kurvenkoordinate und reduzierte generische Faser. Beliebige endlich genaue p-adische Koeffizienten werden weiterhin zurückgewiesen.
- Neue Regressionstests prüfen nichtmaximale Ordnungen, gemischte Verzweigung, die Normalisierung einer Spitze über dem vollen Ganzheitsring, Strukturabbildungen, lokale Irreduzibilität und Ressourcengrenzen. Den ausgeführten Teststand dokumentiert `test_results.txt`.

Die Benennung und Wahl der zusätzlichen Erzeuger ist nicht kanonisch. Beispielsweise kann bei `y^2=9(x^3-x)` nun `z=2y/3` ausgegeben werden. Die gewohnte Normalisierungskoordinate ist dann `t=y-z=y/3`; dieselbe Normalisierung wird präsentiert.
