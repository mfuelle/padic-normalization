# Affine Normalisierung über p-adischen Ganzheitsringen

`padic_normalization.py` ist ein SageMath-Modul. Es berechnet eine **endliche affine Präsentation der Normalisierung einschließlich der Einbettung**. Der Algorithmus rechnet exakt über den ganzen Zahlen, nicht mit gerundeten p-adischen Koeffizienten.

Die separate [Implementierungsdokumentation](IMPLEMENTIERUNG.md) erläutert die internen Datenstrukturen, den Grauert–Remmert-Schritt, die Singular-Anbindung, die vollständige Relationsberechnung und den Validierungsumfang.

## Installation und Aufruf

Die Dateien in dasselbe Verzeichnis legen. In einem Sage-Notebook, dessen Arbeitsverzeichnis dieses Verzeichnis ist:

```sage
from padic_normalization import normalize_affine_padic

R.<x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2 - 9*(x^3-x), p=3, verbose=True)
print(N)
N.relations               # exakte ganzzahlige Gleichungen
N.fraction_images         # Bilder aller Erzeuger im totalen Quotientenring
B = N.padic_ring(prec=40) # Sage-Quotientenring über Z_3
```

Alternativ `load("/voller/pfad/padic_normalization.py")` in Sage benutzen. Das Modul ist normales Python mit Sage-Imports; es benötigt keinen Sage-Präparser. Es läuft nicht in einer gewöhnlichen Python-Installation ohne Sage.

## Unterstützte Eingaben

### K = Q_p

* `f` liegt in einem Polynomring in genau zwei Variablen über `ZZ` oder `QQ`.
* Alle Koeffizienten sind p-integral. Rationale Nenner prim zu p sind erlaubt.
* Der p-adische Inhalt hat Bewertung höchstens eins; bei Bewertung eins wird die vertikale Komponente explizit abgespalten.
* Der primitive Teil von `f` ist quadratfrei über Q. Reduzible, aber reduzierte Ringe sind erlaubt; die Normalisierung kann ein Produkt sein und wird trotzdem als ein affiner Quotientenring präsentiert.

Diese Version verarbeitet also **exakt beschriebene reduzierte Modelle mit nichtkonstanter, von null verschiedener Gleichung**, nicht beliebige endlich genaue Elemente von `Qp(p)` oder `Zp(p)`. Solche Eingaben werden ausdrücklich zurückgewiesen. Nichtreduzierte Ringe werden zurückgewiesen. Eine reduzierte vertikale Komponente wird erhalten und als Produktfaktor ausgegeben.

### Erweiterungen K/Q_p

Optional unterstützt das Modul `extension=g`, wobei `g` ein normiertes Polynom in `ZZ[T]` vom Grad mindestens zwei ist, das **über Q_p irreduzibel** ist. Es muss weder Eisenstein noch modulo p irreduzibel sein. Insbesondere werden auch gemischt verzweigte Erweiterungen erfasst, wenn sie durch ein solches einzelnes Polynom beschrieben sind.

Die lokale Irreduzibilität wird exakt geprüft: Nach der Irreduzibilitätsprüfung über Q muss der Zahlkörper `F=QQ[T]/(g)` genau eine Primstelle über p besitzen. Eisenstein-Polynome und modulo p irreduzible Polynome behalten ihren bisherigen direkten Berechnungsweg.

Dann ist `f` ein Polynom in **drei** Variablen über `ZZ`, in der Reihenfolge `(a,x,y)`. Die erste Variable ist ein Koeffizientengenerator mit `g(a)=0`. Zusätzlich muss `f`, nach Reduktion modulo `g`, normiert in x oder y sein. Diese hinreichende Voraussetzung garantiert ein flaches globales Ausgangsmodell. Die generische Faser muss wiederum reduziert sein.

```sage
V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2 - 3*(x^3-x), 3, extension=T^2-3)
print(N)
print(N.fraction_images)
B = N.padic_ring(40)  # Präsentation über O_{Q_3(sqrt(3))}
```

Im allgemeinen Fall ist `Z_p[a]` möglicherweise nur eine echte Unterordnung von `O_K`. Das Modul berechnet deshalb eine vollständige Ganzheitsbasis von F und beginnt die Normalisierung über `O_F[x,y]/(f)`. Nach dem Koeffizientenwechsel nach Z_p ist dies das gewünschte Modell über O_K.

```sage
V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2-x^3, 2, extension=T^2-5)

O = N.coefficient_ring()        # exakte Präsentation des vollen O_F
B0 = N.ring_over_integers()     # Normalisierung als affine O_F-Algebra
A0, B0, phi = N.normalization_map()
B = N.padic_ring(prec=40)       # gemeinsame Z_2-Präsentation samt O_K-Koeffizientenrelationen
print(N.coefficient_model['basis'])
```

Hier liegt `(1+a)/2` im Ganzheitsring, aber nicht in `Z_2[a]`. Für das gemischt verzweigte Beispiel `g=T^4-14*T^2+9` kann derselbe Aufruf mit diesem Erweiterungspolynom verwendet werden; es hat bei p=2 Verzweigungsindex 2 und Restklassengrad 2.

Mit `integral_basis=True` lässt sich die Ganzheitsbasis-Präsentation auch bei den bisherigen Eisenstein- und unverzweigten Fällen ausdrücklich anfordern. Ohne diese Option bleibt ihre bisherige Ausgabeform erhalten. Konkrete Erzeuger und Standardbasen sind nicht kanonisch und können sich ändern.

Die Eingabesyntax für `f` bleibt auf ganzzahlige Koeffizienten in `(a,x,y)` beschränkt. Beliebige Koeffizienten aus O_K, etwa `(1+a)/2` unmittelbar in der Eingabegleichung, endlich genaue p-adische Koeffizienten und Eingaben als Erweiterungsturm werden nicht unterstützt. Eine bloße Irreduzibilität von `g` über Q reicht nicht: beispielsweise wird `T^2-17` bei p=2 zurückgewiesen, da es über Q_2 zerfällt.

## Bedeutung der Ausgabe

* `N.ambient_ring`, `N.ideal`, `N.relations`: exakte Präsentation über ZZ. Für Q_p interpretiert man dieselben Gleichungen über Z_p. Bei einer Erweiterung gehört die Gleichung `g(a)=0` zum ganzzahligen Modell; über O_K ersetzt man a durch den gewählten Koeffizientengenerator.
* `N.integer_ring()`: der exakte Sage-Quotientenring des globalen Modells.
* `A, B, phi = N.integer_map()`: exakte globale Modelle und Einbettung. Die Quelle verwendet weiterhin die ursprüngliche Ordnung `ZZ[a]`; bei nichtmaximaler Ordnung ist sie nicht das volle O_F-Ausgangsmodell. Die ursprünglichen Variablen bleiben erhalten.
* `N.coefficient_model`: im Ganzheitsbasis-Modus Daten zur vollen Koeffizientenordnung, einschließlich `basis`, `ramification_index` und `residue_degree`; andernfalls `None`.
* `N.coefficient_ring()` und `N.ring_over_integers()`: exakte Präsentationen von O_F und der normalisierten affinen O_F-Algebra. `coefficient_ring(prec=40)` gibt die Koeffizientenpräsentation über endlich genau dargestelltem Z_p zurück. `ring_over_integers(prec=40)` liefert die unten beschriebene p-adische Gesamtpräsentation in einem einzigen Polynomring.
* `A, B, phi = N.normalization_map()`: die exakten Gesamtmodelle über ZZ und die gesuchte Einbettung. Im Ganzheitsbasis-Modus enthalten beide Modelle den vollen Ganzheitsring als Koeffizientenvariablen; die Quelle beruht auf `N.integral_original_ideal`.
* `C, B, iota = N.coefficient_map(prec=40)`: Koeffizientenring, Gesamtpräsentation und die explizite Strukturabbildung von C nach B. Ohne `prec` sind alle drei Daten exakt über ZZ. Die definierenden Relationen werden vorab exakt geprüft. `iota` ist eine aufrufbare `CoefficientPresentationMap`, die auf Polynomrepräsentanten arbeitet, kein natives Sage-Homset-Element.
* `N.fraction_images`: Darstellung jedes neuen Erzeugers als Bruch in den ursprünglichen Variablen. Die Brüche sind **modulo den ursprünglichen Relationen** zu verstehen; sie müssen nicht im umgebenden rationalen Funktionenkörper gleich sein. Ihre Nenner sind Nichtnullteiler im Ausgangsring.
* `N.steps`: die sukzessiv adjungierten Brüche, jeweils im Variablensatz des vorhergehenden Schritts.
  Bei einer vertikalen Komponente enthält der letzte Eintrag stattdessen den Namen des trennenden Idempotenten.
* `N.certificate`: Abbruch durch das Grauert–Remmert-Kriterium, nicht bloß durch Glattheit der generischen Faser.
* `N.padic_ring(prec=40)`: im bisherigen direkten Weg ein Quotientenring über Sages p-adischem Ganzheitsring. Im Ganzheitsbasis-Modus eine **Gesamtpräsentation in einem einzigen Polynomring über Z_p**, mit a und den Ganzheitsbasiselementen als zusätzlichen Variablen samt Koeffizientenrelationen. Ihr `base_ring()` ist Z_p; die volle O_K-Algebrastruktur wird durch `coefficient_model` beziehungsweise `coefficient_ring(prec=40)` beschrieben. Sage unterstützt die Arithmetik über einem verschachtelten p-adischen Quotienten-Koeffizientenring nicht ausreichend. Es wird deshalb keine native O_K-DVR-Klasse vorgetäuscht. **Die exakte Spezifikation bleibt `N.relations`**; das p-adische Objekt ersetzt keinen exakten Gleichheitstest.

Das Programm minimiert die affine Präsentation nicht. Es kann mehrere zusätzliche Variablen und redundante Relationen ausgeben, auch wenn eine wesentlich kürzere isomorphe Präsentation existiert.

Auch über dem exakten, verschachtelten Koeffizientenring implementiert Sage nicht jede Quotientenoperation. Für Idealrechnungen und Gleichheitstests ist `integer_ring()` beziehungsweise `ideal` im gemeinsamen ganzzahligen Polynomring maßgeblich.

## Algorithmus und Korrektheit

Zunächst werden nur p-adische Einheiten aus rationalen Nennern und dem Inhalt von f entfernt. Dadurch erhält man ein exaktes ganzzahliges Modell A0. Für die bisherigen direkten Erweiterungsfälle ist A0 durch `(g(a), f(a,x,y))` gegeben. Im Ganzheitsbasis-Modus werden zusätzlich die Basisvariablen und die Multiplikationstabelle von O_F eingeführt. Damit beginnt die Schleife über dem vollen Ganzheitsring; die Daten der ursprünglichen Ordnung bleiben für die Bruchdarstellungen erhalten.

Für `f=p*f0` mit primitivem f0 wird zunächst f0 normalisiert. Die vollständige Normalisierung ist dann das Produkt mit `F_p[x,y]`. Ein zusätzlicher Idempotent e präsentiert dieses Produkt als affinen Ring; in den ursprünglichen Variablen gilt `e=p/(p+f0)`. Der nachfolgende iterative Algorithmus arbeitet mit dem flachen horizontalen Modell.

Ein Nichtnullteiler d wird aus dem relativen Jacobianideal gewählt. Im Erweiterungsfall wird zusätzlich mit g'(a) multipliziert. Die Quadratfreiheitsprüfung und ein ggT-Test auf der generischen Faser verifizieren die Nichtnullteilereigenschaft. Nach Invertieren von d ist A0 glatt über ZZ und deshalb normal.

Für den jeweils aktuellen Ring Ai werden berechnet:

```
Ji = radical(d*Ai)
Ci = (d*Ji : Ji)
Ai+1 = (1/d)*Ci = End_Ai(Ji).
```

`radicalZ` aus Singulars `primdecint.lib` berechnet das Radikal über ZZ. Der Idealquotient wird durch eine eigene Syzygienmatrix über ZZ bestimmt; das gewöhnliche `quotient` wird hier bewusst nicht verwendet. Die Singular-Bibliothek dokumentiert Einschränkungen dieser Funktion über ZZ.

Die neuen Erzeuger sind die vom Computer berechneten Elemente `u/d`. Ihre vollständigen Relationen werden durch Eliminieren einer Hilfsvariablen aus

```
(alte Relationen, d*z_i-u_i, 1-s*d)
```

bestimmt. Dies ist eine Sättigung nach d. Nur die Gleichungen `d*z_i-u_i` zu verwenden wäre falsch, weil dabei zusätzliche Komponenten auf d=0 entstehen können.

Wenn `Ci = d*Ai`, liefert das Grauert–Remmert-Kriterium Normalität: Ji ist radikal, enthält den Nichtnullteiler d und sein Verschwindungsort enthält den nichtnormalen Ort. Die letzte Eigenschaft bleibt erhalten, weil jeder Erweiterungsschritt nach Invertieren von d ein Isomorphismus ist. Ansonsten wächst Ai echt innerhalb der endlichen Normalisierung. Exzellenz und Noetherianität garantieren daher theoretische Terminierung.

Die Normalisierung wird zuerst über ZZ berechnet. Lokalisierung an p und der reguläre Koeffizientenwechsel nach Z_p erhalten die Normalität und liefern die Normalisierung des gewünschten Modells. Im Erweiterungsfall gilt `O_F tensor_Z Z_p = O_K`, weil genau eine Primstelle über p vorliegt. Es handelt sich um einen **Koeffizientenwechsel**, nicht um die adische Vervollständigung des ganzen Polynomrings.

## Laufzeit und Grenzen

Die mitgelieferte Testsuite wird mit `sage -python test_padic_normalization.py` gestartet. Die Datei `examples.sage` enthält direkt ausführbare Beispiele (`sage examples.sage`).

Die aktuelle Fassung besteht alle 21 Regressionstests unter SageMath 10.10. Die vollständigen lokalen Modul- und Testdateien wurden am 13. November 2026 über SageMathCell ausgeführt; SHA-256-Prüfsummen sichern die Zuordnung zur ausgelieferten Fassung. Geprüft wurden insbesondere nichtmaximale Ordnungen, gemischte Verzweigung, lokale Irreduzibilität und die neuen Ausgabeformen. Das Protokoll steht in `test_results.txt`.

Dies ist eine direkte Referenzimplementierung. Radikale und Eliminations-Gröbnerbasen über ZZ können sehr teuer sein. Sie normalisiert gegebenenfalls auch an anderen Primzahlen, deren Korrekturen nach dem Wechsel nach O_K redundant sind.

`max_steps=20` und `max_variables=24` sind anpassbare Ressourcenwächter. Beim Überschreiten wird `NormalizationLimitError` ausgelöst; ein Zwischenring wird nie als Normalisierung zurückgegeben. Einzelne Gröbnerbasis- oder Radikalberechnungen haben keinen eingebauten Zeitabbruch. Für große Beispiele ist ein separater Sage-Prozess mit äußerem Zeit-/Speicherlimit empfehlenswert.

Wie jedes CAS-Ergebnis setzt das Zertifikat die Korrektheit der verwendeten exakten Singular-Algorithmen voraus. Die Implementierung überprüft zusätzlich die Idealquotienten-Inklusion und die Graphrelationen. Die enthaltenen Regressionstests ersetzen keine unabhängige Begutachtung für Forschungsergebnisse.

## Quellen

* [Greuel–Laplagne–Seelisch: Normalization of Rings](https://arxiv.org/abs/0904.3561)
* [Singular: primdecint.lib](https://github.com/Singular/Singular/blob/spielwiese/Singular/LIB/primdecint.lib)
* [Sage: Singular-Schnittstelle](https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/singular.html)
* [Sage: Zahlkörper, Ganzheitsbasen und Primzerlegung](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html)
* [Stacks Project: Modelle und Vervollständigung](https://stacks.math.columbia.edu/tag/0C2T)
