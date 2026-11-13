# Implementierungsdokumentation: affine Normalisierung über p-adischen Ganzheitsringen

Diese Dokumentation beschreibt den Aufbau von [`padic_normalization.py`](padic_normalization.py), seine mathematischen Invarianten und die Umsetzung in SageMath und Singular. Die [`README.md`](README.md) enthält die Benutzungsanleitung; hier steht im Vordergrund, wie die Berechnung funktioniert und welche Aussagen die Ausgabe rechtfertigt.

## 1. Ziel und Darstellungsprinzip

Gesucht ist die Normalisierung eines reduzierten Rings

$$
A=\mathcal O_K[x,y]/(f)
$$

in seinem totalen Quotientenring. Bei einem Integritätsbereich ist dies der gewöhnliche Quotientenkörper; bei mehreren Komponenten erlaubt die Implementierung auch eine Normalisierung, die ein Produkt von Ringen ist.

Das Ergebnis besteht aus einer endlichen Präsentation

$$
\overline A\cong
\mathcal O_K[x,y,z_1,\ldots,z_m]/I
$$

und den Bildern der zusätzlichen Erzeuger als Brüche in den ursprünglichen Variablen. Die ursprünglichen Variablen bleiben erhalten, sodass die Normalisierungseinbettung explizit ist.

**Die Normalisierungsrechnung verwendet ausschließlich exakte Arithmetik.** Sie findet über `ZZ` statt. Erst die optionale Methode `padic_ring(prec=40)` überträgt das Ergebnis in einen Sage-Ring mit endlich genauer p-adischer Koeffizientendarstellung.

Die Implementierung ist daher kein Algorithmus für beliebige endlich genau bekannte p-adische Koeffizienten. Sie verarbeitet die unten beschriebenen exakten Modelle. Sie minimiert weder die Anzahl der Erzeuger noch die ausgegebenen Relationen.

## 2. Dateien und Abhängigkeiten

| Datei | Aufgabe |
|---|---|
| `padic_normalization.py` | Öffentliche API, Vorverarbeitung, Normalisierungsschleife und Singular-Anbindung |
| `examples.sage` | Aufrufbeispiele einschließlich der verzweigten quadratischen Erweiterung |
| `test_padic_normalization.py` | Algebraische Regressionstests mit `unittest` |
| `README.md` | Installation, Benutzung und mathematischer Überblick |
| `test_results.txt` | Protokoll und genaue Einordnung der bisherigen Validierung |
| `IMPLEMENTIERUNG.md` | Diese interne Dokumentation |

Benötigt werden SageMath und dessen Singular-Schnittstelle `sage.interfaces.singular.singular`. Das Backend lädt Singulars Bibliothek `primdecint.lib` für `radicalZ`. Die Python-Datei braucht keinen Sage-Präparser, aber eine Python-Umgebung mit Sage; ein gewöhnliches CPython ohne Sage reicht nicht aus.

## 3. Vertrag der öffentlichen Funktion

```python
normalize_affine_padic(
    f, p, extension=None,
    max_steps=20, max_variables=24, verbose=False, integral_basis=False
)
```

Der Rückgabewert ist ein `PadicNormalization`-Objekt. Nicht erfüllte Voraussetzungen führen zu einer Exception; die Funktion gibt keine lediglich vermutete Normalisierung zurück.

### 3.1 Grundfall: K = Q_p

Die Eingabe `f` muss ein von null verschiedenes, nichtkonstantes Polynom in genau zwei Variablen über `ZZ` oder `QQ` sein. Die Variablennamen sind frei wählbar. Es gelten folgende Bedingungen:

- `p` ist eine Primzahl.
- Jeder rationale Koeffizient ist p-integral.
- Der primitive Teil von `f` ist über `QQ` quadratfrei.
- Die minimale p-adische Bewertung der von null verschiedenen Koeffizienten ist null oder eins.

Bewertung null führt zur eigentlichen Normalisierungsschleife. Bewertung eins wird durch Abspalten einer vertikalen Komponente behandelt, siehe Abschnitt 8. Bewertung mindestens zwei und wiederholte Faktoren werden als nichtreduzierte Eingaben zurückgewiesen. Die Reduktion eines nichtreduzierten Rings wird nicht stillschweigend vorgenommen.

### 3.2 Erweiterungsfall

Für `extension=g` wird ein normiertes Polynom `g` in `ZZ[T]` vom Grad mindestens zwei verlangt, das über `QQ_p` irreduzibel ist. Auch Polynome, die weder Eisenstein noch modulo p irreduzibel sind, sind zugelassen. Der Zahlkörper F und seine eindeutige p-adische Vervollständigung K erfüllen

$$
F=\mathbb Q[a],\quad g(a)=0,\qquad
\mathcal O_K=\mathcal O_F\otimes_{\mathbb Z}\mathbb Z_p.
$$

`f` liegt dann in einem Polynomring über `ZZ` mit genau drei Variablen. Ihre Reihenfolge hat eine feste Bedeutung: zuerst der Koeffizientengenerator `a`, danach die beiden Kurvenkoordinaten. Rationale Koeffizienten mit nichttrivialen Nennern sind in diesem Modus nicht zugelassen.

Nach Einsetzen von `a` in den Zahlkörper `QQ[T]/(g)` muss `f` nichtkonstant, quadratfrei und in mindestens einer Kurvenkoordinate normiert sein. Diese Normiertheit ist eine hinreichende Bedingung für die benötigte globale Flachheit; sie ist keine notwendige Bedingung für die Existenz einer Normalisierung.

Eine Darstellung durch einen Erweiterungsturm ist weiterhin nicht vorgesehen. Gemischt verzweigte Körper werden jedoch über ein einzelnes geeignetes Polynom unterstützt. Beliebige Koeffizienten aus O_K können noch nicht direkt in `f` eingesetzt werden: die Eingabe bleibt ein ganzzahliges Polynom in `a,x,y`.

Eisenstein-Polynome und modulo p irreduzible Polynome verwenden standardmäßig den bisherigen direkten Weg. Für alle anderen Polynome wird der volle Ganzheitsring anhand einer Ganzheitsbasis präsentiert. `integral_basis=True` fordert diese Präsentation auch für die bisherigen Fälle an.

## 4. Vorverarbeitung und Testelement

### 4.1 Exaktes ganzzahliges Modell

`_exact_polynomial` prüft Variablenzahl, Koeffizientenring, Nichtkonstanz und p-Integralität. Insbesondere weist es Polynome über endlich genauen `Qp`- oder `Zp`-Ringen zurück.

Im Grundfall beseitigt `_prepare` rationale Nenner und teilt durch den ganzzahligen Inhalt. Im flachen Fall sind diese Faktoren Einheiten in `Z_p`. Die Gleichung ändert daher den gewünschten Ring über `Z_p` nicht. Das entstehende primitive Polynom `f_0` definiert ein torsionsfreies, also flaches Modell

$$
A_0=\mathbb Z[x,y]/(f_0).
$$

Die Faktorisierung über `QQ` prüft die Reduziertheit der generischen Faser. Zusammen mit der Torsionsfreiheit folgt daraus die Reduziertheit von `A_0`.

Im direkten Erweiterungsfall ist das ganzzahlige Modell

$$
A_0=\mathbb Z[a,x,y]/(g(a),f_0).
$$

Die Normiertheit von `g` und von `f_0` in einer Kurvenkoordinate liefert die Torsionsfreiheit über `ZZ`. Potenzen von `a` werden exakt modulo `g(a)` reduziert.

### 4.2 Lokale Irreduzibilität und Ganzheitsbasis

`_prepare` prüft zunächst die Irreduzibilität von `g` über Q. Wenn weder das Eisenstein-Kriterium noch die Irreduzibilität modulo p greifen, berechnet es anschließend `F.primes_above(p)`. Es muss genau eine Primstelle P geben und `e(P/p)*f(P/p)=deg(g)` gelten. Dies ist eine exakte Prüfung, dass `F tensor_Q Q_p` ein Körper ist. Ein Polynom wie `T^2-17`, das über Q irreduzibel ist, über Q_2 aber zerfällt, wird zurückgewiesen.

`_integral_coefficient_model` berechnet mit `F.integral_basis()` eine vollständige Z-Basis `b_0,...,b_(n-1)` von O_F. Es wird nicht angenommen, dass der erste Basisvektor eins ist. Rationale lineare Algebra liefert ganzzahlige Koordinaten für eins, a und jedes Produkt zweier Basisvektoren:

$$
1=\sum_k u_kb_k,\qquad a=\sum_k v_kb_k,
\qquad b_i b_j=\sum_k c_{ijk}b_k.
$$

Daraus entsteht der Koeffizientenring

$$
C=\mathbb Z[a,B_0,\ldots,B_{n-1}]/
\left(1-\sum_k u_kB_k,\ a-\sum_k v_kB_k,
\ B_iB_j-\sum_k c_{ijk}B_k\right).
$$

Dieser Ring ist exakt O_F: Die Relationen reduzieren jedes Polynom auf eine Z-Linearkombination der Basisvariablen. Die Auswertung in O_F ist surjektiv, und die lineare Unabhängigkeit der Ganzheitsbasis beweist Injektivität.

Die ursprünglichen Variablen `a,x,y` bleiben erhalten. Die zusätzlichen Variablen heißen `ib0`, `ib1`, …; vorhandene Namen werden übersprungen. Die Multiplikationsrelationen und die ursprüngliche Kurvengleichung werden in einem einzigen ganzzahligen Polynomring zusammengeführt. Die Normalisierungsschleife startet somit bei `O_F[x,y]/(f)`. Die Normiertheit von f garantiert auch hier Flachheit und die Einbettung des ursprünglichen Ordnungsmodells.

Für die Bruchdarstellungen werden alle Basisvektoren als rationale Polynome in a gespeichert. Ihre ganzzahligen Nenner sind Nichtnullteiler im flachen Ausgangsmodell. Die Basisergänzung ist Vorverarbeitung; sie erscheint nicht als Grauert–Remmert-Schritt in `steps`, zählt aber gegen `max_variables`.

### 4.3 Wahl von d

Die Schleife benötigt einen Nichtnullteiler `d`, sodass `A_0[1/d]` normal ist. Die Implementierung wählt ihn aus einem Jacobischen Ideal, sodass diese Lokalisierung sogar glatt über `ZZ` ist.

Im Grundfall werden die von null verschiedenen Ableitungen von `f_0` nach ihrem Totalgrad sortiert. Für kleine ganze Zahlen `c` werden Linearkombinationen getestet. Bei zwei Ableitungen haben sie die Form

$$
d=h_0+c h_1.
$$

Akzeptiert wird eine Kombination mit `gcd(f_0,d)=1` über `QQ`. Die Quadratfreiheit in Charakteristik null stellt sicher, dass keine irreduzible Komponente beide partiellen Ableitungen teilt. Die im Code verwendete endliche Suche findet deshalb eine zulässige Kombination. Die Torsionsfreiheit über `ZZ` überträgt die Nichtnullteilereigenschaft auf das ganzzahlige Modell.

Im Erweiterungsfall erfolgt die entsprechende ggT-Prüfung über `QQ[T]/(g)`. Für eine zulässige Kombination `h` der Ableitungen nach den Kurvenkoordinaten wird

$$
d=g'(a)h
$$

verwendet. Dies ist eine Linearkombination der passenden Jacobischen Minoren der beiden Gleichungen `g(a), f_0`. Der Faktor `g'(a)` berücksichtigt im direkten Weg auch mögliche Nichtnormalität des globalen Koeffizientenmodells an anderen Primstellen. Im Ganzheitsbasis-Modus bleibt dasselbe d zulässig: Nach Invertieren von `g'(a)` ist die ursprüngliche Ordnung bereits étale über ZZ und stimmt mit ihrem Ganzheitsabschluss überein. Das zusätzliche Invertieren von h liefert weiterhin einen glatten, also normalen offenen Teil.

Während der gesamten Iteration bleibt `d` das Bild dieses ursprünglichen Elements. Es wird nur in die jeweils größere Präsentation übertragen.

## 5. Normalisierungsschleife

### 5.1 Mathematischer Schritt

Sei der aktuelle Ring `A_i` reduziert und `d` der gewählte Nichtnullteiler. Setze

$$
J_i=\sqrt{dA_i}.
$$

Dann enthält `J_i` einen Nichtnullteiler, und seine Nullstellenmenge enthält den Nichtnormalitätsort. Das Grauert–Remmert-Kriterium liefert

$$
A_i\text{ normal}
\quad\Longleftrightarrow\quad
\operatorname{End}_{A_i}(J_i)=A_i.
$$

Als Unterring des totalen Quotientenrings lässt sich der Endomorphismenring schreiben als

$$
A_{i+1}=\operatorname{End}_{A_i}(J_i)
=\frac{1}{d}(dJ_i:_{A_i}J_i).
$$

Diese Erweiterung ist endlich, integral und birational. Wenn sie echt größer ist, werden Erzeuger `u_j/d` adjungiert und sämtliche Relationen berechnet. Ist sie nicht größer, endet die Schleife mit dem Normalitätsnachweis.

### 5.2 Übersetzung in Ideale eines Polynomrings

Der Code arbeitet nicht direkt mit Idealen eines Quotientenrings. Er hält eine Präsentation `A_i=R/I` mit einem Polynomring `R` über `ZZ` und rechnet mit Urbildern in `R`:

| Codevariable | Bedeutung im umgebenden Polynomring |
|---|---|
| `I` | Definierendes Ideal von `A_i` |
| `J` | `sqrt(I + (d))`, Urbild von `J_i` |
| `H` | `I + d*J`, Urbild von `d*J_i` |
| `C` | `(H:J)`, Urbild von `(d*J_i:J_i)` |
| `old` | Anfangs `I + (d)`, anschließend um bereits ausgewählte Zähler erweitert |

Die Erzeuger von `C` werden nacheinander modulo `old` reduziert. Ein von null verschiedener Rest wird als neuer Zähler gespeichert und zu `old` hinzugefügt. Dies vermeidet Erzeuger, die durch frühere Zähler bereits erfasst sind; eine minimale Erzeugermenge wird dadurch nicht garantiert.

Bleibt die Zählerliste leer, gilt

$$
C=I+(d),
$$

also das Abbruchkriterium. Andernfalls liefert die Liste Erzeuger der nächsten endlichen Algebra. `max_steps` wird erst geprüft, nachdem festgestellt wurde, dass eine weitere Vergrößerung erforderlich ist.

### 5.3 Invarianten und Terminierung

Nach jedem erfolgreichen Schritt gilt:

1. Die aktuelle Algebra ist eine endliche Unteralgebra der Normalisierung des Ausgangsmodells.
2. Alle neuen Elemente sind als Brüche mit Nichtnullteilernennern im ursprünglichen Ring beschrieben.
3. Nach Invertieren von `d` ist die Algebra unverändert und normal.
4. Die Präsentation enthält den vollständigen Relationskern der adjungierten Brüche.

Endlich erzeugte `ZZ`-Algebren sind exzellent; ihre Normalisierung ist im reduzierten Fall endlich. Die aufsteigende Folge liegt damit in einem endlichen Modul über dem noetherschen Ausgangsring und stabilisiert. Dies ist die theoretische Terminierungsaussage. Sie liefert keine brauchbare allgemeine Laufzeitschranke für die Gröbnerbasisrechnungen.

## 6. Singular-Backend: `_Engine`

### 6.1 Standardbasen und Radikale

`std` verwendet den eingebetteten Singular-Helfer `pn_std`. Dieser berechnet eine Standardbasis über `ZZ`, vereinheitlicht die Vorzeichen der führenden Koeffizienten und markiert das Ergebnis als Standardbasis.

`radical` ruft `radicalZ` aus `primdecint.lib` auf. Es wird tatsächlich das Radikal im ganzzahligen Polynomring berechnet. Ein Radikal nur über `QQ` würde vertikale Nichtnormalität nicht erfassen.

`nf` berechnet Normalformen über Singulars `reduce`; `subset` prüft Idealinklusionen durch Normalformen der Erzeuger.

### 6.2 Idealquotient aus Syzygien

`pn_colon(H,J)` implementiert `(H:J)` mit einer Modul-Syzygienrechnung über `ZZ`. Schreibe

$$
H=(h_1,\ldots,h_m),\qquad J=(j_1,\ldots,j_r).
$$

Zunächst reduziert Python die Erzeuger von J modulo I und entfernt Nullreste. Da I in H liegt, ändern diese Schritte den Idealquotienten nicht. H wird vor der Syzygienrechnung durch eine Standardbasis ersetzt.

Der Helfer beginnt mit K=(1) und verarbeitet die Erzeuger von J einzeln. Für K=(k_1,...,k_s) berechnet er die Syzygien der Liste `(j_i*k_1,...,j_i*k_s,h_1,...,h_m)`. Eine Syzygie erfüllt

$$
j_i\sum_\ell u_\ell k_\ell+\sum_k v_k h_k=0.
$$

Die aus diesen Syzygien gebildeten Elemente `sum(u_l*k_l)` erzeugen genau `K intersect (H:j_i)`. Nach jedem Durchgang wird K durch die Standardbasis dieses Ideals ersetzt. Am Ende ist K=(H:J). Für das Nullideal J bleibt K=(1). Dieses sequenzielle Verfahren vermeidet die große Blockmatrix, die beim gleichzeitigen Erzwingen aller Bedingungen entstehen würde.

Anschließend prüft Python zusätzlich `K*J ⊆ H` für das berechnete Ideal `K`. Diese Prüfung erkennt bestimmte Backend- oder Konvertierungsfehler. Sie prüft nur diese Inklusion; die Vollständigkeit des Idealquotienten beruht auf der korrekten Syzygienrechnung.

### 6.3 Vollständige Relationen durch Sättigung

Werden die Brüche `u_1/d,...,u_s/d` adjungiert, reichen die Gleichungen

$$
dz_j-u_j=0
$$

allein im Allgemeinen nicht aus: Sie können zusätzliche Komponenten über `d=0` erzeugen. Der korrekte Relationskern ist

$$
I_{\mathrm{neu}}=(I,dz_1-u_1,\ldots,dz_s-u_s):d^\infty.
$$

`graph` berechnet diesen Kern durch Einführung einer Hilfsvariablen `t` und Elimination:

$$
I_{\mathrm{neu}}=
(I,dz_1-u_1,\ldots,dz_s-u_s,1-td)
\cap R[z_1,\ldots,z_s].
$$

Der verwendete Block-Termorder stellt `t` in einen vorderen lexikographischen Block. Aus der Standardbasis werden die Polynome ausgewählt, die nicht von `t` abhängen. Nach Entfernung der Hilfsvariablen wird erneut eine Standardbasis berechnet.

Die zusätzliche Python-Prüfung stellt sicher, dass die alten Relationen und alle Gleichungen `dz_j-u_j` im neuen Ideal liegen. Dass der gesamte Kern berechnet wurde, folgt aus der Eliminationsrechnung; die Inklusionsprüfung allein wäre dafür kein Nachweis.

### 6.4 Schnittstellendetails

- Vor jedem Singular-Aufruf wird der passende Ring als aktueller Ring gesetzt.
- Ergebnisse werden mit `sage_poly(R)` ausdrücklich in den vorgesehenen Sage-Ring zurückkonvertiert.
- Die Helfer werden beim Erzeugen einer `_Engine` neu geladen, auch nach einem Neustart des Singular-Prozesses.
- Der Helfertext wird als eine einzige Zeile an `singular.eval` übergeben. Dies vermeidet, dass die Expect-Schnittstelle bei einer noch unvollständigen Prozedurdeklaration wartet. Im Helfertext dürfen deshalb keine `//`-Zeilenkommentare ergänzt werden.
- Die Eliminationsvariable heißt zunächst `pn_inverse`; bei Namenskollisionen werden Unterstriche angehängt. Neue Algebraerzeuger heißen `z1`, `z2`, …; belegte Namen werden übersprungen.
- Die gemeinsame Singular-Sitzung besitzt veränderlichen Zustand. Parallele Aufrufe in derselben Sitzung sind nicht abgesichert; für parallele Rechnungen sollten getrennte Prozesse verwendet werden.

## 7. Ergebnisobjekt und Abbildungen

### 7.1 Attribute

| Attribut | Inhalt |
|---|---|
| `p` | Primzahl als Sage-Ganzzahl |
| `original_ideal` | Ideal des vorbereiteten ganzzahligen Ausgangsmodells |
| `integral_original_ideal` | Ausgangsideal nach eventueller Ergänzung der Ganzheitsbasis |
| `coefficient_model` | Daten des vollen Koeffizientenrings im Ganzheitsbasis-Modus, sonst `None` |
| `ambient_ring` | Ganzzahliger Polynomring der abschließenden Präsentation |
| `ideal` | Abschließendes Relationsideal |
| `relations` | Tupel der Erzeuger von `ideal` |
| `fraction_images` | Dictionary: Variablenname → Bruch in den ursprünglichen Variablen |
| `steps` | Tupel mit den Daten der Vergrößerungsschritte |
| `extension` | Exaktes Polynom `g` oder `None` |
| `certificate` | Text zur erreichten mathematischen Abbruchbedingung |

`original_ideal` kann von der wörtlich eingegebenen Gleichung abweichen, weil Nenner und am Ort `p` invertierbare Inhaltsfaktoren entfernt werden. Im Erweiterungsfall bezeichnet es weiterhin das Modell über der ursprünglichen Ordnung `ZZ[a]`. `integral_original_ideal` enthält im Ganzheitsbasis-Modus dagegen das tatsächliche Ausgangsmodell über O_F.

`coefficient_model` enthält `ambient_ring`, `ideal` und `names` der Koeffizientenpräsentation sowie `number_field`, `basis`, `ramification_index` und `residue_degree`. Alle Basisdaten und Relationskoeffizienten sind exakt.

Ein gewöhnlicher Eintrag in `steps` hat die Form

```python
{
    'variables': ('z1', ...),
    'numerators': (u1, ...),
    'denominator': d,
}
```

Die Polynome dieses Eintrags gehören zum Polynomring **vor** dem jeweiligen Schritt. Bei einer vertikalen Komponente kommt abschließend ein Eintrag `{'vertical_idempotent': name}` hinzu. Deshalb ist `len(steps)` dort nicht die reine Zahl der Grauert–Remmert-Vergrößerungen.

### 7.2 Bruchdarstellungen

Zur Buchführung verwendet der Code den rationalen Funktionenkörper des ursprünglichen umgebenden Polynomrings. Neue Brüche entstehen durch Einsetzen der schon bekannten Erzeugerbilder. Dieser Hilfskörper ist nicht mit dem totalen Quotientenring des Ausgangsrings gleichzusetzen.

Die Ausdrücke in `fraction_images` sind als Bruchrepräsentanten **modulo `original_ideal`** zu verstehen. Beispielsweise stellen bei `y^2=x^3` die Ausdrücke `x^2/y` und `y/x` dasselbe Element dar, obwohl sie im freien rationalen Funktionenkörper verschieden sind. Ein gewöhnlicher Vergleich der beiden Python-Bruchobjekte prüft diese Gleichheit nicht.

Zur Prüfung einer solchen Gleichheit kann man nach Beseitigung der Nichtnullteilernenner den Zähler der Differenz modulo `original_ideal` reduzieren. Der Code erzeugt die Brüche konstruktiv mit zulässigen Nennern.

### 7.3 Methoden

`integer_ring()` konstruiert den exakten Sage-Quotientenring `ambient_ring/ideal`.

`integer_map()` liefert `(A0, B0, phi)`: das vorbereitete ganzzahlige Ordnungsmodell, seine Normalisierung und die Einbettung. Bei nichtmaximaler ursprünglicher Ordnung ist A0 nicht das volle O_F-Modell. Da die ursprünglichen Variablen am Anfang des Variablensatzes bleiben, werden ihre Restklassen auf die gleichnamigen Restklassen in `B0` abgebildet.

`coefficient_ring()` liefert im Ganzheitsbasis-Modus den exakten Quotientenring C aus Abschnitt 4.2. Mit `prec=N` werden die Koeffizienten der Multiplikationstabelle nach `Zp(p,N)` übertragen. Dies ist eine Präsentation von O_K, keine native p-adische DVR-Klasse mit Bewertungs- oder Uniformisierer-API. Ohne Ganzheitsbasis-Modus verlangt die Methode einen erneuten Aufruf mit `integral_basis=True`.

`ring_over_integers()` präsentiert die Normalisierung als affine Algebra über dem exakten Koeffizientenring C. `_relative_ring` verschiebt dafür die Koeffizientenvariablen a und die Ganzheitsbasisvariablen in den Basisring. Zuvor werden die Relationen im ganzzahligen Polynomring exakt modulo dem Koeffizientenideal reduziert; exakte Nullrelationen werden entfernt.

`normalization_map()` liefert die exakten Modelle `(A,B,phi)` in gemeinsamen Polynomringen über ZZ, wobei die Quelle aus `integral_original_ideal` konstruiert wird. Im Ganzheitsbasis-Modus enthalten beide Ringe die Koeffizientenvariablen und definieren somit O_F-Algebren; ihre Sage-Basisringe sind ZZ. Dies ist die Einbettung mit vollem Ganzheitsring in der Quelle. Ihr Koeffizientenwechsel nach Z_p ergibt die gesuchte Normalisierungseinbettung über O_K.

`coefficient_map(prec=None)` liefert `(C,B,iota)` mit dem vollen Koeffizientenring, der Gesamtpräsentation in einem einzigen Polynomring und der Strukturabbildung. Ohne Präzisionsargument ist B das exakte `integer_ring()`; mit `prec=N` werden C und B über `Zp(p,N)` präsentiert. Die Bilder der Koeffizientenerzeuger sind die entsprechenden Variablen von B.

Beide Abbildungsmethoden prüfen die erforderlichen Relationsinklusionen im exakten ganzzahligen Gesamtmodell. `normalization_map()` erzeugt anschließend einen gewöhnlichen Sage-Ringhomomorphismus zwischen den ganzzahligen Gesamtmodellen. `coefficient_map()` verwendet dagegen die kleine aufrufbare Klasse `CoefficientPresentationMap` mit `domain()`, `codomain()` und `im_gens()`. Diese wertet Polynomrepräsentanten aus und erzeugt das Ergebnis ohne Quotientenreduktion. Sie ist kein Element eines nativen Sage-Homsets. Dies vermeidet die Normalformrechnungen, die Sage auch beim Erzeugen von Homomorphismen mit deaktivierter Relationsprüfung verlangen kann. Weitergehende Rechnungen mit verschachtelten oder p-adischen Quotienten können dennoch an Sage-Grenzen stoßen; das gemeinsame ganzzahlige Modell bleibt dafür die Referenz.

`padic_ring(prec=40)` verwendet im bisherigen direkten Erweiterungsweg `Zp(p,prec).extension(g)` und ersetzt a durch den Generator des Koeffizientenrings. Im Ganzheitsbasis-Modus liefert es dagegen wie `ring_over_integers(prec=40)` eine **Gesamtpräsentation in einem einzigen Polynomring** `Zp[a,x,y,ib0,...,z1,...]/I` mit allen Koeffizientenrelationen. Deren Sage-`base_ring()` ist Z_p. Die zusätzliche O_K-Algebrastruktur ist durch `coefficient_model` und `coefficient_ring(prec=40)` ausdrücklich beschrieben. Dies umgeht fehlende Sage-Arithmetik über verschachtelten p-adischen Quotientenringen.

Im direkten Erweiterungsweg reduziert die Methode die Relationen vor dem endlich genauen Koeffizientenwechsel exakt modulo `g(a)` und entfernt exakt verschwindende Polynome. So wird die exakte Nullrelation `g(a)=0` nicht als bloßer Ausdruck `O(p^N)` in das Relationsideal übernommen. In dieser Gesamtpräsentation des Ganzheitsbasis-Modells bleiben a und die Basisvariablen Polynomvariablen, sodass deren Gleichungen als Relationen erhalten bleiben. Es wird kein exakter Nulltest aus endlich genauen p-adischen Daten abgeleitet.

Das zurückgegebene p-adische Objekt ist ein Sage-Quotientenring. Die Verfügbarkeit weiterer Rechenoperationen, insbesondere Gröbnerbasisoperationen über diesem Koeffizientenring, ist eingeschränkt. Die exakte Spezifikation und der Normalitätsnachweis bleiben beim ganzzahligen Ergebnisobjekt.

`certificate` ist eine Beschreibung der erreichten Bedingung, kein separat serialisiertes oder unabhängig formal verifiziertes Beweisobjekt. Die Korrektheit setzt die verwendeten exakten Sage-/Singular-Operationen voraus.

## 8. Vertikale Komponente im Grundfall

Bei p-adischem Inhalt eins hat die Eingabe, bis auf einen invertierbaren Faktor, die Form `p*f_0`, wobei `f_0` primitiv ist. Dann gilt

$$
\overline{\mathbb Z_p[x,y]/(pf_0)}
\cong
\overline{\mathbb Z_p[x,y]/(f_0)}
\times\mathbb F_p[x,y].
$$

Die Funktion normalisiert zunächst den horizontalen Faktor rekursiv. `_with_vertical_component` stellt anschließend das Produkt als einen einzigen affinen Quotientenring dar.

Seien `H_l=0` seine horizontalen Relationen und `z_j` die zusätzlich eingeführten horizontalen Erzeuger. Die Produktpräsentation verwendet einen neuen Erzeuger `e` mit

$$
e^2=e,\qquad p(1-e)=0,\qquad eH_l=0,\qquad (1-e)z_j=0.
$$

Die ursprünglichen Variablen `x,y` bleiben auf beiden Faktoren erhalten. Auf dem horizontalen Faktor ist `e=1`, auf dem vertikalen Faktor `e=0`. Der Name lautet zunächst `vertical_split` und wird bei Kollisionen angepasst.

Die Einbettung wird ebenfalls explizit beschrieben:

$$
e=\frac{p}{p+f_0}.
$$

Ein horizontaler Erzeuger mit Bruchdarstellung `u/v` wird durch

$$
\frac{pu}{pv+f_0}
$$

ersetzt. Dieser Bruch ist horizontal gleich `u/v` und vertikal gleich null. Seine Nenner sind auf allen Komponenten Nichtnullteiler: horizontal wirkt `p*v`, vertikal das von null verschiedene Polynom `f_0 mod p`.

Diese Zusatzbehandlung gilt nur für den Grundfall `K=Q_p`. Der Erweiterungsmodus verlangt das oben beschriebene flache Modell.

## 9. Warum die ganzzahlige Rechnung das p-adische Problem löst

Die Implementierung normalisiert ein ganzzahliges Modell global. Damit können auch Erzeuger entstehen, die nur an Primstellen ungleich `p` benötigt werden. Sie werden nach dem Koeffizientenwechsel redundant, verfälschen das Ergebnis aber nicht.

Der Koeffizientenwechsel `ZZ → Z_p` ist flach und regulär. Normalität bleibt unter diesem regulären Basiswechsel erhalten. Die endliche birationale Normalisierungseinbettung bleibt endlich und ist nach Invertieren des Nichtnullteilers `d` weiterhin ein Isomorphismus. Damit ist die erhaltene normale Algebra die Normalisierung des p-adischen Ausgangsmodells.

Im direkten Erweiterungsweg ist `Z_p[a]/(g(a))` wegen der Eisenstein- beziehungsweise Irreduzibilitätsbedingung modulo p bereits O_K. Im allgemeinen Weg wird stattdessen der volle Ganzheitsring O_F verwendet. Die exakt geprüfte Eindeutigkeit der Primstelle über p ergibt `O_F tensor_Z Z_p = O_K`. Dieselbe Basiswechselargumentation liefert damit auch bei nichtmaximaler ursprünglicher Ordnung und gemischter Verzweigung die gewünschte O_K-Algebra.

Gemeint ist hier der Wechsel der **Koeffizienten** in einem Polynomring. Es wird weder ein Potenzreihenring noch eine Tate-Algebra konstruiert. Auch die vollständige p-adische Vervollständigung des gesamten affinen Rings ist nicht das Ziel der Funktion.

## 10. Beispiel: Q_3(sqrt(3))

```sage
from padic_normalization import normalize_affine_padic

V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(
    y^2 - 3*(x^3-x), p=3, extension=T^2-3
)

print(N.relations)
print(N.fraction_images)
A0, B0, phi = N.integer_map()
B = N.padic_ring(prec=40)
```

Im dokumentierten Rechenlauf wurde ein Erzeuger `z1=y/a` gefunden. Die ausgegebene ganzzahlige Präsentation hatte die Relationen

$$
az_1-y,\quad y^2-3z_1^2,\quad ay-3z_1,
\quad a^2-3,\quad x^3-z_1^2-x.
$$

Nach Interpretation von `a` als Koeffizient und Elimination von `y` erhält man die kürzere Darstellung

$$
\overline A\cong\mathcal O_K[x,z]/(z^2-x^3+x),
\qquad y\longmapsto az.
$$

Diese Vereinfachung erläutert die Ausgabe; eine automatische Minimierung ist nicht Teil des Moduls. Erzeugernamen, Vorzeichen und Reihenfolge einer Standardbasis sollten nicht als kanonisch behandelt werden.

## 11. Fehlerbehandlung und Rechenaufwand

| Exception | Bedeutung |
|---|---|
| `TypeError` | Ungeeigneter Koeffizientenring oder nicht unterstützte Erweiterungsbeschreibung |
| `ValueError` | Verletzte Eingabevoraussetzung, etwa Nichtreduziertheit, nichtprime Zahl oder fehlende Normiertheit |
| `NormalizationLimitError` | Eine weitere Vergrößerung würde eine konfigurierte Ressourcengrenze überschreiten |
| `ArithmeticError` | Testelement konnte nicht konstruiert werden oder eine interne Konsistenzprüfung schlug fehl |

Weitere Exceptions des Sage-/Singular-Backends werden nicht generell abgefangen.

`max_steps` begrenzt die Zahl der eigentlichen Vergrößerungsschritte. Auch bei `max_steps=0` wird geprüft, ob das Ausgangsmodell bereits normal ist. `max_variables` wird vor dem Ergänzen der Ganzheitsbasis, vor dem Adjungieren neuer Erzeuger und beim Abspalten der vertikalen Komponente geprüft. Es ist keine allgemeine Speichergrenze und zählt insbesondere die vorübergehende Eliminationsvariable nicht mit.

Es gibt keinen Zeitabbruch innerhalb einer einzelnen Radikal-, Syzygien- oder Standardbasisrechnung. Für feste Laufzeit- oder Speichergrenzen ist eine Ausführung in einem separaten Prozess erforderlich.

Die dominierenden Kosten entstehen durch `radicalZ`, Modulsyzygien und die Eliminationsbasis der Graphpräsentation. Zahlreiche zusätzliche Variablen, hohe Grade und globale Beiträge anderer Primstellen können die Berechnung stark vergrößern. Die gegebenen Grenzen sind Schutzschranken, keine Zusage, dass eine Eingabe darunter schnell bearbeitet wird.

## 12. Validierung und Wartung

Die Tests prüfen algebraische Eigenschaften der Ergebnisse, nicht ausschließlich die textuelle Form der Standardbasis:

| Testfall | Geprüfte Eigenschaft |
|---|---|
| Glattes Modell `y^2-x` | Kein zusätzlicher Erzeuger erforderlich |
| Spitze `y^2-x^3` | Erzeuger mit `t^2=x`, `t^3=y`; Einbettung der alten Variablen |
| `y^2-x^5` | Mehrere Vergrößerungsschritte und erwartete Parametrisierung |
| `y^2-9(x^3-x)` bei `p=3` | Vertikale Nichtnormalität trotz glatter generischer Faser |
| `xy` | Nichttriviales Idempotent und Trennung der Komponenten |
| Rationale p-integrale Koeffizienten | Korrekte Behandlung von Nennern prim zu `p` |
| Verzweigte und unverzweigte Erweiterung | Relationen und p-adischer Koeffizientenwechsel |
| `T^2-5`, `T^2-8` bei p=2 | Ergänzung fehlender ganzzahliger Koeffizientenelemente |
| `T^4-14*T^2+9` bei p=2 | Gemischte Verzweigung mit e=f=2 |
| `T^2-17` bei p=2 | Zurückweisung trotz Irreduzibilität über Q |
| `3(y^2-x^3)` | Produktzerlegung, Erzeugerbilder und ursprüngliche Gleichung |
| Ungültige Eingaben und Ressourcengrenzen | Erwartete Exceptions |

Am 13. November 2026 bestanden alle 21 Tests unter SageMath 10.10 über die SageMathCell-Oberfläche. Dabei wurden die vollständigen ausgelieferten Modul- und Testdateien ausgeführt, keine verkürzte Nachimplementierung. Für den Transport wurden ihre Bytes komprimiert und in der Sage-Sitzung unverändert dekomprimiert und mit Python kompiliert; der Test-Runner lud die resultierenden Module. Die SHA-256-Prüfsummen im [`test_results.txt`](test_results.txt) stimmen mit den lokalen Dateien überein. Eine lokale Sage-Installation stand nicht zur Verfügung; zusätzlich wurden die lokalen Python-Dateien syntaktisch geprüft.

Zum Testen der ausgelieferten Dateien in einer Sage-Installation dient:

```sh
sage -python test_padic_normalization.py
```

Bei Änderungen sind insbesondere die ganzzahlige Arithmetik, die Sättigung, die Nichtnullteilereigenschaft von `d` und die exakte Behandlung von Nullrelationen vor dem p-adischen Koeffizientenwechsel zu erhalten. Ein Austausch dieser Schritte durch Rechnungen allein über `QQ` oder durch endlich genaue Nulltests würde die mathematische Aussage verändern.

## 13. Weiterführende Quellen

- Greuel, Laplagne, Seelisch: [Normalization of Rings](https://arxiv.org/abs/0904.3561), zum Grauert–Remmert-Verfahren und zur algorithmischen Normalisierung.
- Singular: [`primdecint.lib`](https://github.com/Singular/Singular/blob/spielwiese/Singular/LIB/primdecint.lib), insbesondere die ganzzahlige Radikalberechnung `radicalZ`.
- SageMath: [The Singular interface](https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/singular.html), zur verwendeten Prozessschnittstelle.
- SageMath: [Number fields](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html), insbesondere Ganzheitsbasen und Primzerlegung.
