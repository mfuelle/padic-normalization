# Implementation documentation: affine normalization over p-adic rings of integers

This document describes the structure of [`padic_normalization.py`](padic_normalization.py), its mathematical invariants, and its implementation in SageMath and Singular. The [`user guide`](USER_GUIDE.md) explains how to use it; this document focuses on how the computation works and which claims the output supports.

## 1. Objective and representation

The objective is to compute the normalization of a reduced ring

$$
A=\mathcal O_K[x,y]/(f)
$$

in its total quotient ring. For an integral domain, this is the ordinary fraction field; with multiple components, the implementation also allows a normalization that is a product of rings.

The result consists of a finite presentation

$$
\overline A\cong
\mathcal O_K[x,y,z_1,\ldots,z_m]/I
$$

and the images of the additional generators as fractions in the original variables. The original variables are retained, making the normalization embedding explicit.

**The normalization computation uses exact arithmetic throughout.** It takes place over `ZZ`. Only the optional method `padic_ring(prec=40)` transfers the result into a Sage ring with finite-precision p-adic coefficients.

The implementation is therefore not an algorithm for arbitrary p-adic coefficients known only to finite precision. It processes the exact models described below. By default it removes individually redundant normalization generators, but it does not guarantee minimum cardinality or minimize the returned relations.

## 2. Files and dependencies

| File | Purpose |
|---|---|
| `padic_normalization.py` | Public API, preprocessing, normalization loop, and Singular interface |
| `examples.sage` | Usage examples, including the ramified quadratic extension |
| `test_padic_normalization.py` | Algebraic regression tests using `unittest` |
| `README.md` | Project overview and quick start |
| `USER_GUIDE.md` | Installation, usage, and mathematical overview |
| `test_results.txt` | Report and precise scope of validation to date |
| `IMPLEMENTATION.md` | This internal documentation |

SageMath and its Singular interface, `sage.interfaces.singular.singular`, are required. The backend loads Singular's `primdecint.lib` library for `radicalZ`. The Python file does not need the Sage preparser, but it requires a Python environment with Sage; ordinary CPython without Sage is insufficient.

## 3. Public function contract

```python
normalize_affine_padic(
    f, p, extension=None,
    max_steps=20, max_variables=24, verbose=False, integral_basis=False,
    simplify=True, max_new_generators=None
)
```

The return value is a `PadicNormalization` object. Unmet preconditions raise an exception; the function does not return a merely conjectural normalization.

### 3.1 Base case: K = Q_p

The input `f` must be a nonzero, nonconstant polynomial in exactly two variables over `ZZ` or `QQ`. Variable names can be chosen freely. The following conditions apply:

- `p` is a prime number.
- Every rational coefficient is p-integral.
- The primitive part of `f` is squarefree over `QQ`.
- The minimum p-adic valuation of the nonzero coefficients is zero or one.

Valuation zero leads to the normalization loop itself. Valuation one is handled by splitting off a vertical component; see Section 8. Valuation at least two and repeated factors are rejected as nonreduced inputs. The reduction of a nonreduced ring is not performed silently.

### 3.2 Extension case

For `extension=g`, the input must be a monic polynomial `g` in `ZZ[T]` of degree at least two that is irreducible over `QQ_p`. Polynomials that are neither Eisenstein nor irreducible modulo p are also allowed. The number field F and its unique p-adic completion K satisfy

$$
F=\mathbb Q[a],\quad g(a)=0,\qquad
\mathcal O_K=\mathcal O_F\otimes_{\mathbb Z}\mathbb Z_p.
$$

In this case, `f` belongs to a polynomial ring over `ZZ` with exactly three variables. Their order has a fixed meaning: first the coefficient generator `a`, then the two curve coordinates. Rational coefficients with nontrivial denominators are not allowed in this mode.

After substituting `a` into the number field `QQ[T]/(g)`, `f` must be nonconstant, squarefree, and monic in at least one curve coordinate. Being monic is a sufficient condition for the required global flatness; it is not a necessary condition for the existence of a normalization.

A representation by an extension tower remains unsupported. Fields with mixed ramification are supported through a single suitable polynomial. Arbitrary coefficients from O_K cannot yet be inserted directly into `f`: the input remains an integer polynomial in `a,x,y`.

Eisenstein polynomials and polynomials irreducible modulo p use the existing direct path by default. For all other polynomials, the full ring of integers is presented using an integral basis. `integral_basis=True` requests this presentation for the existing cases as well.

## 4. Preprocessing and test element

### 4.1 Exact integer model

`_exact_polynomial` checks the number of variables, the coefficient ring, nonconstancy, and p-integrality. In particular, it rejects polynomials over finite-precision `Qp` or `Zp` rings.

In the base case, `_prepare` clears rational denominators and divides by the integer content. In the flat case, these factors are units in `Z_p`. The equation therefore defines the same desired ring over `Z_p`. The resulting primitive polynomial `f_0` defines a torsion-free, hence flat model

$$
A_0=\mathbb Z[x,y]/(f_0).
$$

Factorization over `QQ` checks that the generic fibre is reduced. Together with torsion-freeness, this implies that `A_0` is reduced.

In the direct extension case, the integer model is

$$
A_0=\mathbb Z[a,x,y]/(g(a),f_0).
$$

Since `g` is monic and `f_0` is monic in a curve coordinate, the model is torsion-free over `ZZ`. Powers of `a` are reduced exactly modulo `g(a)`.

### 4.2 Local irreducibility and integral basis

`_prepare` first checks that `g` is irreducible over Q. If neither the Eisenstein criterion nor irreducibility modulo p applies, it then computes `F.primes_above(p)`. There must be exactly one prime P, with `e(P/p)*f(P/p)=deg(g)`. This is an exact check that `F tensor_Q Q_p` is a field. A polynomial such as `T^2-17`, which is irreducible over Q but splits over Q_2, is rejected.

`_integral_coefficient_model` uses `F.integral_basis()` to compute a full Z-basis `b_0,...,b_(n-1)` of O_F. It does not assume that the first basis vector is one. Rational linear algebra gives integer coordinates for one, a, and every product of two basis vectors:

$$
1=\sum_k u_kb_k,\qquad a=\sum_k v_kb_k,
\qquad b_i b_j=\sum_k c_{ijk}b_k.
$$

These define the coefficient ring

$$
C=\mathbb Z[a,B_0,\ldots,B_{n-1}]/
\left(1-\sum_k u_kB_k,\ a-\sum_k v_kB_k,
\ B_iB_j-\sum_k c_{ijk}B_k\right).
$$

This ring is exactly O_F: the relations reduce every polynomial to a Z-linear combination of the basis variables. Evaluation in O_F is surjective, and linear independence of the integral basis proves injectivity.

The original variables `a,x,y` are retained. The additional variables are named `ib0`, `ib1`, …; existing names are skipped. The multiplication relations and the original curve equation are combined in a single integer polynomial ring. The normalization loop thus starts from `O_F[x,y]/(f)`. The monicity of f also guarantees flatness and the embedding of the original order model here.

For the fraction representations, all basis vectors are stored as rational polynomials in a. Their integer denominators are nonzerodivisors in the flat starting model. Adding the integral basis is preprocessing; it does not appear as a Grauert–Remmert step in `steps`, but it counts towards `max_variables`.

### 4.3 Choosing d

The loop requires a nonzerodivisor `d` such that `A_0[1/d]` is normal. The implementation chooses it from a Jacobian ideal so that this localization is even smooth over `ZZ`.

In the base case, the nonzero derivatives of `f_0` are sorted by total degree. Linear combinations are tested for small integers `c`. With two derivatives, they have the form

$$
d=h_0+c h_1.
$$

A combination is accepted if `gcd(f_0,d)=1` over `QQ`. Squarefreeness in characteristic zero ensures that no irreducible component divides both partial derivatives. The finite search used in the code therefore finds a suitable combination. Torsion-freeness over `ZZ` transfers the nonzerodivisor property to the integer model.

In the extension case, the corresponding gcd check is performed over `QQ[T]/(g)`. For a suitable combination `h` of derivatives with respect to the curve coordinates, the algorithm uses

$$
d=g'(a)h.
$$

This is a linear combination of the relevant Jacobian minors of the two equations `g(a), f_0`. On the direct path, the factor `g'(a)` also accounts for possible nonnormality of the global coefficient model at other primes.

In integral basis mode, the coefficient ring is already the normal ring O_F. The factor `g'(a)` is therefore removed: inverting h makes the curve smooth over O_F, which suffices for normality. The implementation also computes the relative Jacobian ideal of the curve over O_F. If its standard basis contains a nonzero integer, that integer is used as d instead. Its principal open is likewise smooth over O_F. A constant test element can greatly reduce the cost of `radicalZ`, since the test ideal then has no characteristic-zero component. If no constant is found, h is retained.

Throughout the iteration, `d` remains the image of this original element. It is simply transferred into each larger presentation.

## 5. Normalization loop

### 5.1 Mathematical step

Let the current ring `A_i` be reduced and let `d` be the chosen nonzerodivisor. Set

$$
J_i=\sqrt{dA_i}.
$$

Then `J_i` contains a nonzerodivisor, and its vanishing locus contains the nonnormal locus. The Grauert–Remmert criterion gives

$$
A_i\text{ normal}
\quad\Longleftrightarrow\quad
\operatorname{End}_{A_i}(J_i)=A_i.
$$

As a subring of the total quotient ring, the endomorphism ring can be written as

$$
A_{i+1}=\operatorname{End}_{A_i}(J_i)
=\frac{1}{d}(dJ_i:_{A_i}J_i).
$$

This extension is finite, integral, and birational. If it is strictly larger, generators `u_j/d` are adjoined and all relations are computed. If it is not larger, the loop terminates with a proof of normality.

When `max_new_generators` is a positive integer, only that many fractions are
selected in one iteration. The next ring can then be a proper intermediate
algebra between A_i and its endomorphism ring. It is still a strict finite
integral birational enlargement. The full normality criterion is checked on
each subsequent iteration; a batch limit never permits early certification.

### 5.2 Translation into ideals of a polynomial ring

The code does not work directly with ideals of a quotient ring. It maintains a presentation `A_i=R/I` with a polynomial ring `R` over `ZZ` and works with preimages in `R`:

| Code variable | Meaning in the ambient polynomial ring |
|---|---|
| `I` | Defining ideal of `A_i` |
| `J` | `sqrt(I + (d))`, the preimage of `J_i` |
| `H` | `I + d*J`, the preimage of `d*J_i` |
| `C` | `(H:J)`, the preimage of `(d*J_i:J_i)` |
| `old` | Initially `I + (d)`, then enlarged by the numerators already selected |

The generators of `C` are reduced modulo `old` one at a time. A nonzero remainder is stored as a new numerator and added to `old`. This avoids generators already accounted for by earlier numerators; it does not guarantee a minimal generating set.

If the numerator list remains empty, then

$$
C=I+(d),
$$

which is the termination criterion. Otherwise, the list provides generators of the next finite algebra. `max_steps` is checked only after determining that another enlargement is required.

### 5.3 Invariants and termination

After each successful step:

1. The current algebra is a finite subalgebra of the normalization of the starting model.
2. All new elements are represented as fractions with nonzerodivisor denominators in the original ring.
3. After inverting `d`, the algebra is unchanged and normal.
4. The presentation contains the full relation kernel of the adjoined fractions.

Finitely generated `ZZ`-algebras are excellent; in the reduced case, their normalization is finite. The ascending chain therefore lies in a finite module over the Noetherian starting ring and stabilizes. This is the theoretical termination statement. It does not provide a useful general runtime bound for the Gröbner basis computations.

### 5.4 Removing redundant generators

After each enlargement, `simplify=True` calls `_simplify_presentation`. The
variables of `integral_original_ideal` are protected: these include the
original curve variables and, where applicable, the coefficient basis
variables. Other generators are considered in order, oldest first.

For a candidate z, the ideal is transported into an integer polynomial ring
with a leading one-variable lexicographic block for z. The remaining variables
form a degree reverse lexicographic block. The routine computes a standard
basis G and the normal form P of z. If P still involves z, the generator is
retained. Otherwise it verifies `z-P` reduces to zero modulo G and substitutes
P into all relations, constructing a new standard basis in the smaller ring.

This detects whether the image of z belongs to the subalgebra generated by
the remaining variables. If `z-P` belongs to I, substitution induces an
isomorphism

$$
\mathbb Z[\mathbf w,z]/I
\cong \mathbb Z[\mathbf w]/I\big|_{z=P(\mathbf w)}.
$$

No new fraction or coefficient inversion is introduced. In particular, a
relation `2*z-x` is not sufficient to eliminate z over ZZ. The current test
element d is transported by the same substitution, and the fraction images
of retained variables are kept. Thus normality, the birational embedding,
and the loop invariants are preserved.

A failed redundancy check need not be repeated during the same pass: removing
other generators only shrinks the candidate subalgebra. Once the pass finishes,
no single unprotected variable can be removed over this exact global model.
This is an irredundant subset of the existing generators, not a globally
minimal generating set. Changing generators jointly, or allowing additional
coefficient units after p-adic base change, may permit further reductions.

`max_variables` is checked before constructing the enlarged graph, so it also
bounds the temporary presentation before this simplification. It does not
count auxiliary elimination variables. `max_new_generators=1` limits each
enlargement to one new variable but does not bound the total number required.

## 6. Singular backend: `_Engine`

### 6.1 Standard bases and radicals

`std` uses the embedded Singular helper `pn_std`. It computes a standard basis over `ZZ`, normalizes the signs of the leading coefficients, and marks the result as a standard basis.

`radical` calls `radicalZ` from `primdecint.lib`. It computes the radical in the integer polynomial ring itself. A radical computed only over `QQ` would miss vertical nonnormality.

`nf` computes normal forms using Singular's `reduce`; `subset` checks ideal inclusions using normal forms of the generators.

### 6.2 Ideal quotient by integer elimination

`_Engine.colon(H,J)` computes `(H:J)` over `ZZ`. Write

$$
H=(h_1,\ldots,h_m),\qquad J=(j_1,\ldots,j_r).
$$

First, Python reduces the generators of J modulo I and removes zero remainders. Since I is contained in H, these steps do not change the ideal quotient. H is replaced by a standard basis.

Starting with K=(1), the method processes each nonzero generator j of J by computing

$$
L=H\cap jK,\qquad K_{\mathrm{new}}=L/j=K\cap(H:j).
$$

The intersection is obtained by eliminating an auxiliary variable u from `u*H + (1-u)*j*K`. Each resulting generator is divided exactly by j, and the identity `j*quotient == generator` is checked. This division is valid because the ambient ring `ZZ[variables]` is a domain. After each iteration, K is replaced by a standard basis. At the end, K=(H:J); for the zero ideal J, K remains (1).

All ideal computations remain over ZZ. Extending H to QQ would lose integer torsion and is not allowed. This method replaces the integer syzygy routine used before November 14, 2026, which exhausted memory on the reported degree-eight extension example.

Python additionally checks `K*J ⊆ H`. This detects certain backend or conversion errors; completeness relies on the elimination computation and exact division.

### 6.3 Full relations by saturation

When adjoining the fractions `u_1/d,...,u_s/d`, the equations

$$
dz_j-u_j=0
$$

alone are generally insufficient: they can introduce additional components over `d=0`. The correct relation kernel is

$$
I_{\mathrm{new}}=(I,dz_1-u_1,\ldots,dz_s-u_s):d^\infty.
$$

`graph` computes this kernel by introducing an auxiliary variable `t` and eliminating it:

$$
I_{\mathrm{new}}=
(I,dz_1-u_1,\ldots,dz_s-u_s,1-td)
\cap R[z_1,\ldots,z_s].
$$

The block term order places `t` in a leading lexicographic block. The polynomials independent of `t` are selected from the standard basis. After removing the auxiliary variable, another standard basis is computed.

Before constructing this elimination ideal, common polynomial factors in each `u_j/d` are cancelled exactly. The graph uses these reduced numerators and denominators and inverts the least common multiple of the reduced denominators. This gives the same relation kernel: every divisor of the original nonzerodivisor d is a nonzerodivisor as well. Avoiding removable factors substantially reduces memory use in some examples.

The additional Python check ensures that the old relations and all equations `dz_j-u_j` belong to the new ideal. Completeness of the computed kernel follows from the elimination computation; the inclusion check alone would not prove it.

### 6.4 Interface details

- Before each Singular call, the appropriate ring is set as the current ring.
- Results are explicitly converted back into the intended Sage ring using `sage_poly(R)`.
- The helpers are loaded when an `_Engine` is created and reloaded if Singular was restarted between operations. A restart during an operation raises `NormalizationBackendError`, naming the interrupted operation. The result is discarded; the computation is not silently retried or certified.
- The helper text is passed to `singular.eval` as a single line. This prevents the Expect interface from waiting on an incomplete procedure declaration. Consequently, `//` line comments must not be added to the helper text.
- The elimination variable is initially named `pn_inverse`; underscores are appended in case of name collisions. New algebra generators are named `z1`, `z2`, …; names already in use are skipped.
- The shared Singular session has mutable state. Parallel calls within the same session are not protected; separate processes should be used for parallel computations.

## 7. Result object and maps

### 7.1 Attributes

| Attribute | Contents |
|---|---|
| `p` | Prime number as a Sage integer |
| `original_ideal` | Ideal of the prepared integer starting model |
| `integral_original_ideal` | Starting ideal after adding the integral basis, if applicable |
| `coefficient_model` | Full coefficient ring data in integral basis mode; otherwise `None` |
| `ambient_ring` | Integer polynomial ring of the final presentation |
| `ideal` | Final relation ideal |
| `relations` | Tuple of generators of `ideal` |
| `fraction_images` | Dictionary: variable name → fraction in the original variables |
| `steps` | Tuple containing the enlargement step data |
| `simplifications` | Sequential records of removed variable names and their polynomial replacement expressions |
| `extension` | Exact polynomial `g` or `None` |
| `certificate` | Text describing the mathematical termination condition reached |

`original_ideal` may differ from the equation entered literally, because denominators and content factors invertible at `p` are removed. In the extension case, it still denotes the model over the original order `ZZ[a]`. In integral basis mode, `integral_original_ideal` instead contains the actual starting model over O_F.

`coefficient_model` contains `ambient_ring`, `ideal`, and `names` for the coefficient presentation, together with `number_field`, `basis`, `ramification_index`, and `residue_degree`. All basis data and relation coefficients are exact.

An ordinary entry in `steps` has the form

```python
{
    'variables': ('z1', ...),
    'numerators': (u1, ...),
    'denominator': d,
}
```

The polynomials in this entry belong to the polynomial ring **before** the corresponding step. With a vertical component, a final entry `{'vertical_idempotent': name}` is added. In that case, `len(steps)` is therefore not simply the number of Grauert–Remmert enlargements.

Simplification does not rewrite `steps`: it is a history and can mention
removed generators. Each `simplifications` entry has keys `variable` and
`expression`; the expression belongs to the smaller ring immediately after
that removal, not necessarily the final ring. Later removal substitutions
must be applied in sequence if a final-ring expression is needed. For a
vertical product, the earlier records describe the horizontal normalization.

### 7.2 Fraction representations

For bookkeeping, the code uses the rational function field of the original ambient polynomial ring. New fractions are obtained by substituting the generator images already known. This auxiliary field must not be confused with the total quotient ring of the starting ring.

The expressions in `fraction_images` are fraction representatives **modulo `original_ideal`**. For example, when `y^2=x^3`, the expressions `x^2/y` and `y/x` represent the same element, although they differ in the free rational function field. An ordinary comparison of the two Python fraction objects does not test this equality.

To check such an equality, clear the nonzerodivisor denominators and reduce the numerator of the difference modulo `original_ideal`. The code constructs the fractions with valid denominators.

### 7.3 Methods

`integer_ring()` constructs the exact Sage quotient ring `ambient_ring/ideal`.

`simplified(verbose=False)` applies the same redundancy removal to an existing
result and returns a new `PadicNormalization` object. It retains the original
source ideals, coefficient model, historical steps, and certificate, and
extends the simplification history. The original object is unchanged. This
supports postprocessing a calculation performed with `simplify=False`.

`integer_map()` returns `(A0, B0, phi)`: the prepared integer order model, its normalization, and the embedding. If the original order is nonmaximal, A0 is not the full O_F-model. Since the original variables remain at the beginning of the variable list, their residue classes map to the residue classes with the same names in `B0`.

In integral basis mode, `coefficient_ring()` returns the exact quotient ring C from Section 4.2. With `prec=N`, the coefficients of the multiplication table are transferred to `Zp(p,N)`. This is a presentation of O_K, not a native p-adic DVR class with a valuation or uniformizer API. Outside integral basis mode, the method requires a new normalization call with `integral_basis=True`.

`ring_over_integers()` presents the normalization as an affine algebra over the exact coefficient ring C. For this purpose, `_relative_ring` moves the coefficient variable a and the integral basis variables into the base ring. First, the relations in the integer polynomial ring are reduced exactly modulo the coefficient ideal; exact zero relations are removed.

`normalization_map()` returns the exact models `(A,B,phi)` in flattened polynomial rings over ZZ, with the source constructed from `integral_original_ideal`. In integral basis mode, both rings include the coefficient variables and thus define O_F-algebras; their Sage base rings are ZZ. This is the embedding with the full ring of integers in the source. Its base change to Z_p gives the desired normalization embedding over O_K.

`coefficient_map(prec=None)` returns `(C,B,iota)`, consisting of the full coefficient ring, the flattened presentation in a single polynomial ring, and the structure map. Without a precision argument, B is the exact `integer_ring()`; with `prec=N`, C and B are presented over `Zp(p,N)`. The images of the coefficient generators are the corresponding variables of B.

Both map methods check the required relation inclusions in the exact flattened integer model. `normalization_map()` then creates an ordinary Sage ring homomorphism between the flattened integer models. `coefficient_map()` instead uses the small callable class `CoefficientPresentationMap`, with `domain()`, `codomain()`, and `im_gens()`. It evaluates polynomial representatives and constructs the result without quotient reduction. It is not an element of a native Sage Homset. This avoids normal form computations that Sage may require even when creating homomorphisms with relation checking disabled. Further computations with nested or p-adic quotients may still encounter Sage limitations; the common integer model remains the reference for these computations.

On the existing direct extension path, `padic_ring(prec=40)` uses `Zp(p,prec).extension(g)` and replaces a by the generator of the coefficient ring. In integral basis mode, like `ring_over_integers(prec=40)`, it instead returns a **flattened presentation in a single polynomial ring**, `Zp[a,x,y,ib0,...,z1,...]/I`, containing all coefficient relations. Its Sage `base_ring()` is Z_p. The additional O_K-algebra structure is explicitly described by `coefficient_model` and `coefficient_ring(prec=40)`. This avoids unsupported Sage arithmetic over nested p-adic quotient rings.

On the direct extension path, the method reduces the relations exactly modulo `g(a)` before the finite-precision base change and removes exactly vanishing polynomials. Thus, the exact zero relation `g(a)=0` is not added to the relation ideal as a mere expression `O(p^N)`. In the flattened integral basis presentation, a and the basis variables remain polynomial variables, so their equations remain as relations. No exact zero test is inferred from finite-precision p-adic data.

The returned p-adic object is a Sage quotient ring. Further arithmetic operations, particularly Gröbner basis operations over this coefficient ring, have limited availability. The exact specification and the proof of normality remain with the integer result object.

`certificate` describes the condition reached; it is not a separately serialized or independently formally verified proof object. Correctness assumes the validity of the exact Sage/Singular operations used.

## 8. Vertical component in the base case

When the p-adic content has valuation one, the input has the form `p*f_0` up to a unit, with `f_0` primitive. Then

$$
\overline{\mathbb Z_p[x,y]/(pf_0)}
\cong
\overline{\mathbb Z_p[x,y]/(f_0)}
\times\mathbb F_p[x,y].
$$

The function first normalizes the horizontal factor recursively. `_with_vertical_component` then presents the product as a single affine quotient ring.

Let `H_l=0` be its horizontal relations and let `z_j` be the additional horizontal generators. The product presentation uses a new generator `e` with

$$
e^2=e,\qquad p(1-e)=0,\qquad eH_l=0,\qquad (1-e)z_j=0.
$$

The original variables `x,y` are retained on both factors. On the horizontal factor, `e=1`; on the vertical factor, `e=0`. The name is initially `vertical_split` and is adjusted in case of collisions.

The embedding is also described explicitly:

$$
e=\frac{p}{p+f_0}.
$$

A horizontal generator represented by the fraction `u/v` is replaced by

$$
\frac{pu}{pv+f_0}.
$$

This fraction equals `u/v` horizontally and zero vertically. Its denominator is a nonzerodivisor on every component: it restricts to `p*v` horizontally and to the nonzero polynomial `f_0 mod p` vertically.

This additional handling applies only to the base case `K=Q_p`. Extension mode requires the flat model described above.

## 9. Why the integer computation solves the p-adic problem

The implementation normalizes an integer model globally. It may therefore introduce generators needed only at primes other than `p`. They become redundant after base change but do not invalidate the result.

The base change `ZZ → Z_p` is flat and regular. Normality is preserved under this regular base change. The finite birational normalization embedding remains finite and is still an isomorphism after inverting the nonzerodivisor `d`. The resulting normal algebra is therefore the normalization of the p-adic starting model.

On the direct extension path, `Z_p[a]/(g(a))` is already O_K by the Eisenstein condition or irreducibility modulo p. The general path uses the full ring of integers O_F instead. The exactly verified uniqueness of the prime above p gives `O_F tensor_Z Z_p = O_K`. The same base change argument therefore yields the desired O_K-algebra even when the original order is nonmaximal or the ramification is mixed.

Here, the change concerns the **coefficients** of a polynomial ring. Neither a power series ring nor a Tate algebra is constructed. The full p-adic completion of the entire affine ring is not the objective of the function either.

## 10. Example: Q_3(sqrt(3))

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

The documented computation found a generator `z1=y/a`. The returned integer presentation had the relations

$$
az_1-y,\quad y^2-3z_1^2,\quad ay-3z_1,
\quad a^2-3,\quad x^3-z_1^2-x.
$$

Interpreting `a` as a coefficient and eliminating `y` gives the shorter presentation

$$
\overline A\cong\mathcal O_K[x,z]/(z^2-x^3+x),
\qquad y\longmapsto az.
$$

This simplification explains the output; automatic minimization is not part of the module. Generator names, signs, and the order of a standard basis should not be treated as canonical.

## 11. Error handling and computational cost

| Exception | Meaning |
|---|---|
| `TypeError` | Unsuitable coefficient ring or unsupported extension description |
| `ValueError` | An input precondition is violated, such as reducedness, primality, or monicity |
| `NormalizationLimitError` | Another enlargement would exceed a configured resource limit |
| `NormalizationBackendError` | Singular restarted during an operation; no certified result is returned |
| `ArithmeticError` | The test element could not be constructed, or an internal consistency check failed |

Other exceptions from the Sage/Singular backend are not caught in general.

`max_steps` bounds the number of actual enlargement steps. Even with `max_steps=0`, the algorithm checks whether the starting model is already normal. `max_variables` is checked before adding the integral basis, before adjoining new generators, and when splitting off the vertical component. It is not a general memory limit and, in particular, does not count the temporary elimination variable.

There is no timeout within an individual radical, intersection, or standard basis computation. Fixed time or memory limits require execution in a separate process. If Linux kills Singular for exceeding available memory, Sage may restart it and report a missing helper function. The engine detects the session change and reports the backend interruption explicitly. Raising `max_steps` or `max_variables` does not provide more memory.

The dominant costs come from `radicalZ`, ideal intersection elimination, and the elimination basis of the graph presentation. Numerous additional variables, high degrees, and global contributions from other primes may greatly increase the computation. The configured limits are safeguards, not a promise that an input below them will be processed quickly.

## 12. Validation and maintenance

The tests check algebraic properties of the results, rather than only the textual form of the standard basis:

| Test case | Property checked |
|---|---|
| Smooth model `y^2-x` | No additional generator required |
| Cusp `y^2-x^3` | Generator with `t^2=x`, `t^3=y`; embedding of the original variables |
| `y^2-x^5` | Multiple enlargement steps and the expected parametrization |
| `y^2-9(x^3-x)` at `p=3` | Vertical nonnormality despite a smooth generic fibre |
| `xy` | Nontrivial idempotent and separation of components |
| Rational p-integral coefficients | Correct handling of denominators coprime to `p` |
| Ramified and unramified extensions | Relations and p-adic base change |
| `T^2-5`, `T^2-8` at p=2 | Addition of missing integral coefficient elements |
| `T^4-14*T^2+9` at p=2 | Mixed ramification with e=f=2 |
| `T^2-17` at p=2 | Rejection despite irreducibility over Q |
| `3(y^2-x^3)` | Product decomposition, generator images, and original equation |
| Invalid inputs and resource limits | Expected exceptions |

On November 13, 2026, all 21 tests passed under SageMath 10.10 through the SageMathCell interface. The complete module and test files of that version were executed, rather than an abbreviated reimplementation. For transport, their bytes were compressed, then decompressed unchanged and compiled with Python in the Sage session; the test runner loaded the resulting modules. No local Sage installation was available at that time; the local Python files also passed a syntax check.

On November 14, the backend update passed 24 tests under SageMath 10.10 in Arch Linux/WSL. The generator simplification update subsequently passed all 27 tests. Additional coverage includes nonlinear generator elimination with transport of the test element, refusal to divide by an integer nonunit, simplification during and after multi-step normalization, preservation of the original embedding, and bounded batches. [`test_results.txt`](test_results.txt) records hashes for the corresponding versions.

The degree-eight extension example reported by the user was reproduced during the backend update. Kernel logs confirmed out-of-memory termination of Singular. With that update, its first enlargement completed; the next radical computation was stopped to limit test time. This is a partial performance check, not a completed normalization or a full regression-test pass for that input. The large example was not rerun after adding generator simplification.

To test the distributed files in a Sage installation, run:

```sh
sage -python test_padic_normalization.py
```

When making changes, preserve in particular the integer arithmetic, saturation, the nonzerodivisor property of `d`, and the exact handling of zero relations before p-adic base change. Replacing these steps with computations only over `QQ` or with finite-precision zero tests would change the mathematical claim.

## 13. Further references

- Greuel, Laplagne, Seelisch: [Normalization of Rings](https://arxiv.org/abs/0904.3561), for the Grauert–Remmert method and algorithmic normalization.
- Singular: [`primdecint.lib`](https://github.com/Singular/Singular/blob/spielwiese/Singular/LIB/primdecint.lib), particularly the integer radical computation `radicalZ`.
- SageMath: [The Singular interface](https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/singular.html), for the process interface used here.
- SageMath: [Number fields](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html), particularly integral bases and prime decomposition.
