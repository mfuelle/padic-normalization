# Affine normalization over p-adic rings of integers

`padic_normalization.py` is a SageMath module. It computes a **finite affine presentation of the normalization, including the embedding**. The algorithm uses exact arithmetic over the integers, rather than rounded p-adic coefficients.

The separate [implementation documentation](IMPLEMENTATION.md) explains the internal data structures, the Grauert–Remmert step, the Singular interface, the computation of the full relation ideal, and the scope of validation.

## Installation and usage

Place the files in the same directory. In a Sage notebook with that directory as its working directory:

```sage
from padic_normalization import normalize_affine_padic

R.<x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2 - 9*(x^3-x), p=3, verbose=True)
print(N)
N.relations               # Exact equations over the integers
N.fraction_images         # Images of all generators in the total quotient ring
B = N.padic_ring(prec=40) # Sage quotient ring over Z_3
```

Alternatively, use `load("/full/path/padic_normalization.py")` in Sage. The module is ordinary Python with Sage imports; it does not require the Sage preparser. It does not run in an ordinary Python installation without Sage.

## Supported inputs

### K = Q_p

* `f` belongs to a polynomial ring in exactly two variables over `ZZ` or `QQ`.
* All coefficients are p-integral. Rational denominators coprime to p are allowed.
* The p-adic content has valuation at most one; when its valuation is one, the vertical component is split off explicitly.
* The primitive part of `f` is squarefree over Q. Reducible but reduced rings are allowed; the normalization may be a product and is still presented as a single affine quotient ring.

This version therefore handles **exactly specified reduced models with a nonzero, nonconstant defining equation**, rather than arbitrary finite-precision elements of `Qp(p)` or `Zp(p)`. Such inputs are explicitly rejected. Nonreduced rings are rejected. A reduced vertical component is retained and returned as a product factor.

### Extensions K/Q_p

The module optionally supports `extension=g`, where `g` is a monic polynomial in `ZZ[T]` of degree at least two that is **irreducible over Q_p**. It need not be Eisenstein or irreducible modulo p. In particular, extensions with mixed ramification are covered when described by a single such polynomial.

Local irreducibility is checked exactly: after checking irreducibility over Q, the number field `F=QQ[T]/(g)` must have exactly one prime above p. Eisenstein polynomials and polynomials irreducible modulo p retain their existing direct computation path.

In this case, `f` is a polynomial in **three** variables over `ZZ`, in the order `(a,x,y)`. The first variable is a coefficient generator satisfying `g(a)=0`. In addition, after reduction modulo `g`, `f` must be monic in x or y. This sufficient condition guarantees a flat global starting model. The generic fibre must again be reduced.

```sage
V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2 - 3*(x^3-x), 3, extension=T^2-3)
print(N)
print(N.fraction_images)
B = N.padic_ring(40)  # Presentation over O_{Q_3(sqrt(3))}
```

In the general case, `Z_p[a]` may be a proper suborder of `O_K`. The module therefore computes a full integral basis of F and starts the normalization from `O_F[x,y]/(f)`. After base change to Z_p, this is the desired model over O_K.

```sage
V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2-x^3, 2, extension=T^2-5)

O = N.coefficient_ring()        # Exact presentation of the full O_F
B0 = N.ring_over_integers()     # Normalization as an affine O_F-algebra
A0, B0, phi = N.normalization_map()
B = N.padic_ring(prec=40)       # Flattened Z_2 presentation including the O_K coefficient relations
print(N.coefficient_model['basis'])
```

Here `(1+a)/2` belongs to the ring of integers but not to `Z_2[a]`. For the example with mixed ramification, `g=T^4-14*T^2+9`, use the same call with that extension polynomial; at p=2 it has ramification index 2 and residue degree 2.

The option `integral_basis=True` explicitly requests the integral basis presentation for the existing Eisenstein and unramified cases as well. Without this option, their existing output format is retained. The specific generators and standard bases are not canonical and may change.

The input syntax for `f` remains restricted to integer coefficients in `(a,x,y)`. Arbitrary coefficients from O_K, such as `(1+a)/2` directly in the input equation, finite-precision p-adic coefficients, and extension towers are not supported. Irreducibility of `g` over Q alone is insufficient: for example, `T^2-17` is rejected at p=2 because it splits over Q_2.

## Interpreting the output

* `N.ambient_ring`, `N.ideal`, `N.relations`: an exact presentation over ZZ. For Q_p, interpret the same equations over Z_p. For an extension, the equation `g(a)=0` belongs to the integer model; over O_K, replace a by the chosen coefficient generator.
* `N.integer_ring()`: the exact Sage quotient ring of the global model.
* `A, B, phi = N.integer_map()`: exact global models and the embedding. The source still uses the original order `ZZ[a]`; when this order is nonmaximal, it is not the full O_F starting model. The original variables are retained.
* `N.coefficient_model`: in integral basis mode, data for the full coefficient order, including `basis`, `ramification_index`, and `residue_degree`; otherwise `None`.
* `N.coefficient_ring()` and `N.ring_over_integers()`: exact presentations of O_F and the normalized affine O_F-algebra. `coefficient_ring(prec=40)` returns the coefficient presentation over finite-precision Z_p. `ring_over_integers(prec=40)` returns the flattened p-adic presentation in a single polynomial ring described below.
* `A, B, phi = N.normalization_map()`: the exact flattened models over ZZ and the desired embedding. In integral basis mode, both models include the full ring of integers through coefficient variables; the source is based on `N.integral_original_ideal`.
* `C, B, iota = N.coefficient_map(prec=40)`: the coefficient ring, the flattened presentation, and the explicit structure map from C to B. Without `prec`, all three are exact over ZZ. The defining relations are checked exactly in advance. `iota` is a callable `CoefficientPresentationMap` operating on polynomial representatives, rather than a native Sage Homset element.
* `N.fraction_images`: each new generator represented as a fraction in the original variables. These fractions are interpreted **modulo the original relations**; they need not be equal in the ambient rational function field. Their denominators are nonzerodivisors in the starting ring.
* `N.steps`: the successively adjoined fractions, each expressed in the variables of the preceding step.
  When a vertical component is present, the last entry instead contains the name of the separating idempotent.
* `N.simplifications`: successive removal records with `variable` and `expression`. Each expression belongs to the smaller polynomial ring immediately after that removal. Later substitutions may be needed to express it in the final variables. `steps` remains a history of enlargements, so it may mention variables subsequently removed.
* `N.certificate`: termination by the Grauert–Remmert criterion, not merely by smoothness of the generic fibre.
* `N.padic_ring(prec=40)`: on the existing direct path, a quotient ring over Sage's p-adic ring of integers. In integral basis mode, a **flattened presentation in a single polynomial ring over Z_p**, with a and the integral basis elements as additional variables and with their coefficient relations. Its `base_ring()` is Z_p; the full O_K-algebra structure is described by `coefficient_model` or `coefficient_ring(prec=40)`. Sage does not sufficiently support arithmetic over a nested p-adic quotient coefficient ring. The result is therefore not presented as a native O_K DVR class. **The exact specification remains `N.relations`**; the p-adic object does not replace an exact equality test.

By default, the program removes individually redundant normalization generators after every enlargement. It does not search for a globally smallest set of algebra generators or minimize the list of defining relations. A shorter isomorphic presentation may still exist.

Sage does not implement every quotient operation over the exact nested coefficient ring either. For ideal computations and equality tests, `integer_ring()` or `ideal` in the common integer polynomial ring is the reference.

## Controlling the number of generators

Use `simplify=True` (the default) to remove redundant generators after each
enlargement. A generator z is removed only after finding and verifying an exact
relation `z=P(other variables)`. The algorithm substitutes P into all defining
relations and into its current test element. The original variables and all
coefficient presentation variables are protected, so the original embedding
and coefficient structure maps remain available.

For example, the normalization of `y^2-x^5` can acquire intermediate elements
corresponding to `t^3=y/x` and then `t=y/x^2`. Once t is present, the older
element is redundant because `t^3=x*t`. The simplified output retains only
one extra parameter, with `x=t^2` and `y=t^5`.

```sage
R.<x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2-x^5, 5, simplify=True,
                           max_new_generators=1, max_steps=100,
                           verbose=True)
print(N.ambient_ring.variable_names())
print(N.simplifications)

# Alternatively, simplify only after normalization has finished.
raw = normalize_affine_padic(y^2-x^5, 5, simplify=False)
small = raw.simplified(verbose=True)  # raw is unchanged
```

`max_new_generators` is a positive integer or `None` (the default). With a
finite value, only that many of the computed fractions are adjoined per
enlargement. Each resulting ring is still a finite integral birational
extension, and the normality criterion is checked again on the next iteration.
This can require more iterations, so `max_steps` may need to be increased.

`max_variables` remains a hard guard on the total presentation size **before**
removing redundancies. It includes original variables, coefficient variables,
and normalization generators, but not temporary elimination variables. It
does not promise that normalization is possible within the chosen bound.
Leave room for at least one new variable even if an older variable is expected
to disappear immediately afterwards.

The simplifier is conservative: it works over the exact global integer model,
not just its p-adic base change. It never divides by a nonunit integer to remove
a variable. Thus it can miss a relation available only after localizing away
from p. Among the current unprotected generators, no single generator can be
deleted over that global model when the check finishes. This is **not minimum
cardinality**: replacing several generators by new combinations may allow a
smaller set, and the order of deletions may affect the result. The coefficient
presentation itself is not minimized.

Each redundancy check needs an elimination standard basis and may be costly.
Smaller presentations can make later normalization steps cheaper, but this
extra work is not guaranteed to reduce total running time.

## Algorithm and correctness

First, only p-adic units are removed from rational denominators and the content of f. This produces an exact integer model A0. For the existing direct extension cases, A0 is defined by `(g(a), f(a,x,y))`. Integral basis mode additionally introduces basis variables and the multiplication table of O_F. The loop therefore starts over the full ring of integers; the original order data are retained for the fraction representations.

For `f=p*f0` with primitive f0, the algorithm first normalizes f0. The full normalization is then its product with `F_p[x,y]`. An additional idempotent e presents this product as an affine ring; in the original variables, `e=p/(p+f0)`. The iterative algorithm below operates on the flat horizontal model.

A nonzerodivisor d is chosen from the relative Jacobian ideal. On the direct extension path, it is additionally multiplied by g'(a), so that inverting d makes A0 smooth over ZZ and hence normal. In integral basis mode, O_F is already normal and this extra factor is unnecessary: smoothness over O_F suffices. If a nonzero integer is found in the relative Jacobian ideal, it is preferred as the test element; otherwise the curve derivative combination is used. Squarefreeness and a gcd test verify the nonzerodivisor property of the derivative combination.

For each current ring Ai, the algorithm computes:

```
Ji = radical(d*Ai)
Ci = (d*Ji : Ji)
Ai+1 = (1/d)*Ci = End_Ai(Ji).
```

With `max_new_generators` set, Ai+1 may instead be a proper intermediate
algebra obtained by adjoining only the selected fractions from this
endomorphism ring. Termination still requires the full normality criterion.
Removing redundant variables changes only the presentation of Ai+1.

`radicalZ` from Singular's `primdecint.lib` computes the radical over ZZ. The ideal quotient is computed by integer ideal intersections and exact division: starting with K=(1), each nonzero generator j of Ji replaces K by `(H intersect j*K)/j`, where H is the preimage of d*Ji. Intersections are computed by elimination over ZZ. The ordinary `quotient` is not used; the Singular library documents limitations of that function over ZZ.

The new generators are the computer-generated elements `u/d`. Their full relations are determined by eliminating an auxiliary variable from

```
(old relations, d*z_i-u_i, 1-s*d)
```

This is saturation with respect to d. Using only the equations `d*z_i-u_i` would be incorrect because additional components on d=0 may appear.

The implementation first cancels common polynomial factors in each fraction. It then uses the reduced graph equations and saturates by the least common multiple of their denominators. This preserves the full relation kernel while avoiding unnecessary elimination work.

When `Ci = d*Ai`, the Grauert–Remmert criterion establishes normality: Ji is radical, contains the nonzerodivisor d, and its vanishing locus contains the nonnormal locus. The last property is preserved because each extension step becomes an isomorphism after inverting d. Otherwise, Ai strictly increases within the finite normalization. Excellence and Noetherianity therefore guarantee theoretical termination.

The normalization is first computed over ZZ. Localization at p and the regular base change to Z_p preserve normality and yield the normalization of the desired model. In the extension case, `O_F tensor_Z Z_p = O_K` because there is exactly one prime above p. This is a **change of coefficients**, not the adic completion of the entire polynomial ring.

## Performance and limitations

Run the supplied test suite with `sage -python test_padic_normalization.py`. The file `examples.sage` contains directly executable examples (`sage examples.sage`).

The current version passes all 27 regression tests under SageMath 10.10 in Arch Linux/WSL (November 14, 2026). The tests cover nonmaximal orders, mixed ramification, local irreducibility, output formats, integer ideal quotients, Singular restarts, generator removal, and bounded batches. `test_results.txt` records the current source hashes and retains earlier reports. On distributions without `sage -python`, run the suite with a Python installation that provides `sage.all`.

This is a direct reference implementation. Radicals and elimination Gröbner bases over ZZ can be very expensive. It may also normalize at other primes, whose corrections become redundant after base change to O_K.

`max_steps=20` and `max_variables=24` are configurable resource guards. Exceeding them raises `NormalizationLimitError`; an intermediate ring is never returned as the normalization. Individual Gröbner basis or radical computations have no built-in timeout. For large examples, a separate Sage process with an external time or memory limit is recommended.

As with any CAS result, the certificate assumes the correctness of the exact Singular algorithms used. The implementation additionally checks the ideal quotient inclusion and the graph relations. The supplied regression tests do not replace independent review for research results.

### Singular restarts and memory use

`NormalizationBackendError` means that Singular restarted during an operation and the computation was interrupted. No normalization is returned. A common cause is the operating system terminating Singular because it ran out of memory; a backend error alone does not establish that cause. On Linux, check the kernel log for an out-of-memory message. Increasing `max_steps` or `max_variables` does not increase available memory.

The November 14 update avoids the memory-intensive syzygy path, simplifies the graph fractions, and improves the test element in integral basis mode. Large radical computations can still take substantial time. Restart your Sage session after updating the module so that Python imports the new code.

For the reported input `(t^2-s^3)^2*(1+s)-2` over the extension defined by `T^8+4*T^6+2*T^4+28*T^2+1` at p=2, the backend update completed the first enlargement, including the previously problematic ideal quotient and graph computation. The following radical computation was stopped to limit test time. A full normalization of this example has **not** yet been verified. The large example was not rerun after adding generator simplification.

## References

* [Greuel–Laplagne–Seelisch: Normalization of Rings](https://arxiv.org/abs/0904.3561)
* [Singular: primdecint.lib](https://github.com/Singular/Singular/blob/spielwiese/Singular/LIB/primdecint.lib)
* [Sage: The Singular interface](https://doc.sagemath.org/html/en/reference/interfaces/sage/interfaces/singular.html)
* [Sage: Number fields, integral bases, and prime decomposition](https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/number_field/number_field.html)
* [Stacks Project: Models and completion](https://stacks.math.columbia.edu/tag/0C2T)
