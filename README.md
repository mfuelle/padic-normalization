# p-adic Normalization

Exact normalization of affine curve models over p-adic integer rings using
SageMath and Singular. The implementation computes finite affine presentations
and explicit generator images using an iterative Grauert–Remmert algorithm.
General finite extensions are supported through exact local irreducibility
checks and integral bases, including mixed ramification and nonmaximal initial
orders.

**Created with OpenAI Codex.** The implementation, documentation, and regression
tests were developed with Codex from user-supplied mathematical requirements.

## Quick start

Run the following in SageMath from this repository's directory:

```sage
from padic_normalization import normalize_affine_padic

V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)

# K = Q_2(sqrt(5)); Z_2[a] is not the full ring of integers.
N = normalize_affine_padic(y^2-x^3, p=2, extension=T^2-5)

print(N.relations)        # Exact defining relations over ZZ
print(N.fraction_images)  # Generator images in the original total quotient ring

C, B, iota = N.coefficient_map(prec=40)
```

Here `C` presents the full coefficient ring and `iota` gives its structure map
into the normalized algebra `B`. For general extensions, the p-adic output is
a flattened presentation over `Zp` with an explicit `O_K`-algebra structure,
rather than a native Sage p-adic DVR coefficient class. The exact integer
relations remain authoritative; finite-precision output is a convenience.

## Scope

- For `K = QQ_p`, use a nonconstant polynomial in two variables over `ZZ` or
  `QQ`, with p-integral coefficients and squarefree primitive part. A single
  factor of `p` is supported as a vertical component.
- For extensions, use a monic polynomial in `ZZ[T]` that is irreducible over
  `QQ_p`. The curve equation must have integer coefficients in `(a,x,y)`, be
  monic in a curve coordinate, and have reduced generic fibre.
- Arbitrary finite-precision p-adic input coefficients are not supported.
  Quotient-ring arithmetic in Sage has limitations; the flattened exact
  integer presentation is the computational reference.

## Documentation and tests

- [Detailed user guide](USER_GUIDE.md)
- [Implementation and mathematics](IMPLEMENTATION.md)
- [Examples](examples.sage)
- [Changelog](CHANGELOG.md)
- [Test report and source hashes](test_results.txt)

All 21 regression tests passed under SageMath 10.10 through SageMathCell on
November 13, 2026, using the complete source files. To run them locally:

```sh
sage -python test_padic_normalization.py
```

SageMath and its Singular backend are required; ordinary Python alone is not
sufficient.
