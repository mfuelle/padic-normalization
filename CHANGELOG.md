# Changes on November 14, 2026

- `verbose=True` now prints every normalization adjunction before graph elimination, in both current and original variables, with immediate flushing. Completion, coefficient basis elements, the vertical idempotent, and final certification are also reported.
- Added exact removal of individually redundant normalization generators after every enlargement (`simplify=True`, the default), plus `N.simplified()` for postprocessing. Original and coefficient variables are retained. Removal expressions are recorded in `N.simplifications`.
- Added `max_new_generators` to bound the number of fractions adjoined per step. This can increase the number of steps; `max_variables` still applies before simplification. These options do not promise a globally minimal generating set.
- Added regression coverage for nonlinear substitutions, integer nonunits, multi-step cusp simplification, postprocessing, preserved embeddings, and one-generator batches.
- Replaced integer syzygy ideal quotients with exact integer intersection elimination and division. The previous route exhausted memory on the reported degree-eight extension example.
- Cancelled common polynomial factors before graph saturation, retaining the full relation kernel with smaller elimination inputs.
- In integral basis mode, removed the unnecessary coefficient-order derivative factor from the test element and prefer a nonzero integer from the relative Jacobian ideal when available.
- Added `NormalizationBackendError` for a Singular restart during an operation. Helpers are reloaded after a restart between operations.
- Added tests for integer torsion in ideal quotients and controlled Singular restarts. Large examples remain subject to radical-computation time and memory costs.

# Changes on November 13, 2026

- General monic extension polynomials over ZZ are now allowed, provided they are irreducible over Q_p. This is checked exactly using prime ideal decomposition in the number field. In particular, this includes extensions of Q_2 whose defining polynomial is neither Eisenstein nor irreducible modulo 2.
- A full integral basis and multiplication table extend the nonmaximal order ZZ[a] when needed. Curve normalization then starts over the full ring of integers.
- New methods: `coefficient_ring`, `ring_over_integers`, `normalization_map`, and `coefficient_map`; `integral_basis=True` requests the new output path for the existing extension types as well.
- The general p-adic output path uses a flattened presentation over Z_p with an explicit O_K-algebra structure. It does not require a native Sage DVR class for the general extension polynomial. Sage's limitations in quotient arithmetic are documented in the [user guide](USER_GUIDE.md).
- The ideal quotient is computed using sequential syzygy computations. Redundant generators are first removed modulo the relation ideal. This avoids the large block matrices that caused Singular crashes with the new test models.
- The existing input restrictions for the curve equation remain: integer coefficients in a,x,y, monicity in a curve coordinate, and a reduced generic fibre. Arbitrary finite-precision p-adic coefficients are still rejected.
- New regression tests check nonmaximal orders, mixed ramification, normalization of a cusp over the full ring of integers, structure maps, local irreducibility, and resource limits. The executed test suite is documented in `test_results.txt`.

The names and choices of additional generators are not canonical. For example, `y^2=9(x^3-x)` may now produce `z=2y/3`. The familiar normalization coordinate is then `t=y-z=y/3`; this presents the same normalization.
