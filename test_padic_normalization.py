"""Run: sage -python test_padic_normalization.py

Tests compare with independently known normalizations, exercise multiple
iterations, and verify rejection of invalid input. They use exact arithmetic.
"""
import unittest
from sage.all import ZZ, QQ, Zp, PolynomialRing
from padic_normalization import normalize_affine_padic, NormalizationLimitError


class NormalizationTests(unittest.TestCase):
    def setUp(self):
        self.R = PolynomialRing(ZZ, names=('x', 'y'))
        self.x, self.y = self.R.gens()

    def test_smooth(self):
        N = normalize_affine_padic(self.y**2-self.x, 3)
        self.assertEqual(len(N.steps), 0)
        self.assertEqual(N.padic_ring().base_ring().prime(), 3)

    def test_cusp(self):
        N = normalize_affine_padic(self.y**2-self.x**3, 3)
        B = N.integer_ring()
        x, y, t = B.gens()
        self.assertEqual(t**2, x)
        self.assertEqual(t**3, y)
        A, B, phi = N.integer_map()
        self.assertEqual(phi(A.gen(0)), B.gen(0))
        self.assertEqual(phi(A.gen(1)), B.gen(1))

    def test_vertical_non_normality(self):
        from sage.all import FunctionField
        from padic_normalization import _Engine
        N = normalize_affine_padic(self.y**2-9*(self.x**3-self.x), 3)
        E = _Engine()
        R = N.ambient_ring
        # y/3 belongs to B iff y belongs to 3B. Avoid assumptions about
        # which integral linear combination of generators represents it.
        self.assertEqual(E.nf(R.gen(1), E.std(N.ideal+R.ideal(3))), 0)
        # Conversely every returned generator must lie in Z[u,t]/(t^2-u^3+u).
        F = FunctionField(QQ, 'u')
        u = F.gen()
        P = PolynomialRing(F, 'v')
        v = P.gen()
        K = F.extension(v**2-u**3+u, names='t')
        t = K.gen()
        ev = N.original_ideal.ring().hom([K(u), 3*t], K)
        for q in N.fraction_images.values():
            value = ev(q.numerator())/ev(q.denominator())
            for c in value.list():
                self.assertEqual(c.denominator().degree(), 0)
                polynomial = c.numerator()/c.denominator()[0]
                self.assertTrue(all(a in ZZ for a in polynomial.coefficients()))

    def test_reducible(self):
        N = normalize_affine_padic(self.x*self.y, 3)
        x, y, e = N.integer_ring().gens()
        self.assertEqual(e**2, e)
        self.assertNotEqual(e, 0)
        self.assertNotEqual(e, 1)
        self.assertEqual(x*e, 0)
        self.assertEqual(y*(1-e), 0)

    def test_multiple_steps(self):
        # Z[t^2,t^5] normalizes to Z[t].
        N = normalize_affine_padic(self.y**2-self.x**5, 5)
        self.assertGreaterEqual(len(N.steps), 2)
        B = N.integer_ring()
        x, y = B.gens()[:2]
        candidates = list(B.gens()[2:])
        self.assertTrue(any(t**2 == x and (t**5 == y or t**5 == -y)
                            for t in candidates))

    def test_rational_p_integral_input(self):
        R = self.R.change_ring(QQ)
        x, y = R.gens()
        N = normalize_affine_padic((y**2-x**3)/2, 3)
        xx, yy, t = N.integer_ring().gens()
        self.assertEqual(t**2, xx)
        self.assertEqual(t**3, yy)

    def test_ramified_extension(self):
        V = PolynomialRing(ZZ, 'T')
        T = V.gen()
        R = PolynomialRing(ZZ, names=('a', 'x', 'y'))
        a, x, y = R.gens()
        N = normalize_affine_padic(y**2-3*(x**3-x), 3, extension=T**2-3)
        aa, xx, yy, z = N.integer_ring().gens()
        self.assertEqual(aa*z, yy)
        self.assertEqual(z**2, xx**3-xx)
        self.assertEqual(N.fraction_images['z1'], y/a)
        self.assertEqual(N.padic_ring().base_ring().prime(), 3)

    def test_unramified_extension(self):
        V = PolynomialRing(ZZ, 'T')
        T = V.gen()
        R = PolynomialRing(ZZ, names=('a', 'x', 'y'))
        a, x, y = R.gens()
        N = normalize_affine_padic(y**2-x**3, 3, extension=T**2+1)
        aa, xx, yy, z = N.integer_ring().gens()
        t = aa*z
        self.assertEqual(t**2, xx)
        self.assertEqual(t**3, yy)
        self.assertEqual(N.padic_ring().base_ring().prime(), 3)

    def test_invalid_inputs(self):
        for f,p in [(self.y**2,3), (9*self.y**2-9*self.x,3),
                    (self.y**2-self.x,9), (self.R.zero(),3)]:
            with self.assertRaises(ValueError):
                normalize_affine_padic(f,p)
        P = PolynomialRing(Zp(3), names=('x','y'))
        with self.assertRaises(TypeError):
            normalize_affine_padic(P.gen(1)**2-P.gen(0),3)
        P = self.R.change_ring(QQ)
        with self.assertRaises(ValueError):
            normalize_affine_padic(P.gen(1)**2-P.gen(0)/3,3)

    def test_vertical_component(self):
        N = normalize_affine_padic(3*(self.y**2-self.x**3),3)
        B = N.integer_ring()
        x,y,t,e = B.gens()
        self.assertEqual(e**2,e)
        self.assertNotEqual(e,0)
        self.assertNotEqual(e,1)
        self.assertEqual(3*(1-e),0)
        self.assertEqual((1-e)*t,0)
        self.assertEqual(e*(t**2-x),0)
        self.assertEqual(e*(t**3-y),0)
        N.integer_map()  # Also checks that the original equation maps to zero.
        # Every returned relation must hold in the original total quotient ring.
        F = next(iter(N.fraction_images.values())).parent()
        ev = N.ambient_ring.hom(list(N.fraction_images.values()),F)
        for r in N.relations:
            self.assertEqual(N.original_ideal.reduce(ev(r).numerator()),0)

    def test_guards(self):
        with self.assertRaises(NormalizationLimitError):
            normalize_affine_padic(self.y**2-self.x**3,3,max_steps=0)
        with self.assertRaises(NormalizationLimitError):
            normalize_affine_padic(self.y**2-self.x**3,3,max_variables=2)


class GeneralExtensionTests(unittest.TestCase):
    """General defining polynomials, including genuinely mixed extensions."""

    def setUp(self):
        self.V = PolynomialRing(ZZ, 'T')
        self.T = self.V.gen()
        self.R = PolynomialRing(ZZ, names=('a', 'x', 'y'))
        self.a, self.x, self.y = self.R.gens()

    def _smooth(self, g, **kwds):
        return normalize_affine_padic(self.y**2-self.x, 2,
                                      extension=g, **kwds)

    def _zero(self, N, value):
        # Keep algebraic verification in the flattened exact ZZ presentation;
        # Sage need not implement Groebner bases over nested quotient rings.
        from padic_normalization import _Engine
        self.assertEqual(_Engine().nf(value, N.ideal), 0)

    def _coefficient_element(self, N, value):
        from sage.all import matrix, vector
        model = N.coefficient_model
        F = model['number_field']
        n = F.degree()
        basis = model['basis']
        M = matrix(QQ, [[F(b)[i] for i in range(n)] for b in basis])
        row = vector(QQ, [F(value)[i] for i in range(n)])
        coefficients = tuple(ZZ(c) for c in row*M.inverse())
        R = N.ambient_ring
        variables = dict(zip(R.variable_names(), R.gens()))
        return sum((c*variables[name] for c, name in
                    zip(coefficients, model['names'][1:])), R.zero())

    def _check_coefficient_table(self, N):
        model = N.coefficient_model
        F = model['number_field']
        C = model['ambient_ring']
        ev = C.hom([F.gen()] + list(model['basis']), F)
        self.assertTrue(all(ev(h) == 0 for h in model['ideal'].gens()))
        self._zero(N, self._coefficient_element(N, F.one())-1)
        self._zero(N, self._coefficient_element(N, F.gen())-
                   N.ambient_ring.gen(0))

    def test_nonmaximal_unramified_quadratic_at_two(self):
        N = self._smooth(self.T**2-5)
        self.assertEqual(len(N.steps), 0)
        self.assertEqual(N.coefficient_model['ramification_index'], 1)
        self.assertEqual(N.coefficient_model['residue_degree'], 2)
        F = N.coefficient_model['number_field']
        beta = self._coefficient_element(N, (1+F.gen())/2)
        self._zero(N, beta**2-beta-1)
        self._zero(N, 2*beta-1-N.ambient_ring.gen(0))
        self._check_coefficient_table(N)

        O = N.coefficient_ring()
        self.assertEqual(O.base_ring(), ZZ)
        B = N.ring_over_integers()
        self.assertEqual(B.base_ring(), O)
        self.assertEqual(B.ngens(), 2)
        A, B, phi = N.normalization_map()
        self.assertEqual(A.base_ring(), B.base_ring())
        self.assertEqual(A.ngens(), N.integral_original_ideal.ring().ngens())
        self.assertEqual(phi(A.gen(1)), B.gen(1))
        self.assertEqual(phi.domain(), A)
        self.assertEqual(phi.codomain(), B)
        Op = N.coefficient_ring(prec=20)
        self.assertEqual(Op.base_ring().prime(), 2)
        Bp = N.padic_ring(prec=20)
        self.assertEqual(Bp.base_ring().prime(), 2)
        self.assertEqual(Bp.ngens(), N.ambient_ring.ngens())

    def test_nonmaximal_ramified_quadratic_at_two(self):
        N = self._smooth(self.T**2-8)
        self.assertEqual(len(N.steps), 0)
        self.assertEqual(N.coefficient_model['ramification_index'], 2)
        self.assertEqual(N.coefficient_model['residue_degree'], 1)
        F = N.coefficient_model['number_field']
        pi = self._coefficient_element(N, F.gen()/2)
        self._zero(N, pi**2-2)
        self._zero(N, 2*pi-N.ambient_ring.gen(0))
        self._check_coefficient_table(N)

    def test_coefficient_structure_maps(self):
        N = self._smooth(self.T**2-5)
        for prec in (None, 20):
            C, B, iota = N.coefficient_map(prec)
            self.assertIs(iota.domain(), C)
            self.assertIs(iota.codomain(), B)
            # Use representatives so no finite-precision quotient reduction
            # or exact equality decision is delegated to Sage.
            q = C.cover_ring().gen(0)
            image = iota(q).lift()
            self.assertEqual(image, B.cover_ring().gen(0))
            self.assertEqual(len(iota.im_gens()), C.ngens())

    def test_maximal_order_non_eisenstein_at_two(self):
        N = self._smooth(self.T**2+1)
        self.assertEqual(len(N.steps), 0)
        a = N.ambient_ring.gen(0)
        pi = a-1
        self._zero(N, pi**2+2*pi+2)
        self.assertEqual(N.coefficient_model['ramification_index'], 2)
        self._check_coefficient_table(N)

    def test_mixed_quartic_at_two(self):
        # a = sqrt(2) + sqrt(5): an unramified quadratic times a ramified one.
        N = self._smooth(self.T**4-14*self.T**2+9)
        self.assertEqual(len(N.steps), 0)
        self.assertEqual(N.coefficient_model['ramification_index'], 2)
        self.assertEqual(N.coefficient_model['residue_degree'], 2)
        F = N.coefficient_model['number_field']
        a = F.gen()
        pi = self._coefficient_element(N, (a**3-11*a)/6)
        beta = self._coefficient_element(N, (1+(17*a-a**3)/6)/2)
        self._zero(N, pi**2-2)
        self._zero(N, beta**2-beta-1)
        self._zero(N, pi+2*beta-1-N.ambient_ring.gen(0))
        self._check_coefficient_table(N)
        self.assertEqual(N.padic_ring(prec=20).ngens(), N.ambient_ring.ngens())

    def test_cusp_over_full_integer_ring(self):
        from padic_normalization import _Engine
        N = normalize_affine_padic(self.y**2-self.x**3, 2,
                                  extension=self.T**2-5)
        self.assertGreaterEqual(len(N.steps), 1)
        self._check_coefficient_table(N)
        R = N.ambient_ring
        a, x, y = R.gens()[:3]
        E = _Engine()
        # t=y/x belongs to B exactly when y lies in xB, since x is regular.
        self.assertEqual(E.nf(y, E.std(N.ideal+R.ideal(x))), 0)

        # Independently check B is contained in O_F[t] by evaluating every
        # returned generator at x=t^2, y=t^3. Together these prove B=O_F[t].
        F = N.coefficient_model['number_field']
        U = PolynomialRing(F, 't')
        t = U.gen()
        Q = U.fraction_field()
        old = N.original_ideal.ring()
        ev = old.hom([Q(F.gen()), Q(t**2), Q(t**3)], Q)
        OF = F.maximal_order()
        for q in N.fraction_images.values():
            value = ev(q.numerator())/ev(q.denominator())
            self.assertEqual(value.denominator().degree(), 0)
            for c in U(value).coefficients():
                OF(c)  # Raises if a polynomial coefficient is not integral.
        A, B, phi = N.normalization_map()
        self.assertEqual(A.ngens(), N.integral_original_ideal.ring().ngens())
        self.assertGreater(B.ngens(), A.ngens())
        self.assertEqual(phi.domain(), A)

    def test_reject_locally_or_globally_reducible_extension(self):
        for g in (self.T**2-17, self.T**2-1):
            with self.assertRaises(ValueError):
                self._smooth(g)

    def test_integral_basis_name_collision(self):
        R = PolynomialRing(ZZ, names=('a', 'ib0', 'y'))
        a, x, y = R.gens()
        N = normalize_affine_padic(y**2-x, 2, extension=self.T**2-5)
        names = N.ambient_ring.variable_names()
        self.assertEqual(len(set(names)), len(names))
        self.assertNotIn('ib0', N.coefficient_model['names'])
        self.assertEqual(N.ring_over_integers().ngens(), 2)
        self._check_coefficient_table(N)

    def test_integral_basis_guard(self):
        with self.assertRaises(NormalizationLimitError):
            self._smooth(self.T**2-5, max_variables=4)

    def test_explicit_integral_basis_on_fast_path(self):
        N = self._smooth(self.T**2-2, integral_basis=True)
        self.assertIsNotNone(N.coefficient_model)
        self.assertEqual(len(N.steps), 0)
        self.assertEqual(N.coefficient_model['ramification_index'], 2)
        self._check_coefficient_table(N)
        self.assertEqual(N.ring_over_integers().ngens(), 2)
        self.assertEqual(N.padic_ring(prec=20).base_ring().prime(), 2)


if __name__ == '__main__':
    unittest.main()
