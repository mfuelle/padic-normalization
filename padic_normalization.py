"""Exact normalization of reduced affine plane models over Z_p.

Run with SageMath (not ordinary CPython).  See README.md for the precise
input contract, mathematics, examples, and limitations.  No p-adic rounding
is used by the algorithm.  The backend is Singular's integer arithmetic.
"""
from sage.all import ZZ, QQ, GF, PolynomialRing, NumberField, Zp, TermOrder, matrix, vector
from sage.interfaces.singular import singular


class NormalizationLimitError(RuntimeError):
    """A resource guard fired. No normalization is returned as certified."""


class NormalizationBackendError(RuntimeError):
    """Singular restarted during an operation; no result is certified."""


class CoefficientPresentationMap:
    """Verified coefficient inclusion, evaluated on unreduced representatives.

    This lightweight callable deliberately avoids Sage's quotient homset,
    which requests unavailable Groebner bases even with check=False.
    """
    def __init__(self, source, target, indices):
        self._source, self._target = source, target
        self._indices = tuple(indices)

    def domain(self):
        return self._source

    def codomain(self):
        return self._target

    def im_gens(self):
        B = self._target
        return tuple(B.element_class(B, B.cover_ring().gen(i), reduce=False)
                     for i in self._indices)

    def __call__(self, value):
        C, B = self._source, self._target
        if hasattr(value, 'parent') and value.parent() is C:
            value = value.lift()
        polynomial = C.cover_ring()(value)
        R = B.cover_ring()
        ev = C.cover_ring().hom([R.gen(i) for i in self._indices], R)
        return B.element_class(B, ev(polynomial), reduce=False)

    def __repr__(self):
        return 'Verified coefficient inclusion (evaluation on polynomial representatives)'


_SINGULAR_HELPERS = r'''
LIB "primdecint.lib";
proc pn_std(ideal I)
{
  ideal G = std(I);
  int i;
  for (i=1; i<=size(G); i++)
  { if (leadcoef(G[i])<0) { G[i]=-G[i]; } }
  attrib(G,"isSB",1);
  return(G);
}
'''


class _Engine:
    def __init__(self):
        self._load_helpers()

    def _load_helpers(self):
        # Expect's line-by-line evaluator otherwise waits at an unfinished
        # procedure declaration. This helper string contains no // comments.
        singular.eval(' '.join(_SINGULAR_HELPERS.splitlines()))
        self._session = singular._session_number

    def _execute(self, name, R, args, convert):
        # A deliberate restart between operations is safe: all input data
        # live in Sage and can be converted into the new Singular session.
        singular.pid()  # Start the process if the user called singular.quit().
        if self._session != singular._session_number:
            self._load_helpers()
        session = singular._session_number
        message = ('Singular restarted during %s in a polynomial ring with %s '
                   'variables. The interrupted computation is not certified; '
                   'no normalization has been returned.' % (name, R.ngens()))
        try:
            R._singular_().set_ring()
            result = getattr(singular, name)(*[a._singular_() for a in args])
            if session != singular._session_number:
                raise NormalizationBackendError(message)
            result = convert(result)
        except Exception as exc:
            if session != singular._session_number:
                raise NormalizationBackendError(message) from exc
            raise
        if session != singular._session_number:
            raise NormalizationBackendError(message)
        return result

    def call(self, name, R, *args):
        return self._execute(name, R, args,
                             lambda result: R.ideal([g.sage_poly(R) for g in result]))

    def std(self, I):
        return self.call('pn_std', I.ring(), I)

    def radical(self, I):
        return self.std(self.call('radicalZ', I.ring(), I))

    def colon(self, H, J):
        """Compute (H:J) by integer elimination, without integer syzygies.

        Starting with K=(1), replace K by (H intersect j*K)/j for each
        nonzero generator j of J. Since ZZ[variables] is a domain this is
        exactly K intersect (H:j). Division is exact, never saturation or
        extension of an ideal to QQ (which would lose integer torsion).
        """
        R = H.ring()
        aux = 'pn_intersection'
        while aux in R.variable_names():
            aux += '_'
        S = PolynomialRing(ZZ, names=(aux,)+R.variable_names(),
                           order=TermOrder('lex', 1)+TermOrder('degrevlex', R.ngens()))
        into = R.hom(S.gens()[1:], S)
        back = S.hom([R.zero()]+list(R.gens()), R)
        z = S.gen(0)
        K = R.ideal(1)
        for j in J.gens():
            if not j:
                continue
            G = self.std(S.ideal([z*into(h) for h in H.gens()] +
                                [(1-z)*into(j*k) for k in K.gens()]))
            intersection = [back(h) for h in G.gens() if h.degree(z) == 0]
            quotients = [R(h/j) for h in intersection]
            if any(h != j*q for h, q in zip(intersection, quotients)):
                raise ArithmeticError('Ideal intersection division was not exact.')
            K = self.std(R.ideal(quotients))
        if not self.subset(K*J, H):
            raise ArithmeticError('The computed ideal quotient failed verification.')
        return K

    def nf(self, f, G):
        R = G.ring()
        return self._execute('reduce', R, (f, G, ZZ(5)),
                             lambda result: result.sage_poly(R))

    def subset(self, I, J):
        G = self.std(J)
        return all(self.nf(f, G) == 0 for f in I.gens())

    def graph(self, I, d, numerators, new_names):
        """Kernel of z_i -> numerator_i/d, by integer elimination.

        Saturation is essential: the unsaturated graph has spurious
        components on d=0.  All computations here remain over ZZ.
        """
        R = I.ring()
        names = R.variable_names() + tuple(new_names)
        S = PolynomialRing(ZZ, names=names, order='degrevlex')
        aux = 'pn_inverse'
        while aux in names:
            aux += '_'
        E = PolynomialRing(ZZ, names=(aux,) + names,
                           order=TermOrder('lex', 1) + TermOrder('degrevlex', len(names)))
        toE = R.hom(E.gens()[1:1+R.ngens()], E)
        # Cancel common polynomial factors before introducing the inverse.
        # This preserves u/d in the total quotient ring: every divisor of
        # the nonzerodivisor d is itself a nonzerodivisor. Keeping removable
        # factors can cause a severe memory blow-up in integer elimination.
        fractions = [R.fraction_field()(u)/d for u in numerators]
        denominator = R.one()
        for q in fractions:
            denominator = denominator.lcm(R(q.denominator()))
        rels = [toE(f) for f in I.gens()]
        rels += [toE(R(q.denominator()))*E.gen(1+R.ngens()+i)-toE(R(q.numerator()))
                 for i, q in enumerate(fractions)]
        rels.append(1-E.gen(0)*toE(denominator))
        G = self.std(E.ideal(rels))
        back = E.hom([S.zero()] + list(S.gens()), S)
        out = [back(f) for f in G.gens() if f.degree(E.gen(0)) == 0]
        return self.std(S.ideal(out))


class PadicNormalization:
    """An exact finite affine presentation, interpreted over O_K.

    `integer_ring()` is an exact Sage quotient ring. `padic_ring()` is
    a finite-precision display/base-change convenience, not the certificate.
    Original variables are preserved, so the normalization map is explicit.
    """
    def __init__(self, p, original_ideal, ideal, images, steps, extension, certificate):
        self.p = ZZ(p)
        self.original_ideal = original_ideal
        self.ideal = ideal
        self.ambient_ring = ideal.ring()
        self.relations = tuple(ideal.gens())
        self.fraction_images = dict(zip(self.ambient_ring.variable_names(), images))
        self.steps = tuple(steps)
        self.extension = extension
        self.certificate = certificate
        self.coefficient_model = None
        self.integral_original_ideal = original_ideal

    def __repr__(self):
        base = 'Z_%s' % self.p if self.extension is None else 'O_K (a root of %s)' % self.extension
        return ('Normalization over %s\nVariables: %s\nRelations: %s\n'
                'Normalization steps: %s\nCertificate: %s' %
                (base, ', '.join(self.ambient_ring.variable_names()),
                 self.relations, len(self.steps), self.certificate))

    def integer_ring(self):
        return self.ambient_ring.quotient(self.ideal)

    def integer_map(self):
        """Return (original exact model, normalized model, inclusion)."""
        A = self.original_ideal.ring().quotient(self.original_ideal)
        B = self.integer_ring()
        phi = A.hom(list(B.gens())[:A.ngens()], B)
        return A, B, phi

    def padic_ring(self, prec=40):
        """Return a finite-precision presentation of the normalized O_K-algebra.

        Sage's support for arithmetic/Groebner bases in this quotient may
        be limited. Use the exact integer presentation for certificates.
        For integral-basis models the returned base_ring() is Zp, with
        coefficient variables retained; coefficient_map(prec) specifies
        the full O_K algebra structure. Fast paths use native O_K objects.
        """
        if self.coefficient_model is not None:
            return self.ring_over_integers(prec=prec)
        O = Zp(self.p, prec=prec)
        names = self.ambient_ring.variable_names()
        if self.extension is None:
            P = PolynomialRing(O, names=names)
            mapper = self.ambient_ring.hom(list(P.gens()), P)
        else:
            V = PolynomialRing(O, 'T')
            O = O.extension(V(self.extension), names=names[0])
            P = PolynomialRing(O, names=names[1:])
            mapper = self.ambient_ring.hom([P(O.gen())]+list(P.gens()), P)
        exact_relations = self.relations
        if self.extension is not None:
            ga = self.extension(self.ambient_ring.gen(0))
            exact_relations = [f.reduce([ga]) for f in exact_relations]
        # Discard only exact zeros, before conversion to finite precision.
        exact_relations = [f for f in exact_relations if f != 0]
        # Do not try to decide exact zero from a finite precision coefficient.
        return P.quotient(P.ideal([mapper(f) for f in exact_relations]))

    def coefficient_ring(self, prec=None):
        """Full O_F, or its coefficient base change to Z_p (finite precision).

        For general extensions this is an explicit quotient presentation,
        not Sage's native p-adic DVR class. Fast-path results have no such
        model; use integral_basis=True to request one for every extension.
        """
        if self.coefficient_model is None:
            raise ValueError('Use integral_basis=True to request a coefficient presentation.')
        C = self.coefficient_model['ambient_ring']
        if prec is None:
            return C.quotient(self.coefficient_model['ideal'])
        P = PolynomialRing(Zp(self.p, prec=prec), names=C.variable_names())
        ev = C.hom(list(P.gens()), P)
        return P.quotient(P.ideal([ev(h) for h in self.coefficient_model['ideal'].gens()]))

    def _relative_ring(self, I, O):
        """Move coefficient variables into O, removing exact zero relations first."""
        model = self.coefficient_model
        R = I.ring()
        positions = dict(zip(model['names'], O.gens()))
        curve_names = tuple(n for n in R.variable_names() if n not in positions)
        P = PolynomialRing(O, names=curve_names)
        variables = dict(zip(curve_names, P.gens()))
        images = [P(positions[n]) if n in positions else variables[n]
                  for n in R.variable_names()]
        toR = model['ambient_ring'].hom([R.gen(R.variable_names().index(n))
                                       for n in model['names']], R)
        E = _Engine()
        coefficient_ideal = E.std(R.ideal([toR(h) for h in model['ideal'].gens()]))
        exact = [E.nf(h, coefficient_ideal) for h in I.gens()]
        # Build via monomials: generic quotient coefficient rings need not
        # implement all coercion/homomorphism machinery used by Singular.
        rels = []
        for h in exact:
            if h == 0:
                continue
            value = P.zero()
            for powers, c in h.dict().items():
                term = P(c)
                for z, e in zip(images, powers):
                    term *= z**e
                value += term
            rels.append(value)
        return P.quotient(P.ideal(rels))

    def ring_over_integers(self, prec=None):
        """Normalized affine ring over the full coefficient-ring presentation.

        prec=None is exact over O_F. With prec=N return a FLATTENED
        Zp(p,N)-presentation, since nested quotient arithmetic over Zp is
        not implemented by Sage. coefficient_map(N) supplies its O_K
        algebra structure. The exact data remain authoritative.
        """
        if prec is not None:
            if self.coefficient_model is None:
                raise ValueError('Use integral_basis=True to request a coefficient presentation.')
            P = PolynomialRing(Zp(self.p, prec=prec),
                               names=self.ambient_ring.variable_names())
            ev = self.ambient_ring.hom(list(P.gens()), P)
            return P.quotient(P.ideal([ev(h) for h in self.relations if h != 0]))
        return self._relative_ring(self.ideal, self.coefficient_ring(prec))

    def coefficient_map(self, prec=None):
        """Return (coefficient ring, flattened normalized ring, structure map).

        Relations are checked in the exact integer model. The callable map
        evaluates polynomial representatives without unsupported reductions;
        it is a CoefficientPresentationMap, not a native Sage homset element.
        """
        C = self.coefficient_ring(prec)
        B = self.integer_ring() if prec is None else self.ring_over_integers(prec)
        names = self.ambient_ring.variable_names()
        indices = [names.index(n) for n in self.coefficient_model['names']]
        ev = self.coefficient_model['ambient_ring'].hom(
            [self.ambient_ring.gen(i) for i in indices], self.ambient_ring)
        J = self.ambient_ring.ideal([ev(h) for h in self.coefficient_model['ideal'].gens()])
        if not _Engine().subset(J, self.ideal):
            raise ArithmeticError('The coefficient structure map failed exact verification.')
        return C, B, CoefficientPresentationMap(C, B, indices)

    def normalization_map(self):
        """Exact flattened O_F-models (A, B, inclusion).

        Unlike integer_map(), the source uses the full ring of integers.
        Both Sage base rings are ZZ, with coefficient variables retained.
        Coefficient base change to Z_p gives the requested normalization map.
        """
        source = self.integral_original_ideal
        A = source.ring().quotient(source)
        B = self.integer_ring()
        ev = source.ring().hom(list(self.ambient_ring.gens())[:source.ring().ngens()],
                               self.ambient_ring)
        if not _Engine().subset(self.ambient_ring.ideal([ev(h) for h in source.gens()]), self.ideal):
            raise ArithmeticError('The normalization map failed exact verification.')
        return A, B, A.hom(list(B.gens())[:A.ngens()], B)


def _exact_polynomial(f, p, nvars):
    R = f.parent()
    if not hasattr(R, 'ngens') or R.ngens() != nvars:
        raise ValueError('Expected a polynomial in exactly %s variables.' % nvars)
    if R.base_ring() not in (ZZ, QQ):
        raise TypeError('Use exact ZZ or QQ coefficients; finite-precision p-adics are rejected.')
    if f == 0 or f.total_degree() <= 0:
        raise ValueError('Expected a nonzero, nonconstant polynomial.')
    if any(QQ(c).valuation(p) < 0 for c in f.coefficients()):
        raise ValueError('Every rational coefficient must be p-integral.')
    Q = PolynomialRing(QQ, names=R.variable_names())
    return Q(f)


def _prepare(f, p, extension):
    """Produce a reduced, Z-flat integer model and a test element d."""
    if extension is None:
        fq = _exact_polynomial(f, p, 2)
        den = ZZ(fq.denominator())
        R = PolynomialRing(ZZ, names=fq.parent().variable_names())
        fi = R(den*fq)
        content = ZZ(0)
        for c in fi.coefficients():
            content = content.gcd(ZZ(c))
        if content.valuation(p) > 0:
            raise ValueError('f is divisible by p: the model is not flat over Z_p.')
        fi = R(fi/content)
        q = fq.parent()(fi)
        if any(e > 1 for _, e in q.factor()):
            raise ValueError('The generic fibre is nonreduced; normalization of its reduction is not implicit.')
        partials = sorted([q.derivative(v) for v in q.parent().gens()
                           if q.derivative(v) != 0], key=lambda a: a.total_degree())
        for c in range(int(q.total_degree())+2):
            cand = sum((QQ(c)**i)*h for i,h in enumerate(partials))
            if cand and q.gcd(cand).total_degree() == 0:
                d = R(cand)
                return R.ideal(fi), d, None
        raise ArithmeticError('Could not construct a nonzerodivisor in the Jacobian ideal.')

    # Optional monogenic extensions: the first variable of f is the
    # coefficient generator a; the next two are curve coordinates x,y.
    fq = _exact_polynomial(f, p, 3)
    V = PolynomialRing(ZZ, 'T')
    if (not hasattr(extension, 'parent') or
            not hasattr(extension.parent(), 'base_ring') or
            extension.parent().base_ring() not in (ZZ, QQ)):
        raise TypeError('extension must have exact integer coefficients, not p-adic approximations.')
    try:
        g = V(extension)
    except (TypeError, ValueError):
        raise TypeError('extension must be a monic polynomial over ZZ.')
    if g.degree() < 2 or not g.is_monic():
        raise ValueError('extension must be monic of degree at least two.')
    eisenstein = (g[0].valuation(p) == 1 and
                  all(g[i] % p == 0 for i in range(g.degree())))
    unramified = g.change_ring(GF(p)).is_irreducible()
    if not g.is_irreducible():
        raise ValueError('extension must be irreducible over QQ_p (already reducible over QQ).')
    F = NumberField(g, 'aa')
    if not (eisenstein or unramified):
        primes = F.primes_above(p)
        if len(primes) != 1 or (primes[0].ramification_index() *
                                primes[0].residue_class_degree() != g.degree()):
            raise ValueError('extension is reducible over QQ_%s; choose one local factor.' % p)
    U = PolynomialRing(F, names=fq.parent().variable_names()[1:])
    ev = fq.parent().hom([U(F.gen())]+list(U.gens()), U)
    fu = ev(fq)
    if not fu or fu.total_degree() <= 0:
        raise ValueError('f must be nonconstant in the two curve coordinates.')
    # This sufficient condition guarantees Z-flatness of the global model.
    is_monic = False
    for v in U.gens():
        e = fu.degree(v)
        if e > 0 and fu.coefficient({v:e}) == 1:
            is_monic = True
    if not is_monic:
        raise ValueError('For extensions, f must be monic in at least one curve coordinate.')
    if any(e > 1 for _, e in fu.factor()):
        raise ValueError('The generic fibre is nonreduced.')
    # Integral coefficients are required here; rational denominators in
    # number-field coefficients would need a separately localized model.
    if fq.denominator() != 1:
        raise ValueError('For extensions, use integral coefficients in Z[a,x,y].')
    R = PolynomialRing(ZZ, names=fq.parent().variable_names())
    a = R.gen(0)
    ga = g(a)
    fi = R(fq)
    # Reduce powers of a so the monicity checked over F holds over Z[a].
    fi = fi.reduce([ga])
    choices = [fi.derivative(v) for v in R.gens()[1:]]
    choices = sorted([h for h in choices if h], key=lambda a: a.total_degree())
    for c in range(int(fu.total_degree())+2):
        h = sum(ZZ(c)**i*v for i,v in enumerate(choices))
        hu = ev(fq.parent()(h))
        if hu and fu.gcd(hu).total_degree() == 0:
            return R.ideal(ga, fi), g.derivative()(a)*h, g
    raise ArithmeticError('Could not construct a test element for the extension.')


def _integral_coefficient_model(I, d, g, p, max_variables):
    """Replace Z[a] by the full O_F using a multiplication-table presentation.

    Retain a,x,y and adjoin all basis elements (even 1) to avoid assuming
    any particular shape of the integral basis. No finite precision is used.
    """
    F = NumberField(g, 'aa')
    basis = tuple(F.integral_basis())
    n = F.degree()
    def row(v):
        v = F(v)
        return vector(QQ, [v[i] for i in range(n)])
    M = matrix(QQ, [row(b) for b in basis])
    inverse = M.inverse()
    def coords(v):
        return tuple(ZZ(c) for c in row(v)*inverse)
    names = list(I.ring().variable_names())
    basis_names = []
    counter = 0
    while len(basis_names) < n:
        name = 'ib%s' % counter
        counter += 1
        if name not in names:
            basis_names.append(name)
    if len(names)+n > max_variables:
        raise NormalizationLimitError('max_variables reached while adjoining the integral basis.')
    C = PolynomialRing(ZZ, names=(names[0],)+tuple(basis_names))
    a = C.gen(0)
    b = C.gens()[1:]
    def linear(v):
        return sum((c*z for c,z in zip(coords(v), b)), C.zero())
    rels = [C.one()-linear(F.one()), a-linear(F.gen())]
    rels += [b[i]*b[j]-linear(basis[i]*basis[j])
             for i in range(n) for j in range(i,n)]
    E = _Engine()
    J = E.std(C.ideal(rels))
    S = PolynomialRing(ZZ, names=tuple(names+basis_names))
    old = I.ring().hom(list(S.gens())[:len(names)], S)
    coeff = C.hom([S.gen(0)]+list(S.gens()[len(names):]), S)
    newI = E.std(S.ideal([old(h) for h in I.gens()]+[coeff(h) for h in J.gens()]))
    F0 = I.ring().fraction_field()
    images = [F0(v) for v in I.ring().gens()]
    images += [sum((F0(v[i])*F0(I.ring().gen(0))**i for i in range(n)), F0.zero())
               for v in basis]
    prime = F.primes_above(p)[0]
    model = {'ambient_ring':C, 'ideal':J, 'names':C.variable_names(),
             'number_field':F, 'basis':basis,
             'ramification_index':prime.ramification_index(),
             'residue_degree':prime.residue_class_degree()}
    # O_F is already normal. After inverting the curve Jacobian element h,
    # the model is smooth over O_F, hence normal. The factor g'(a) was only
    # needed for the original order Z[a]; retaining it here can make the
    # subsequent saturation needlessly expensive.
    h = I.ring()(d/g.derivative()(I.ring().gen(0)))
    # If the generic fibre is smooth, the relative Jacobian ideal may
    # contain a nonzero integer. Its principal open is normal too, and a
    # constant avoids expensive characteristic-zero elimination in radicalZ.
    derivatives = [old(q.derivative(v)) for q in I.gens()
                   for v in I.ring().gens()[1:] if q.derivative(v)]
    jacobian = E.std(newI+S.ideal(derivatives))
    constants = [ZZ(q) for q in jacobian.gens() if q and q.total_degree() == 0]
    test_element = S(min(abs(c) for c in constants)) if constants else old(h)
    return newI, test_element, images, model


def _with_vertical_component(H, max_variables):
    """Normalize p*f as Normalization(f) x F_p[x,y].

    e is the horizontal idempotent. Extra horizontal generators are extended
    by zero on the vertical component. Original x,y are kept on both factors.
    """
    R0 = H.original_ideal.ring()
    fi = H.original_ideal.gen(0)
    names = H.ambient_ring.variable_names()
    ename = 'vertical_split'
    while ename in names:
        ename += '_'
    if len(names)+1 > max_variables:
        raise NormalizationLimitError('max_variables reached while separating the vertical component.')
    S = PolynomialRing(ZZ, names=names+(ename,), order='degrevlex')
    e = S.gens()[-1]
    into = H.ambient_ring.hom(list(S.gens())[:-1], S)
    rels = [e**2-e, H.p*(1-e)]
    rels += [e*into(h) for h in H.relations]
    rels += [(1-e)*z for z in S.gens()[R0.ngens():-1]]
    E = _Engine()
    I = E.std(S.ideal(rels))
    F0 = R0.fraction_field()
    images = [F0(v) for v in R0.gens()]
    for name in names[R0.ngens():]:
        q = H.fraction_images[name]
        images.append(F0(H.p*q.numerator())/(H.p*q.denominator()+fi))
    images.append(F0(H.p)/(H.p+fi))
    return PadicNormalization(H.p, R0.ideal(H.p*fi), I, images,
                              list(H.steps)+[{'vertical_idempotent':ename}],
                              None, 'Product of a Grauert-Remmert-normalized ring and F_p[x,y]')


def normalize_affine_padic(f, p, extension=None, max_steps=20,
                           max_variables=24, verbose=False, integral_basis=False):
    """Normalize O_K[x,y]/(f), returning an exact affine presentation.

    Required input (base case): f in ZZ[x,y] or QQ[x,y], p-integral
    coefficients, squarefree primitive part over QQ, and p-content at most
    one. A p-factor is separated as a vertical component; p^2 is rejected.

    Optional extension: monic g in ZZ[T], irreducible over QQ_p (checked
    exactly using prime decomposition if needed). Input f in ZZ[a,x,y],
    monic in x or y; a denotes a root of g. General extensions automatically
    use a full integral-basis coefficient presentation. integral_basis=True
    requests that presentation for Eisenstein/unramified fast paths too.

    The returned integer presentation is interpreted by coefficient base
    change to O_K. The algorithm also handles reduced reducible inputs.
    max_steps/max_variables are safety guards, not claims of convergence.
    There is no hard per-Groebner-computation timeout; use process isolation
    for one. A guard raises an exception, never returns a partial answer.
    """
    p = ZZ(p)
    if integral_basis and extension is None:
        raise ValueError('integral_basis=True requires an extension polynomial.')
    if not p.is_prime():
        raise ValueError('p must be prime.')
    if max_steps < 0 or max_variables < 2:
        raise ValueError('Invalid resource guards.')
    if extension is None:
        fq = _exact_polynomial(f, p, 2)
        v = min(QQ(c).valuation(p) for c in fq.coefficients())
        if v > 1:
            raise ValueError('f is divisible by p^2: the ring is nonreduced.')
        if v == 1:
            H = normalize_affine_padic(fq/p, p, max_steps=max_steps,
                                       max_variables=max_variables, verbose=verbose)
            return _with_vertical_component(H, max_variables)
    I, d, extension = _prepare(f, p, extension)
    E = _Engine()
    I = E.std(I)
    original = I
    R0 = I.ring()
    F0 = R0.fraction_field()
    images = [F0(v) for v in R0.gens()]
    coefficient_model = None
    if extension is not None:
        g = extension
        simple = ((g[0].valuation(p) == 1 and
                   all(g[i] % p == 0 for i in range(g.degree()))) or
                  g.change_ring(GF(p)).is_irreducible())
        if integral_basis or not simple:
            I, d, images, coefficient_model = _integral_coefficient_model(
                I, d, g, p, max_variables)
    integral_original = I
    steps = []
    counter = 0
    while True:
        R = I.ring()
        if verbose:
            print('Step %s: %s variables, %s relations' % (len(steps), R.ngens(), len(I.gens())), flush=True)
        J = E.radical(I+R.ideal(d))
        # Generators lying in I impose no colon condition, since I is in H.
        # Reduce before elimination to avoid redundant colon conditions.
        residues = [E.nf(j, I) for j in J.gens()]
        Jmod = R.ideal([j for j in residues if j != 0])
        H = E.std(I+R.ideal([d*j for j in Jmod.gens()]))
        C = E.colon(H, Jmod)
        old = E.std(I+R.ideal(d))
        numerators = []
        for u in C.gens():
            v = E.nf(u, old)
            if v:
                numerators.append(v)
                old = E.std(old+R.ideal(v))
        if not numerators:
            result = PadicNormalization(p, original, I, images, steps, extension,
                                        'Grauert-Remmert: (d*J:J) = d*A')
            result.coefficient_model = coefficient_model
            result.integral_original_ideal = integral_original
            return result
        if len(steps) >= max_steps:
            raise NormalizationLimitError('max_steps reached; no certified result returned.')
        if R.ngens()+len(numerators) > max_variables:
            raise NormalizationLimitError('max_variables reached; no certified result returned.')
        new_names = []
        while len(new_names) < len(numerators):
            counter += 1
            name = 'z%s' % counter
            if name not in R.variable_names():
                new_names.append(name)
        ev = R.hom(images, F0)
        new_images = [ev(u)/ev(d) for u in numerators]
        old_d = d
        newI = E.graph(I, d, numerators, new_names)
        into = R.hom(list(newI.ring().gens())[:R.ngens()], newI.ring())
        # Verify the original equations and graph equations survived.
        checks = [into(h) for h in I.gens()]
        checks += [into(d)*newI.ring().gen(R.ngens()+i)-into(u)
                   for i,u in enumerate(numerators)]
        if not E.subset(newI.ring().ideal(checks), newI):
            raise ArithmeticError('The graph presentation failed verification.')
        steps.append({'variables':tuple(new_names), 'numerators':tuple(numerators),
                      'denominator':old_d})
        images += new_images
        d = into(d)
        I = newI
