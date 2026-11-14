from padic_normalization import normalize_affine_padic

# 1. Nonnormal special fibre over Z_3.
R.<x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2 - 9*(x^3-x), 3, verbose=True)
print(N)
print(N.fraction_images)
B = N.padic_ring(prec=40)
print(B)

# 2. The original example: K = Q_3(sqrt(3)).
V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2 - 3*(x^3-x), 3, extension=T^2-3,
                           verbose=True)
print(N)
print(N.fraction_images)  # z1 = y/a
B = N.padic_ring(prec=40)
print(B)

# 3. A cusp: normalization is also needed on the generic fibre.
R.<x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2-x^3, 5)
print(N)
A0, B0, phi = N.integer_map()
print(phi)

# 4. Reducible: an affine presentation with an idempotent.
N = normalize_affine_padic(x*y, 3)
print(N)

# 5. Neither Eisenstein nor irreducible mod 2: K = Q_2(sqrt(5)).
# Z_2[a] is not maximal; (1+a)/2 is integral.
V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2-x^3, 2, extension=T^2-5)
print(N.coefficient_model['basis'])
print(N.fraction_images)
O = N.coefficient_ring()       # Exact presentation of the full O_F.
B0 = N.ring_over_integers()
A0, B0, phi = N.normalization_map()
print(phi)
B = N.padic_ring(prec=40)      # Flattened over Z_2, including the O_K relations.
print(B)
C, B, iota = N.coefficient_map(prec=40)
print(iota)

# 6. Mixed ramification: e=f=2; a = sqrt(2)+sqrt(5).
N = normalize_affine_padic(y^2-x, 2, extension=T^4-14*T^2+9)
print(N.coefficient_model['ramification_index'])  # 2
print(N.coefficient_model['residue_degree'])      # 2
print(N.ring_over_integers())

# 7. Remove older generators after each step; adjoin at most one at a time.
R.<x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2-x^5, 5, simplify=True,
                           max_new_generators=1, max_steps=100, verbose=True)
print(N.ambient_ring.variable_names())  # x, y, and one normalization parameter
print(N.simplifications)
# For an existing result computed with simplify=False:
# small = N.simplified(verbose=True)
