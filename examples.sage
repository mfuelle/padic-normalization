from padic_normalization import normalize_affine_padic

# 1. Nichtnormale spezielle Faser ueber Z_3.
R.<x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2 - 9*(x^3-x), 3, verbose=True)
print(N)
print(N.fraction_images)
B = N.padic_ring(prec=40)
print(B)

# 2. Dein urspruengliches Beispiel: K = Q_3(sqrt(3)).
V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2 - 3*(x^3-x), 3, extension=T^2-3,
                           verbose=True)
print(N)
print(N.fraction_images)  # z1 = y/a
B = N.padic_ring(prec=40)
print(B)

# 3. Eine Kuspe: Normalisierung auch auf der generischen Faser.
R.<x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2-x^3, 5)
print(N)
A0, B0, phi = N.integer_map()
print(phi)

# 4. Reduzibel: eine affine Praesentation mit einem Idempotenten.
N = normalize_affine_padic(x*y, 3)
print(N)

# 5. Weder Eisenstein noch irreduzibel mod 2: K = Q_2(sqrt(5)).
# Z_2[a] ist nicht maximal; (1+a)/2 ist ganz.
V.<T> = PolynomialRing(ZZ)
R.<a,x,y> = PolynomialRing(ZZ)
N = normalize_affine_padic(y^2-x^3, 2, extension=T^2-5)
print(N.coefficient_model['basis'])
print(N.fraction_images)
O = N.coefficient_ring()       # Exakte Praesentation des vollen O_F.
B0 = N.ring_over_integers()
A0, B0, phi = N.normalization_map()
print(phi)
B = N.padic_ring(prec=40)      # Gemeinsam ueber Z_2, inklusive O_K-Relationen.
print(B)
C, B, iota = N.coefficient_map(prec=40)
print(iota)

# 6. Gemischt verzweigt: e=f=2; a = sqrt(2)+sqrt(5).
N = normalize_affine_padic(y^2-x, 2, extension=T^4-14*T^2+9)
print(N.coefficient_model['ramification_index'])  # 2
print(N.coefficient_model['residue_degree'])      # 2
print(N.ring_over_integers())
