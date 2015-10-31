# We construct the five number fields K of Table 8 using the values of u and z of the table; v and alpha are computed using the expressions of Proposition 4.5.

R.<x> = PolynomialRing(Rationals())
u = 81/320; z = 1; s = 1; v = (1 - z^2*u)/(s^2*u)
K.<ru> = NumberField(x^2 - u)
R.<x> = PolynomialRing(K)
KK.<ra> = NumberField(x^2 - (1-z*ru)/2)
KK.<ra> = KK.absolute_field()
R.<x> = PolynomialRing(KK)
L.<rv> = NumberField(x^2 - v)
LL.<gL> = L.absolute_field()
R.<x> = PolynomialRing(LL)
K.<r40> = NumberField(x^2 + 40)
K.<gK> = K.absolute_field()
pol1 = gK.minpoly()

u = 81/320; z = -16/9; s = 1; v = (1 - z^2*u)/(s^2*u)
R.<x> = PolynomialRing(QQ)
KK.<ra> = NumberField(x^2 - (1-z*ru)/2)
KK.<ra> = KK.absolute_field()
R.<x> = PolynomialRing(KK)
L = KK
LL.<gL> = L.absolute_field()
R.<x> = PolynomialRing(LL)
K.<r40> = NumberField(x^2 + 40)
K.<gK> = K.absolute_field()
pol2 = gK.minpoly()

u = 4/17; z = 1; v = (u^3 - z^2)/(3*s^2)
R.<x> = PolynomialRing(Rationals())
K.<alpha> = NumberField(x^3 - 3*u/4*x - z/4)
R.<x> = PolynomialRing(K)
KK.<ra> = NumberField(x^2 - u)
KK.<ra> = KK.absolute_field()
R.<x> = PolynomialRing(KK)
L.<rv> = NumberField(x^2 - v)
LL.<gL> = L.absolute_field()
R.<x> = PolynomialRing(LL)
K.<r51> = NumberField(x^2 + 51)
K.<gK> = K.absolute_field()
pol3 = gK.minpoly()

u = 4/17; z = -5/4; v = (u^3 - z^2)/(3*s^2)
R.<x> = PolynomialRing(Rationals())
K.<alpha> = NumberField(x^3 - 3*u/4*x - z/4)
R.<x> = PolynomialRing(K)
KK.<ra> = NumberField(x^2 - u)
KK.<ra> = KK.absolute_field()
R.<x> = PolynomialRing(KK)
L.<rv> = NumberField(x^2 - v)
LL.<gL> = L.absolute_field()
R.<x> = PolynomialRing(LL)
K = LL
K.<gK> = K.absolute_field()
pol4 = gK.minpoly()

u = 19; z = 19/2; v = (u^3 - z^2)/(3*s^2)
R.<x> = PolynomialRing(Rationals())
K.<alpha> = NumberField(x^3 - 3*u/4*x - z/4)
R.<x> = PolynomialRing(K)
KK.<ra> = NumberField(x^2 - u)
KK.<ra> = KK.absolute_field()
R.<x> = PolynomialRing(KK)
L = KK;
LL.<gL> = L.absolute_field()
R.<x> = PolynomialRing(LL)
K.<r3> = NumberField(x^2 + 3)
K.<gK> = K.absolute_field()
pol5 = gK.minpoly()

pols = [pol1, pol2, pol3, pol4, pol5];

# We now check that the only quadratic subfield of k_0 contained in K is the one given in the fourth row of Table 8
for p in pols:
    R.<x> = PolynomialRing(Rationals())
    fs = p.factor()
    num_factors = len(fs)
    K.<gK> = NumberField(fs[0][0])
    R.<x> = PolynomialRing(K)
    expected_degree = fs[0][0].degree()
    for i in range(1,num_factors):
        R.<x> = PolynomialRing(K)
        NewK.<gNewK> = K.extension(fs[i][0])
        K.<gK> = NewK.absolute_field()
        expected_degree = expected_degree*fs[i][0].degree()
    assert K.degree() == expected_degree
    # K.<g> = K.galois_closure()
    A.<y> = PolynomialRing(K)
    Ds = [1,-10,-163, -67, 19*43,-3*19, -51]
    for i in range(len(Ds)):
        for j in range(i+1,len(Ds)):
            g = (y^2 - Ds[i]*Ds[j])
            if len(g.roots()) > 0:
                print Ds[i],Ds[j]
                # print p.factor()

