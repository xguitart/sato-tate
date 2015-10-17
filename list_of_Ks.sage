R.<x> = PolynomialRing(Rationals())

pols = [x^2 + 3, (x^2 + 1)*(x^2 - 2), (x^2 + 3)*(x^3 - 2), (x^2 + 2)*(x^4 + 17*x^2 + 68), (x^2 + 3)*(x^6 + 2), (x^2 + 1)*(x^2 - 2)*(x^4 - 3), (x^2 + 1)*(x^2 - 2)*(x^2 - 3)*(x^3 + 3*x - 2), (x^2 + 2)*(x^3 - 7*x + 7)*(x^4 +4*x^2 + 8*x +8), (x^2 + 2)*(x^2 + 11)*(x^3 - 4*x +4)*(x^4 +22*x + 22), (x^2 + 1)*(x^2 - 2), (x^2 + 1)*(x^2 - 2), (x^2 + 3)*(x^6 + 2), (x^2 + 1)*(x^2 - 3)*(x^3 +3*x^2 -1), (x^2 + 1)*(x^2 - 2)*(x^2 - 3), (x^2 + 3)*(x^6 + 2), (x^2 + 1)*(x^2 - 2)*(x^4 - 3), (x^2 + 3), (x^2 + 1)*(x^4 - 2), (x^2 + 3)*(x^3 -3*x +1), (x^2 + 1)*(x^2 - 2), (x^2 + 1)*(x^4 - 2), (x^2 + 2)*(x^2 + 3)*(x^3 - 9*x + 6), (x^2 + 3)*(x^3 - 2), (x^2 + 2)*(x^4 - 14*x^2 + 28*x - 14), (x^2 + 3)*(x^6 - 2), (x^2 + 2)*(x^3 + 5*x + 10)*(x^4 + 4*x^2 + 8*x + 2), (x^2 + 1)*(x^2 - 2), (x^2 + 1)*(x^2 - 2), (x^2 + 1)*(x^2 - 2), (x^4 + 5*x^2 + 5), (x^2 + 1)*(x^2 - 2), (x^2 - 2), (x^3 - 3*x + 1), (x^4 - 5*x^2 + 5), (x^2 - 7)*(x^3 - 7*x - 7), (x^2 + 1), (x^2 + 1)*(x^2 - 2), (x^2 + 3)*(x^3 - 2), (x^2 + 1)*(x^4 - 2), (x^2 + 3)*(x^6 + 2), (x^2 + 1), (x^2 + 1), (x^2 + 1), x^4 - x^3 + 5*x^2 - 5*x + 2];

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
    K.<g> = K.galois_closure()
    A.<y> = PolynomialRing(K)
    Ds = [1,-10,-163, -67, 19*43,-3*19,31*71, -51]
    for i in range(len(Ds)):
        for j in range(i+1,len(Ds)):
            g = (y^2 - Ds[i]*Ds[j])
            if len(g.roots()) > 0:
                print Ds[i],Ds[j]
                print p.factor()
