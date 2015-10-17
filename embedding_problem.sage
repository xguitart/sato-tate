K0.<a> = NumberField(x^4 - x^3 + 5*x^2 - 5*x + 2)
Q = QuadraticForm(K0.trace_pairing([1,a,a^2,a^3]))
Disc = Q.disc()
# the Hasse-Witt invariant only ramifies at 5 and oo
for p in Disc.support():
    print "p = ", p, ": Hasse-Witt invariant = ", Q.hasse_invariant(p)

k = QuadraticField(-40)
Bp = QuaternionAlgebra(k(2), k(Disc))
Bm = QuaternionAlgebra(k(-2), k(Disc))
# we chech that (2, Disc)_k and (-2, Disc)_k are trivial
assert Bp.discriminant() == 1 and Bm.discriminant() == 1
