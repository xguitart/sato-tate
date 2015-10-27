R<x> := PolynomialRing(Rationals());
K<a> := NumberField(x^4 - x^3 + 5*x^2 - 5*x + 2);
ms := [];
for x in [1, a, a^2, a^3] do
	for y in [1, a, a^2, a^3] do
		Append(~ms, Trace(x*y));
end for;
end for;
M := Matrix(4,4,ms);
Q := QuadraticForm(M);
// the Hasse-Witt (or Hasse-Minkowski) invariant only ramifies at 2 and oo
HasseMinkowskiInvariants(M);

// Alternatively, we can diagonalize the quadratic form
R:=Parent(Q);
x1 := R.1; x2 := R.2; x3 := R.3; x4 := R.4;
d := Evaluate(Q,[x1-x2-80*x3+56*x4,4*x2-13*x3-117*x4,-37*x3+9*x4,-26*x4]);
print "Coefficients of the diagonal form", Coefficients(d);
// The result is
//4*x1^2 - 148*x2^2 + 32708*x3^2 - 4420*x4^2
// Therefore the Witt invariant is the element in the Brauer group given by
// (4,-148)(4,32708)(4,-4420)(-148,32708)(-148,-4420)(32708,-4420)
// One easily checks that this ramifies at 2 and oo

Disc := Discriminant(RingOfIntegers(K));
k := NumberField(x^2 + 40);
Bp := QuaternionAlgebra<k | k!2, k!Disc>;
Bm := QuaternionAlgebra< k | k!(-2), k!Disc>;
  // we chech that (2, Disc)_k and (-2, Disc)_k are trivial
assert Norm(Discriminant(Bp)) eq 1 and Norm(Discriminant(Bm)) eq 1;
