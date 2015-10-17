//////////////////////////////////////////////////
// Script for the computations of Lemma 5.7
//////////////////////////////////////////////////

//////////////////////////////////////////////////
////         Some utils
//////////////////////////////////////////////////

// The cocycle c as given by MAGMA takes values in {0,1}; this function evaluates the cocycle at elements s,t of G returning the result as an element of {-1,+1}
ev_c := function(c,s,t)
    if c(<s,t>)[1] eq 0 then return 1; end if;
    if c(<s,t>)[1] eq 1 then return -1; end if;
end function;

// Given an element g an a list ee, returns the index of g in ee (or -1 if g is not in ee)
find_index := function(g,ee)
    for i in [1..#ee] do
        if ee[i] eq g then return i;end if;
    end for;
    return -1;
end function;

// given a subgroup H of S4, returns true if H is contained in A4
is_in_A4 := function(H)
    for h in Set(H) do
        if IsEven(h) eq false then return false; end if;
    end for;
    return true;
end function;

// Given a cocycle c and an element u of G, checks whether c is c-regular
is_c_regular := function(G,c,u)
    for b in Centralizer(G,u) do
        if c(<b,u>) ne c(<u,b>) then 
            return false; 
        end if;
    end for;
    return true;
end function;

/////////////////////////////////////////////////
////      This is the main computation
/////////////////////////////////////////////////
SetSeed(0);
R<t> := PolynomialRing(Rationals());

// We compute the cohomology of S_4
G := SymmetricGroup(4);
L := SubgroupLattice(G);
G_e := IndexedSet(Set(G)); // a list containing the elements of G
Irrs := AbsolutelyIrreducibleModules(G, GF(2));
MM := Irrs[1]; // the G-module Z/2Z
CM := CohomologyModule(G, MM); 
H2 := CohomologyGroup(CM,2); // the cohomology of G with values in Z/2Z
C := IndexedSet(Set(H2)); // these are the four 2-cocycles of H^2(G,Z/2Z), let us name them as c1, c2, c3, c4 (c1 is the trivial class)
c1 := TwoCocycle(CM, C[1]);
c2 := TwoCocycle(CM, C[2]);
c3 := TwoCocycle(CM, C[3]);
c4 := TwoCocycle(CM, C[4]);

// now we do a loop for finding out which are the two non-symmetric classes
c_ns := [];
for c in [c1, c2, c3, c4] do
    b := true;
    for s in G_e do 
	for t in G_e do
	    if s*t eq t*s then
	         b := b and (c(<s,t>) eq c(<t,s>));
            end if;
        end for;
    end for;
    if b eq true then //the cocycle is symmetric
	continue;
    end if;
    if b eq false then // the cocycle is not symmetric
	 Append(~c_ns,c);
    end if;           
end for;
cns1 := c_ns[1]; cns2 := c_ns[2]; // these are the two non-symmetric cocycles

////////////////////////////////////////////////////////////////////////////
// Here we identify cns1 and cns2 in terms of the group extension they correspond
GL23 := SmallGroup(48,28); // GL(2,Z/3Z)
GL23_e := IndexedSet(Set(GL23));
H := SubgroupLattice(GL23)[2]; // this is a normal subgroup of index 2 with quotient isomorphic to S4
S4, pi := quo<GL23 | H>;
S4_e := IndexedSet(Set(S4));
// given s in S4, returns a lift of it to GL(2,Z/3Z); if s = 1 we return the identity.
lift := function(s)
    if s eq Identity(S4) then return Identity(GL23); end if;
    for i in [1..48] do
        if pi(GL23_e[i]) eq s then return GL23_e[i]; end if;
    end for;
    return -1;
end function;
// we construct the cocycle corresponding to GL(2,Z/3Z)
c_GL23 := function(s,t)
    if lift(s)*lift(t)*lift(s*t)^-1 eq Identity(GL23) then return 1; else return -1; end if;
end function;
// c_GL23 is non-trivial on a subgroup of order 2 and length 6 of S4
LL := SubgroupLattice(S4);
C2 := LL[3];
C2_e := IndexedSet(Set(C2));
sigma := C2_e[1];
assert C2_e[2] eq Id(C2);// Sigma is the non-trivial element
s := L[3].1; assert Order(L[3].1) eq 2;// L[3] is the subgroup of G of order 2 and length 6, and s its non-trivial element
assert ev_c(c_ns[1],s,s) ne ev_c(c_ns[2],s,s);
if c_GL23(sigma,sigma) eq ev_c(c_ns[1],s,s) then
    gamma_plus := c_ns[1];
    gamma_minus := c_ns[2];
else
    gamma_plus := c_ns[2];
    gamma_minus := c_ns[2];
end if;
// gamma_plus corresponds to GL(2,Z/3Z)
// gamma_minus corresponds to B_O
//////////////////////////////////////////////////////////////////////////

// for the rest of the code we fix gamma to be gamma_plus 
gamma := gamma_plus;
// uncomment the following line if you want to try gamma = gamma_plus
// gamma := gamma_minus;

// now we find the conjugacy classes of gamma-regular elements
ccne :=[]; // the conjugacy classes of regular elements for gamma
for g in G_e do
    if is_c_regular(G,gamma,g) then
	class :=  Class(G,g);
        is_already := false;
        for a in ccne do // check if the class is already in the list
	    if a eq class then
                is_already := true; 
            end if;
        end for;
        if is_already eq false then
            Append(~ccne,class); end if;
        end if;
end for;

// We check that cnne contains 3 conjugacy classes; therefore dim(Z(M^c[G])) = 3
print "Number of conjugacy classes of c-regular elements =",#ccne;
// we call x, y, and z three representatives of the conjugacy classes
x := IndexedSet(ccne[1])[1]; y:= IndexedSet(ccne[2])[1]; z:= IndexedSet(ccne[3])[1];
RTx := Transversal(G,Centralizer(G,x)); RTy := Transversal(G,Centralizer(G,y)); RTz := Transversal(G,Centralizer(G,z));
// The above are RIGHT transversals, whereas Tx, Ty, and Tz are LEFT transversals; so we just invert the elements...
Tx := []; Ty := []; Tz := [];
for g in RTx do
	Append(~Tx,g^-1);
end for;
for g in RTy do
	Append(~Ty,g^-1);
end for;
for g in RTz do
	Append(~Tz,g^-1);
end for;

// we compute kx, ky, and kz
// to represent an element of M^gamma[G] we just use vectors of length 24; each position of the vector corresponds to an element of G_e
kx := [];
for i in [1..24] do kx[i] := 0; end for;
for g in Tx do
    j :=find_index(g*x*g^-1, G_e);
    kx[j] := kx[j] + ev_c(gamma,g,g^-1)^-1*ev_c(gamma,g,x)*ev_c(gamma,g*x,g^-1); // this is formula (5.11)
end for;
ky := [];
for i in [1..24] do ky[i] := 0; end for;
for g in Ty do 
	j :=find_index(g*y*g^-1, G_e);
        ky[j] := ky[j] + ev_c(gamma,g,g^-1)^-1*ev_c(gamma,g,y)*ev_c(gamma,g*y,g^-1);
end for;
kz := [];
for i in [1..24] do kz[i] := 0; end for;
for g in Tz do
	j :=find_index(g*z*g^-1, G_e);
        kz[j] := kz[j] + ev_c(gamma,g,g^-1)^-1*ev_c(gamma,g,z)*ev_c(gamma,g*z,g^-1);
end for;
kx := Vector(kx); ky:= Vector(ky); kz := Vector(kz);

// We check that the cocycle gamma restricted to all normal subgroups of G is non-trivial (although a proof of this is given in the proof of Lemma 6.6)
for i in [1..#L] do
    H := L[i];
    if IsNormal(G,H) and #H ne 1 then
	// print "Normal subgroup with i = ",i;
        CMH := Restriction(CM, H);
        res := hom<H2 -> CohomologyGroup(CMH,2) | x:-> IdentifyTwoCocycle(CMH,TwoCocycle(CM,x)) >; 
        //print "The restriction of the cocycle is trivial: ", IsTwoCoboundary(CMH,gamma);
        assert IsTwoCoboundary(CMH,gamma) eq false;
    end if;
end for;

// this implements the multiplication in the twisted group algebra; given two vectors v1, v2 representing elements of M^c[G], it returns the vector corresponding to v1*v2 in M^c[G]
multiply_in_twisted_group_algebra := function(v1,v2,c)
    res := [];
    for i in [1..24] do res[i] := 0; end for;
    for i in [1..24] do
        for j in [1..24] do
	    k := find_index(G_e[i]*G_e[j],G_e);
            res[k] := res[k] + ev_c(c,G_e[i],G_e[j])*v1[i]*v2[j];
	end for;
    end for;
    return Vector(res);
end function;

// compute kx^2 and kx^3      
kx2:=multiply_in_twisted_group_algebra(kx,kx,gamma);
kx3:=multiply_in_twisted_group_algebra(kx2,kx,gamma);
// compute ky^2 and kz^3            
ky2:=multiply_in_twisted_group_algebra(ky,ky,gamma);
ky3:=multiply_in_twisted_group_algebra(ky2,ky,gamma);

// this is the expected relation for gamma_plus
// we compute ky^3 - ky, and it turns out to be zero (that is, the vector of zeros)
print "kx satisfies ky3 - 18*ky";
print  ky3 -18*ky;
// this is the expected relation for gamma_minus
// uncomment the following line if you are testing gamma_minus
// print "kx satisfies ky3 + 18*ky";
// print  ky3 +18*ky;
