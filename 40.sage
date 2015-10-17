###################################
# Some utils
##################################

# computes C^sigma, where C is a curve and sigma an automorphism of the base field of C
def conjugate_curve(C,sigma):
    new_as = [sigma(a) for a in C.ainvs()]
    return EllipticCurve(new_as)

# given a polynomial P over a field and an automorphism sigma of the field, computes P^sigma
def conjugate_pol(P,sigma):
    P_sigma = parent(P)(0)
    x = parent(P).gen()
    i = 0
    for c in P.coefficients():
        P_sigma = P_sigma + sigma(c)*x**i
        i += 1
    return P_sigma

####################################
# The main computation
####################################
D = -40
R.<x> = PolynomialRing(QQ)
hcp = hilbert_class_polynomial(D)
M.<gM> = NumberField((x^2-D))
W.<y> = PolynomialRing(M)
H.<g> = NumberField(W(hcp))
H.<g> = H.absolute_field() # the hilbert class field of M
G = H.automorphisms()
W.<x> = PolynomialRing(H)
hcpH = W(hcp) # the hilbert class polynomial, viewed as having coefficients in H
# roots of the hilbert class polynomial, viewed as elements of H
j0= - hcpH.factor()[0][0][0] 
j1= - hcpH.factor()[1][0][0]
# these are two elliptic curves with CM by M
E0 = EllipticCurve_from_j(j0)
E1 = EllipticCurve_from_j(j1)

# Now we find sigma. It is the non-trivial automorphism of H that fixes M;
rD = (x^2 - D).roots()[0][0] # the square root of D, as an element of H
r = (x^2 - hcp.discriminant().squarefree_part()).roots()[0][0] # r generates Q(j)
for s in G:
    if s(rD) == rD and s(r) != r: # sigma should fix the square root of D and it does not fix Q(j)
        sigma = s
        break
# Now we twist E0 by the following d (that we found by a simple search method)
d = -7/4311049539847060068249600*g^3 + 34489/33264271140795216576*g^2 - 4779151372807/107776238496176501706240*g - 502768059522914115499/12474101677798206216
E = E0.quadratic_twist(d)
E_sigma = conjugate_curve(E,sigma)
assert E.is_isogenous(E_sigma) # we check that E is indeed isogenous to its Gal(H/F)-conjugate
# This is an isogeny of degree 2 from E
mu = E.isogenies_prime_degree(2)[0]
# we check that the codomain of the isogeny is isomorphic to E_sigma
assert mu.codomain().is_isomorphic(E_sigma)

# Now we compute the conjugate of mu 
E2 = mu.codomain()
deg = mu.degree()
# sometimes E2 is not equal to E_sigma, but isomorphic; we compose mu with this iso
phi = E2.isomorphism_to(E_sigma)
mu.set_post_isomorphism(phi)
assert mu.codomain().ainvs() == E_sigma.ainvs() # mu is indeed an isogeny between E and E_sigma

# the isogeny is given by an expression of the form
# mu(x,y) = (u(x)/v(x), s(x)/t(x)*y) 
# see, for instance http://math.mit.edu/classes/18.783/LectureNotes6.pdf
# therefore, the x-coordinate of the isogeny is u/v
u = mu[0].numerator()
v = mu[0].denominator()
# mu, up to multiplication by +-1, is determined by the kernel polynomial psi, such that v =  psi**2 (see Remark 6.14 of loc.cit.)
psi = mu.kernel_polynomial()
# we conjugate the coefficients of psi
psi_sigma = conjugate_pol(psi,sigma)

# mu_sigma is the isogeny with kernel polynomial the conjugate of psi
mu_sigma = E_sigma.isogeny(W(psi_sigma).monic())
# once again, sometimes the codomain of mu_sigma is just isomorphic to E (and not equal)
phi_prime =  mu_sigma.codomain().isomorphism_to(E)
mu_sigma.set_post_isomorphism(phi_prime)

# so now we have mu and mu_sigma; let us check a few things
assert mu.domain() == E and mu.codomain() == E_sigma # mu: E --> E_sigma
# observe that here mu has source E and target E_sigma; in the paper we do it the other way around, but there's no problem at all (the mu_sigma in the code is, therefore, what in the paper is called mu_sigma^sigma) 
assert mu_sigma.domain() == E_sigma and mu_sigma.codomain() == E # mu_sigma: E_sigma --> E
# the kernel polynomials are one the conjugate of the other
# this means that mu_sigma equals the conjugate of mu by sigma, perhaps up to multiplication by -1
assert conjugate_pol(mu_sigma.kernel_polynomial(),sigma) == mu.kernel_polynomial()

# now we compute the composition
u_sigma = mu_sigma[0].numerator()
v_sigma = mu_sigma[0].denominator()
psi_mumusigma = W((u_sigma.substitute(x=u/v)/v_sigma.substitute(x=u/v)).denominator()) # this is the denominator of the x-coordinate of the composition (so the kernel polynomial)
psi_mumusigma = psi_mumusigma.monic() # we make it monic because we want to construct the isogeny that has this polynomial as kernel polynomial (and it needs to be monic)
mumusigma = E.isogeny(psi_mumusigma)
# we check if mumusigma coincides with multiplication by the degree
ker_mult_by_m = E.multiplication_by_m_isogeny(deg).kernel_polynomial()
if mumusigma.kernel_polynomial() == ker_mult_by_m:
    print "mumusigma EQUALS  ",deg," (up to +-1)"
else:
    print "mumusigma NOT  ",deg

