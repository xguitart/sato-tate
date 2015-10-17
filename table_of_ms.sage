####################################################
# Script for checking Table 2
###################################################
# Choose the discriminant to be checked
D = -24 # D can be one of -24, -35, -40, -51, -115 
print "Now doing D = ",D
R.<x> = PolynomialRing(QQ)
hcp = hilbert_class_polynomial(D)
K.<r> = NumberField((x^2-D))
W.<y> = PolynomialRing(K)
H.<g> = NumberField(W(hcp))
H.<g> = H.absolute_field()
W.<x> = PolynomialRing(H)
hcpH = W(hcp)
j0= - hcpH.factor()[0][0][0]
j1= - hcpH.factor()[1][0][0]
E0 = EllipticCurve_from_j(j0)
E1 = EllipticCurve_from_j(j1)
G = H.automorphisms()
# We find sigma. It is the non-trivial automorphism of H that fixes K;
rD = (x^2 - D).roots()[0][0] # this is the square root of D
r = (x^2 - hcp.discriminant().squarefree_part()).roots()[0][0] # r generates the ring of integers of Q(j)
for s in G:
    if s(rD) == rD and s(r) != r: # sigma fixes the square root of -D and it does not fix Q(j)
        sigma = s
	break

def conjugate(C,sigma): # this computes C^sigma, where C is a curve
    new_as = [sigma(a) for a in C.ainvs()]
    return EllipticCurve(new_as)

# the next step is to twist E0 by some d in Q(j), in such a way that the resulting curve is isogenous to it Gal(H/F)-conjugate
# for each D we precomputed a d that works, which we now select
ds = [-5/19904993473440411648*g^3 + 6995/3839697815092672*g^2 - 482570375/43651301476843008*g - 13473044164908229/479962226886584,  -1/16685712674460724428800*g^3 - 108/10184150802283157*g^2 + 80530636793/3337142534892144885760*g + 209760385008210112/10184150802283157, -7/4311049539847060068249600*g^3 + 34489/33264271140795216576*g^2 - 4779151372807/107776238496176501706240*g - 502768059522914115499/12474101677798206216, -5/10315810342865769980290891776*g^3 - 31315/7773174026199732636094*g^2 - 31310311587925/3438603447621923326763630592*g + 105835869598207865674827/7773174026199732636094, -1/8757378258340164562358992839332497829068800*g^3 - 4836071/65988635844215408606702962535961747*g^2 - 78138796356403223/1751475651668032912471798567866499565813760*g - 3516652924182489887950815247747708985/197965907532646225820108887607885241 ]
if D == -24:
    d = ds[0]
elif D == -35:
    d = ds[1]
elif D == -40:
    d = ds[2]
elif D == -51:
    d = ds[3]
elif D == -115:
    d =  ds[4]
E = E0.quadratic_twist(d)
E_sigma = conjugate(E,sigma)
assert E_sigma.is_isogenous(E) # indeed, it is isogenous to it sGal(H/F)-conjugate

# we look for isogenies with source E, and check whether one of those has target (isomorphic to) E_sigma
Isog = E.isogenies_prime_degree()
is_work = [] # this will contain the found isogenies 
for iso in Isog:
    Ei = iso.codomain()
    if Ei.is_isomorphic(E_sigma) and Ei.is_isogenous(E):
        is_work.append(iso)
        print "success! isogeny found with degree ",iso.degree()
# for every isogeny mu: E --> E_sigma, say of degree m, we compute sigma(mu)·mu: E --> E and check whether this composition equals multiplication by m 
for mu in is_work:     
    E2 = mu.codomain()
    deg = mu.degree()
    print 'Now doing isogeny of degree ',deg
    # sometimes E2 is not equal to E_sigma, but isomorphic; we compose mu with this iso
    phi = E2.isomorphism_to(E_sigma)
    mu.set_post_isomorphism(phi)
    # the x-coordinate of the isogeny is u/v
    u = mu[0].numerator()
    v = mu[0].denominator()
    # mu is determined by the kernel polynomial
    psi = mu.kernel_polynomial()
    # we conjugate the coefficients of u
    psi_sigma = W(0)
    i = 0;
    for c in psi.coefficients():
        psi_sigma = psi_sigma + sigma(c)*x**i
        i += 1
    # mu_sigma is the isogeny with kernel polynomial the conjugate of psi
    mu_sigma = E_sigma.isogeny(W(psi_sigma).monic())
    # once again, sometimes the codomain of mu_sigma is just isomorphic to E (and not equal)
    phi_prime =  mu_sigma.codomain().isomorphism_to(E)
    mu_sigma.set_post_isomorphism(phi_prime)
    u_sigma = mu_sigma[0].numerator()
    v_sigma = mu_sigma[0].denominator()
    # this is the x-coordinate of the composition
    v_mumusigma = W((u_sigma.substitute(x=u/v)/v_sigma.substitute(x=u/v)).denominator()).monic()
    # so the kernel polynomial is just the square root of the x-coordinate
    psi_mumusigma = v_mumusigma.squarefree_decomposition()[0][0]
    mumusigma = E.isogeny(psi_mumusigma)
    # sometimes the codomain of mumusigma is not equal to E
    phi_prime_prime =  mumusigma.codomain().isomorphism_to(E)
    mumusigma.set_post_isomorphism(phi_prime_prime)
    # we check if mumusigma coincides with multiplication by the degree
    ker_mult_by_m = E.multiplication_by_m_isogeny(deg).kernel_polynomial()
    if mumusigma.kernel_polynomial() == ker_mult_by_m:
        print "mumusigma EQUALS  ",deg
        break
    else:
        print "mumusigma NOT  ",deg
