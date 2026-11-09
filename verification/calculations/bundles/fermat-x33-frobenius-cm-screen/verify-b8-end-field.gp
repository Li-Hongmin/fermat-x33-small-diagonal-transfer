\\ Exact characteristic-zero endomorphism-field certificate for the fixed B8.
\\ Run with PARI/GP 2.17 or later:
\\   gp -q calculations/bundles/fermat-x33-frobenius-cm-screen/verify-b8-end-field.gp

x = 'x;
y = 'y;
z = 'z;

P331 = x^16 + x^15 + 761*x^14 - 13183*x^13 + 147512*x^12 - 9668365*x^11 + 99077160*x^10 - 1981264536*x^9 + 64505725132*x^8 - 655798561416*x^7 + 10854992726760*x^6 - 350620269200215*x^5 + 1770676919700152*x^4 - 52378640371812133*x^3 + 1000812265940691041*x^2 + 435307306210734211*x + 144086718355753023841;

P661 = x^16 + 6*x^15 + 454*x^14 + 9182*x^13 + 716346*x^12 + 23307588*x^11 + 314636724*x^10 + 16817734524*x^9 + 266633006035*x^8 + 11116522520364*x^7 + 137471392086804*x^6 + 6731342847978228*x^5 + 136750422918799386*x^4 + 1158629510490621782*x^3 + 37867323493919914294*x^2 + 330796927261864273326*x + 36442794820015380778081;

Epol = z^4 + 28*z^2 + 64;

root_ratio_data(P) =
{
  my(R, fac, degrees, cyclo);
  R = polresultant(P, subst(P, x, y*x), x);
  fac = factor(R);
  degrees = List();
  cyclo = List();
  for(i = 1, matsize(fac)[1],
    if(type(fac[i,1]) != "t_POL", next);
    for(j = 1, fac[i,2], listput(degrees, poldegree(fac[i,1])));
    my(n = poliscyclo(fac[i,1]));
    if(n, listput(cyclo, n));
  );
  [Vec(degrees), Vec(cyclo)];
}

check(condition, message) = if(!condition, error(message));

subfield_data(P) =
{
  my(S, counts, quartics, octics);
  S = nfsubfields(P);
  counts = vector(17, i, 0);
  quartics = List();
  octics = List();
  for(i = 1, #S,
    my(d = poldegree(S[i][1]));
    counts[d+1]++;
    if(d == 4, listput(quartics, S[i][1]));
    if(d == 8, listput(octics, S[i][1]));
  );
  [counts, Vec(quartics), Vec(octics)];
}

verify_prime(q, P) =
{
  my(rr, sd, Q4, Q8);
  print("q=", q);
  print("irreducible=", polisirreducible(P));
  rr = root_ratio_data(P);
  print("root-ratio factor degrees=", vecsort(rr[1]));
  print("cyclotomic factor orders=", vecsort(rr[2]));
  sd = subfield_data(P);
  print("subfield counts by degree 1..16=", sd[1][2..17]);
  Q4 = sd[2];
  Q8 = sd[3];
  print("quartic count=", #Q4);
  print("quartic field discriminant=", nfdisc(Q4[1]));
  print("quartic isomorphic to E=", nfisisom(Q4[1], Epol) != 0);
  print("octic count=", #Q8);
  print("octic field discriminant=", nfdisc(Q8[1]));
  print("full field discriminant=", nfdisc(P));
  print();
  [Q4[1], Q8[1], nfdisc(Q8[1]), nfdisc(P)];
}

check(polisirreducible(Epol), "E polynomial is reducible");
check(nfdisc(Epol) == 1089, "unexpected E discriminant");

D331 = verify_prime(331, P331);
D661 = verify_prime(661, P661);

check(D331[3] != D661[3], "octic discriminants should differ");
check(D331[4] != D661[4], "degree-16 discriminants should differ");
print("common E-compatible subfield degrees 8 and 16 excluded by discriminants");
print("status: End^0_Qbar(B8) = E, using good-reduction specialization");

nfE = nfinit(Epol);
zE = Mod(z, Epol);
e1 = (zE^3 + 20*zE - 16)/32;
e2 = -(zE^3 + 20*zE + 16)/32;
e3 = -e1;
e4 = -e2;

check(e1^2 + e1 + 1 == 0, "q=331 first normalized determinant is not a cube root");
check(e2^2 + e2 + 1 == 0, "q=331 second normalized determinant is not a cube root");
check(e3^2 - e3 + 1 == 0, "q=661 first normalized determinant is not a sixth root");
check(e4^2 - e4 + 1 == 0, "q=661 second normalized determinant is not a sixth root");

determinant_constants(q, P) =
{
  my(F);
  F = nffactor(nfE, P);
  check(matsize(F)[1] == 4, "Frobenius polynomial does not split into four factors over E");
  for(i = 1, 4,
    check(poldegree(F[i,1]) == 4 && F[i,2] == 1,
      "unexpected E-factor degree or multiplicity")
  );
  vector(4, i, polcoeff(F[i,1], 0)/q^2)
}

count_equal(v, c) =
{
  my(n = 0);
  for(i = 1, #v, if(v[i] == c, n++));
  n;
}

C331 = determinant_constants(331, P331);
C661 = determinant_constants(661, P661);
check(count_equal(C331, e1) == 2, "q=331 first determinant multiplicity is not two");
check(count_equal(C331, e2) == 2, "q=331 second determinant multiplicity is not two");
check(count_equal(C661, e3) == 2, "q=661 first determinant multiplicity is not two");
check(count_equal(C661, e4) == 2, "q=661 second determinant multiplicity is not two");
check(lift(Mod(5, 331)^110) == 31, "unexpected cubic-residue test at q=331");
check(lift(Mod(5, 661)^220) == 364, "unexpected cubic-residue test at q=661");
print("q=331: Frobenius on W_E(B8)(2) has characteristic polynomial (T^2+T+1)^2");
print("q=661: Frobenius on W_E(B8)(2) has characteristic polynomial (T^2-T+1)^2");
print("minimal residue degrees fixing W_E(B8)(2): 3 at q=331; 6 at q=661");
