// SachiRaymond.ipynb
prime := 5;
C1 := CyclicGroup(prime-1);
C2 := CyclicGroup(prime);
prod := prime*(prime-1);
/*
  //rt_of_unity := 3;
  //assert( GCD(rt_of_unity, prime) eq 1 and Order(Integers(prime)!rt_of_unity) eq prime - 1);
  W := WreathProduct(C3,C2);
  D0 := DirectProduct(C3,C2);
  D := [el : el in Subgroups(W : OrderEqual := (#C3*#C2)) | IsIsomorphic(el`subgroup, D0)][1]`subgroup;
  T := Transversal(W,D);
  p := RelativeInvariant(W,D);
  R := Parent(p);
  S<y> := PolynomialRing(R);
  h := &*[y - Evaluate(p, [R.(i^t) : i in [1..Degree(W)]]) : t in T];
  cs := Coefficients(h);
  invs := [];
  for pi in Sym((prime-1)*prime) do
    inv_bool := true;
    for c in cs do 
      if not IsInvariant(c, pi) then
        inv_bool := false;
      end if;
    end for;
    if inv_bool then
      Append(~invs, pi);
    end if;
  end for;

  // group under which coefficients are invariant
  //K := sub< Sym((prime-1)*prime) | invs>; 
  //ReduceGenerators(~K);

  // compute G-relative H-invariant resolvent
  //pp := MultivariatePolynomial(p);
  R0 := Parent(p);
  //R1 := BaseChange(R0,QQ);
  R3 := PolynomialRing(QQ,prime^2);
  y := R3.(prime^2);
  SS<yy> := PolynomialRing(R0);
  print(SS);
  hh := &*[yy - Evaluate(p, [R0.(i^t) : i in [1..Degree(W)]]) : t in T];
  print(1);
  S := PolynomialRing(QQ,six);
  mp_S := hom< S -> R3 | [R3.i : i in [1..six]] >;
  syms := [mp_S(ElementarySymmetricPolynomial(S,i)) : i in [1..six]];
  //cs_hh := Coefficients(hh);
*/

K<z> := CyclotomicField(prime);
KTs<[T]>:= FunctionField(K, prime-1);
KXTs<X> := PolynomialRing(KTs);
g := &*[ (X^prime - &+[KTs.i*z^(j*(i-1)) : i in [1..prime-1]]) : j in [1..prime-1]];
print(g);
QQt<t> := FunctionField(QQ);
QQxt<x> := PolynomialRing(QQt);
mp_cs := hom< KTs -> QQt | [t,1,1,1] >;
mp := hom< KXTs -> QQxt | mp_cs, [x] >;
gev := mp(g);
SetVerbose("GaloisGroup",5);
G := GaloisGroup(gev);
assert #G eq prime^(prime-1)*(prime-1);
P := DirectProduct(C1,C2);
Hs := [];
for el in Subgroups(G : OrderEqual := prod) do
  if IsIsomorphic(el`subgroup,P) then
    Append(~Hs, el);
  end if;
end for;
assert #Hs eq 1;
H := Hs[1];
GaloisSubgroup(gev,H`subgroup); // hangs; try over finite field?

gens := [syms[six+1-i]-temp(elt) : i->elt in Coefficients(g)[1..six]] where temp := hom<KTs -> R3 | [R3.(six+i) : i in [1..prime-1]]>;
mp_x := hom<BaseRing(Parent(hh)) -> R3 | [R3.i : i in [1..six]]>;
mp2 := hom<Parent(hh) -> R3 | mp_x, y>;
gens cat:= [mp2(hh)];

