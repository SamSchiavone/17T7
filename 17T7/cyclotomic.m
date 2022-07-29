function CreateCyclicGroups(prime)
  return CyclicGroup(prime-1), CyclicGroup(prime), prime*(prime-1);
end function;

function WreathProductPolynomial(prime)
  K<z> := CyclotomicField(prime);
  KTs<[T]>:= FunctionField(K, prime-1);
  KXTs<X> := PolynomialRing(KTs);
  g := &*[ (X^prime - &+[KTs.i*z^(j*(i-1)) : i in [1..prime-1]]) : j in [1..prime-1]];
  return g;
end function;

// vals should look like [t,1,1,1], [1,t,1,1], etc.
function SpecializePolynomial(g, vals)
  KXTs<X> := Parent(g);
  KTs<[T]>:= BaseRing(KXTs);
  QQt<t> := FunctionField(QQ);
  QQxt<x> := PolynomialRing(QQt);
  mp_cs := hom< KTs -> QQt | vals >;
  mp := hom< KXTs -> QQxt | mp_cs, [x] >;
  return mp(g);
end function;

function ReducePolynomial(gev, ll)
  QQt<t> := FunctionField(QQ);
  QQxt<x> := PolynomialRing(QQt);
  FFu<u> := FunctionField(GF(ll));
  FFxt<xx> := ChangeRing(QQxt,FFu);
  pi := hom< QQxt -> FFxt | FFxt.1 >;
  return pi(gev);
end function;

function FindGaloisSubfield(gev_FF,prime)
  G := GaloisGroup(gev_FF);
  assert #G eq prime^(prime-1)*(prime-1);
  C1, C2, prod := CreateCyclicGroups(prime);
  P := DirectProduct(C1,C2);
  Hs := [];
  for el in Subgroups(G : OrderEqual := prod) do
    if IsIsomorphic(el`subgroup,P) then
      Append(~Hs, el);
    end if;
  end for;
  assert #Hs eq 1;
  H := Hs[1];
  fH, hH := GaloisSubgroup(gev_FF,H`subgroup); // hangs; try over finite field?
  return fH, hH;
end function;

p := 5;
C1, C2 := CreateCyclicGroups(p);
g := WreathProductPolynomial(p);
KXTs<X> := Parent(g);
KTs<[T]>:= BaseRing(KXTs);
QQt<t> := FunctionField(QQ);
gev := SpecializePolynomial(g, [t,1,1,1]);

polys := [* *];
for ll in PrimesInInterval(3,13) do
  gev_FF := ReducePolynomial(gev,ll);
  fH, _ := FindGaloisSubfield(gev_FF,p);
  Append(~polys, fH);
end for;
