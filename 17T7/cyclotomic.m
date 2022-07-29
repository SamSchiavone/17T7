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
  gev_FF := pi(gev);
  return gev_FF;
end function;

function FindGaloisSubfield(gev_FF,prime)
  G, roots, data := GaloisGroup(gev_FF);
  //printf "G is %o\n", G;
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
  //printf "H is %o\n", H;
  fH, hH := GaloisSubgroup(data,H`subgroup);
  return fH, hH;
end function;

//----------------------------------------------------------------------------------------
// now compute!
//----------------------------------------------------------------------------------------

SetVerbose("GaloisGroup",0);
p := 5;
C1, C2 := CreateCyclicGroups(p);
g := WreathProductPolynomial(p);
KXTs<X> := Parent(g);
KTs<[T]>:= BaseRing(KXTs);
QQt<t> := FunctionField(QQ);
gev := SpecializePolynomial(g, [t,1,1,1]);

polys := [* *];
primes := PrimesInInterval(5,101);
for ll in primes do
  printf "prime for reduction ll = %o\n", ll;
  gev_FF := ReducePolynomial(gev,ll);
  if not IsIrreducible(gev_FF) then
    print "reduced polynomial is reducible";
    continue;
  end if;
  fH, _ := FindGaloisSubfield(gev_FF,p);
  Append(~polys, fH);
end for;

// now CRT
moduli := [Characteristic(Parent(poly)) : poly in polys];
cs := [* Coefficients(el) : el in polys *];
for i := 1 to #cs[1] do
  for el in cs do // cs[j] is a list of polys in u
    //el[i]; // el[i] is a polynomial in u
  end for;
end for;
