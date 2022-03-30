AttachSpec("spec");
QQ := Rationals();
L<nu> := QuadraticField(5);
R<t,x> := PolynomialRing(QQ,2);
// genus 0
//F := t*x^17 + 15/4*t*x^16 + 69/10*t*x^15 + 359/40*t*x^14 + 2943/320*t*x^13 + 48301/6400*t*x^12 + 328079/64000*t*x^11 + 745869/256000*t*x^10 + 27691787/20480000*t*x^9 + 41887701/81920000*t*x^8 + 60286491/409600000*t*x^7 + 41159017/1638400000*t*x^6 - 2535961/13107200000*t*x^5 - 602647659/262144000000*t*x^4 - 2137207039/2621440000000*t*x^3 - 1355028301/10485760000000*t*x^2 - 16350628239/1677721600000000*t*x - 1908029761/6710886400000000*t + 8192/263671875*x^2 + 512/48828125;
// reduced model from Ciaran
F := 390625*t*x^17 + 5859375*t*x^16 + 43125000*t*x^15 + 224375000*t*x^14 + 919687500*t*x^13 + 3018812500*t*x^12 + 8201975000*t*x^11 + 18646725000*t*x^10 + 34614733750*t*x^9 + 52359626250*t*x^8 + 60286491000*t*x^7 + 41159017000*t*x^6 - 1267980500*t*x^5 - 60264765900*t*x^4 - 85488281560*t*x^3 - 54201132040*t*x^2 - 16350628239*t*x - 1908029761*t + 5*x^2 + 27;
S<s> := PolynomialRing(L);
t0s := [a+b*nu : a,b in [-10..10]];
//t0s := [t0 : t0 in t0s | not t0 in [0,1]];
t0s := [t0 : t0 in t0s | not t0 in QQ];
t0s := SetToSequence(SequenceToSet(t0s));
special_t0s := [];
//t0 := nu + 4;
for t0 in t0s do
  printf "evaluating at %o\n", t0;
  Fev := Evaluate(F, [t0, s]);
  Grel := GaloisGroup(Fev);
  assert TransitiveGroupIdentification(Grel) eq 7;
  print "relative extension has Galois group 17T7";
  M := ext< L | Fev>;
  Mabs := AbsoluteField(M);
  Mabs_poly := DefiningPolynomial(Mabs);
  Mabs_poly := Polredbest(Mabs_poly);
  printf "polredbest-ed polynomial for extension over QQ:\n%o\n", Mabs_poly;
  G, seq, data := GaloisGroup(Mabs_poly);
  k := TransitiveGroupIdentification(G);
  if k ne 88 then
    print "\n\nsmaller than generic Galois group found!\n\n";
    print t0;
    print G;
    Append(~special_t0s, [* t0, k *]);
  end if;
end for;
