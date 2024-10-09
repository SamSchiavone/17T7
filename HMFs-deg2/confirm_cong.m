load "2.2.24.1-726.1-i.m";
// "2.2.12.1-578.1-c.m", "2.2.12.1-722.1-i.m", "2.2.8.1-2601.1-j.m", 
// "2.2.8.1-2738.1-e.m", "2.2.12.1-1587.1-i.m", "2.2.24.1-726.1-i.m";

function sigma(aa);
  sigmaF := hom<F -> F | Trace(F.1)-F.1>;
  return ideal<Order(aa) | [sigmaF(a) : a in Generators(aa)]>;
end function;

M := HilbertCuspForms(F, NN);
S := NewSubspace(M);
// SetVerbose("ModFrmHil", 1);

kerf := Kernel(Matrix(ChangeRing(HeckeOperator(S,primes[1]),K)) - heckeEigenvaluesArray[1]);
kerfsigma := Kernel(Matrix(ChangeRing(HeckeOperator(S,sigma(primes[1])),K)) - heckeEigenvaluesArray[1]);

i := 2;
while Dimension(kerf)*Dimension(kerfsigma) gt 1 do
  kerf meet:= Kernel(Matrix(ChangeRing(HeckeOperator(S,primes[i]),K)) - heckeEigenvaluesArray[i]);
  kerfsigma meet:= Kernel(Matrix(ChangeRing(HeckeOperator(S,sigma(primes[i])),K)) - heckeEigenvaluesArray[i]);
  i +:= 1;  
end while;

ZK := Integers(K);
makeintegralodd := function(v);  // make v integral at 2, assuming 2 is inert
  aa := ideal<ZK | Eltseq(v)>;
  return v*2^(-Valuation(aa,2*ZK));
end function;
  
v := makeintegralodd(kerf.1);
vsigma := makeintegralodd(kerfsigma.1);

AutK := Automorphisms(K);
taus := [tau : tau in AutK | tau(K.1) ne K.1 and tau(tau(K.1)) eq K.1];
assert #taus eq 1;
tau := taus[1];
tauv := Vector([tau(vi) : vi in Eltseq(v)]);

// vsigma only defined up to scalar, so line up
j := 1;
while Valuation(v[j],2*ZK) ne 0 do j +:= 1; end while;
assert Valuation(vsigma[j],2*ZK) eq 0;
vsigma *:= tauv[j]/vsigma[j];
assert &and[Valuation(vi,2*ZK) ge 0 : vi in Eltseq(vsigma)];
assert tauv[j]-vsigma[j] eq 0;
k, mk := ResidueClassField(2*ZK);

diffmod2 := [mk(tauv[i]-vsigma[i]) : i in [1..#Eltseq(v)]];
print diffmod2;
if &and[d eq 0 : d in diffmod2] then print "So it holds mod 2!"; end if;
print tauv-vsigma;
if tauv eq vsigma then print "Actually it holds on the nose"; 
else print "but not on the nose"; end if;

