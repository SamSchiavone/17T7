AttachSpec("~/projects/CHIMP/CHIMP.spec");
S<X> := PolynomialRing(Integers());
K<nu> := NumberField(X^2 - X - 1);
square := eval square;
sqrt5 := 2*nu - 1;
R<x> := PolynomialRing(K);

function Fbest(t)
  overK := NumberField(390625*t*x^17 + 5859375*t*x^16 + 43125000*t*x^15 + 224375000*t*x^14 + 919687500*t*x^13 + 3018812500*t*x^12 + 8201975000*t*x^11 + 18646725000*t*x^10 + 34614733750*t*x^9 + 52359626250*t*x^8 + 60286491000*t*x^7 + 41159017000*t*x^6 - 1267980500*t*x^5 - 60264765900*t*x^4 - 85488281560*t*x^3 - 54201132040*t*x^2 - 16350628239*t*x - 1908029761*t + 5*x^2 + 27);
  f := DefiningPolynomial(AbsoluteField(overK));
  f := Polredabs(f : Best:=true);
  return f;
end function;
u := sqrt5 - 5;

s := square/u + 3375/8796093022208;
s := square*u + 3375/8796093022208;
try
    f := Fbest(s);
catch e
    f := false;
    print s, "failed";
end try;
if not f cmpeq false then
    Gt, Gn :=  TransitiveGroupIdentification(GaloisGroup(f));
    if [Gt,Gn] ne [82,34] then
        print s, [Gt, Gn];
    end if;
end if;
exit;
