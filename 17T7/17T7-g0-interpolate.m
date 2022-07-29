Attach("polredabs.m");
QQ := Rationals();
K<nu> := QuadraticField(5);
R0<T,X> := PolynomialRing(K,2);
R<T0,T1,xx> := PolynomialRing(K,3);
QQts<t0,t1> := FunctionField(QQ,2);
S0<x> := PolynomialRing(QQts);
QQt<t> := FunctionField(QQ);
S<z> := PolynomialRing(QQt);
// genus 0
// reduced model from Ciaran
F := 390625*T*X^17 + 5859375*T*X^16 + 43125000*T*X^15 + 224375000*T*X^14 + 919687500*T*X^13 + 3018812500*T*X^12 + 8201975000*T*X^11 + 18646725000*T*X^10 + 34614733750*T*X^9 + 52359626250*T*X^8 + 60286491000*T*X^7 + 41159017000*T*X^6 - 1267980500*T*X^5 - 60264765900*T*X^4 - 85488281560*T*X^3 - 54201132040*T*X^2 - 16350628239*T*X - 1908029761*T + 5*X^2 + 27;
Fwr := Evaluate(F, [T0+T1*nu,xx])*Evaluate(F, [T0-T1*nu,xx]);
Fev := Evaluate(Fwr, [T0, 1, xx]);

mp := hom<R -> S0 | [t0, t1, x] >;
mp_cs := hom< QQts -> QQt | [t,0] >;
mp_one_var := hom< S0 -> S | mp_cs, [z]>;

Fev := mp_one_var(mp(Fev));
SetVerbose("GaloisGroup",5);
//G := GaloisGroup(Fev); // hangs; try over finite field

p := 3;
FFu<u> := FunctionField(GF(p));
S_FF<zz> := ChangeRing(S,FFu);
pi := hom< S -> S_FF | S_FF.1 >;
Fev_FF := pi(Fev);
assert IsIrreducible(Fev_FF);
G := GaloisGroup(Fev_FF);
