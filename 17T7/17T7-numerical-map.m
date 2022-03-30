// intrinsic BelyiMap(Gamma::GrpPSL2Tri : prec := 0, Al := "Default", ExactAl := "AlgebraicNumbers", DegreeBound := 0, precNewto    n := 0, Federalize := true) -> Any, Any, Any

k := 0;
while SkDimension(Gamma,k) lt 2 do
  k +:= 2;
end while;
Sk := PowerSeriesBasis(Gamma, k : dim := 2, Federalize := Federalize);
if precNewton eq 0 then
  precNewton := 200*m0*Round(Log(d));
end if;
//_ := TrianglePhiGenusZeroNumericalBelyiMap(Sk, Gamma : Al := Al, m := DegreeBound, precNewton := precNewton);     

ram := TriangleRamificationValues(Sk, Gamma : NormalizeByTheta := false);
fact := TriangleGenusZeroFactorizationPattern(Gamma);

polsvec := [];
CC := Parent(ram[1][1][1]);
_<xCC> := PolynomialRing(CC);

if Al eq "ByRamification" then
  Theta := TriangleTheta(Sk, Gamma);
else
  Theta := 1;
end if;
