C4 := CyclicGroup(4);
C5 := CyclicGroup(5);
W := WreathProduct(C5,C4);
D0 := DirectProduct(C5,C4);
D := [el : el in Subgroups(W : OrderEqual := 20) | IsIsomorphic(el`subgroup, D0)][1]`subgroup;
T := Transversal(W,D);
p := RelativeInvariant(W,D);
p;
p0 := p;
p := MultivariatePolynomial(p0);
R<[x]> := Parent(p);
p;
//Evaluate(p, [x[i^W.2] : i in [1..20]]);     
S<y> := PolynomialRing(R);
//&*[y - Evaluate(p, [x[i^t] : i in [1..20]]) : t in T]; // hangs forever
// ParametrizeProjectiveSurface to parametrize a rational surface
