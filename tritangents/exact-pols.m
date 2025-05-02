load "sols_refined.m";
#vec;
#vec[1];

// each plane intersects in 3! = 6 points, we only care about the Galois orbit
// which is the hyperplanes, so remove redundancies

h1 := [ComplexField(100) | s[1] : s in vec];
eps := 10^(-80);

uniqs := {[i : i in [1..#h1] | Abs(h1[i]-h1[j]) lt eps] : j in [1..#h1]};
assert &and[#u eq 6 : u in uniqs];
assert #uniqs eq 120;
uniqs := Sort([u[1] : u in uniqs]);

hs := [[vec[i][j] : i in uniqs] : j in [1..3]];

_<xCC> := PolynomialRing(CC);

// rationalize
f := function(x);
  try
    return Roots(PowerRelation(Re(x),1),Rationals())[1][1];
  catch e;
    return 0;
  end try;
end function;

hpols := [];
for i := 1 to 3 do
  hi := &*[xCC-s : s in hs[i]];
  hi := Polynomial([f(c) : c in Coefficients(hi)]);
  Append(~hpols, hi);
  print i;
end for;

hs_refined := [[HenselLift(hpols[i],hs[i][j],100000) : j in [1..#hs[i]]] : i in [1..3]];

// actually, we'll only use one representation (first one) and interpolate the others
// we want to solve beta = sum_i c_i alpha^i with c_i in QQ
// and we know complex approximations v_j(beta) = sum_i c_i v_j(alpha)^i
// rather than linear algebra (inverting Vandermonde matrix), 
// we view this as an interpolation problem f(v_j(alpha)) = v_j(beta)
// so can use Lagrange interpolation

h2interpCC := Interpolation(hs_refined[1],hs_refined[2]);
h2interp := Polynomial([f(c) : c in Coefficients(h2interpCC)]);

h3interpCC := Interpolation(hs_refined[1],hs_refined[3]);
h3interp := Polynomial([f(c) : c in Coefficients(h3interpCC)]);

/*
// exact polynomials
hs := [];
for i := 1 to 3 do
  hi := &*[xCC-s[i] : s in vec];
  hi := Polynomial([f(c) : c in Coefficients(hi)]);
  hifac := Factorization(hi);
  assert #hifac eq 1 and hifac[1][2] eq 6;
  Append(~hs, hifac[1][1]);
  print i;
end for;

xyzs := [];
for i := 4 to 12 do
  hi := &*[xCC-s[i] : s in vec];
  hi := Polynomial([f(c) : c in Coefficients(hi)]);
  hifac := Factorization(hi);
  assert #hifac eq 1 and hifac[1][2] eq 2;
  Append(~xyzs, hifac[1][1]);
  print i;
end for;
*/

// check solutions
K<alpha> := NumberField(hpols[1]);
// K := GF(487);
// alpha := Roots(hpols[1],K)[1][1];
hexact := [alpha,Evaluate(ChangeRing(h2interp,K),alpha),Evaluate(ChangeRing(h3interp,K),alpha)];

R<[x]> := PolynomialRing(K,3);

Igens := [
-8*x[1]^2 + 8*x[1]*x[2] - 34*x[1]*x[3] - 10*x[1] + 17*x[2]^2 - 2*x[2]*x[3] -
          9*x[2] - 28*x[3]^2 - 18*x[3] + 2,
4*x[1]^3 - 6*x[1]^2*x[2] + 12*x[1]^2*x[3] + 2*x[1]^2 - 6*x[1]*x[2]^2 +
          6*x[1]*x[2]*x[3] + 7*x[1]*x[2] - 12*x[1]*x[3]^2 + 4*x[1]*x[3] - 20*x[1] +
          24*x[2]^2*x[3] + 4*x[2]^2 - 13*x[2]*x[3] - 24*x[3]^3 - 8*x[3]^2 - 3*x[3] -
          12,
hexact[1]*x[1] + hexact[2]*x[2] + hexact[3]*x[3] - 1];

r := Resultant(Resultant(Igens[1],Igens[3],x[3]),Resultant(Igens[2],Igens[3],x[3]),x[2]);
assert Degree(r) eq 6;
Factorization(r);

/*

Igens cat:= [Evaluate(ChangeRing(hs[i],k), h[i]) : i in [1..3]];  // first one should be zero
Igens cat:= [Evaluate(ChangeRing(xyzs[i],k), xyz[i]) : i in [1..9]]; 

I := ideal<R | Igens>;
*/

SetVerbose("Groebner", 3);
Dimension(I);
