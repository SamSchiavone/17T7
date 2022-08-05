function IsRational(C)
  QQ := Rationals();
  G := Parent(C[3]);
  Cs := ConjugacyClasses(G);
  i := Index(Cs,C);
  tab := CharacterTable(G);
  for chi in tab do
    if not chi[i] in QQ then
      return false;
    end if;
  end for;
  return true;
end function;

G := TransitiveGroup(17,7);
Cs := ConjugacyClasses(G);
tab := CharacterTable(G);
//Cs_rat := Cs[1..5] cat [Cs[8]];
Cs_rat := [C : C in Cs | IsRational(C)];
inds := [Index(Cs, el) : el in Cs_rat];
H := sub< G | [el[3] : el in Cs_rat] >;
G eq H;
// mass formula: see Serre, Topics in Galois Theory, Theorem 7.2.1 (p. 68)
n := &+[(&*[chi[ind] : ind in inds])/(chi[1]^(#inds-2)) : chi in tab]; // n = number of solutions
n *:= (1/#G)*&*[C[2] : C in Cs_rat];
n;
