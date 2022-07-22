// Example 4.3.4 in Couvillion's thesis
C4 := CyclicGroup(4);
C5 := CyclicGroup(5);
W := WreathProduct(C5,C4);
D0 := DirectProduct(C5,C4);
D := [el : el in Subgroups(W : OrderEqual := 20) | IsIsomorphic(el`subgroup, D0)][1]`subgroup;
T := Transversal(W,D);
p := RelativeInvariant(W,D);
