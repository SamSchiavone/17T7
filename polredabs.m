/* pari interpreter */

intrinsic Polred(f::RngUPolElt) -> RngUPolElt
  { A smallest generating polynomial of the number field, using pari. }
  cmd := Sprintf(
    "{print(Vecrev(Vec(polredabs(Pol(Vecrev(%o))))))}", Coefficients(f));
  s := Pipe("gp -q --default parisizemax=1G", cmd);
  ss := [ StringToInteger(x) : x in Split(s, ", []\n") | x ne "" ];
  return Parent(f) ! ss;
end intrinsic;

intrinsic Polredbest(f::RngUPolElt) -> RngUPolElt
  { A smallest generating polynomial of the number field, using pari. }
  cmd := Sprintf(
   "{print(Vecrev(Vec(polredbest(Pol(Vecrev(%o))))))}",
   Coefficients(f));
  s := Pipe("gp -q --default parisizemax=1G", cmd);
  ss := [ StringToInteger(x)
  : x in Split(s, ", []\n") | x ne "" ];
  return Parent(f) ! ss;
end intrinsic;

intrinsic Polredbestabs(f::RngUPolElt) -> RngUPolElt, SeqEnum, BoolElt
  { A smallest generating polynomial of the number field, using pari.  First polredbest, then polredabs.}

  fbest, fbest_root := Polredabs(f : Best := true);
  fredabs, fredabs_root, bl := Polredabs(fbest);

  K := NumberField(f);
  Kbest := NumberField(fbest);
  iotabest := hom<K -> Kbest | fbest_root>;
  Kredabs := NumberField(fredabs);
  iotaredabs := hom<Kbest -> Kredabs | fredabs_root>;
  iota := iotabest*iotaredabs;  // functional composition is backwards in Magma, for some reason
  return fredabs, Eltseq(iota(K.1)), bl;
end intrinsic;

intrinsic Polredabs_db(f::RngUPolElt) -> RngUPolElt
  {first do best...twice then do abs.}
  return Polred(Polredbest(Polredbest(f)));
end intrinsic;

/* polredabsify base fields */

intrinsic Polredabs_db(K::FldNum) -> FldNum
  {}
  if Degree(K) le 1 then
    assert DefiningPolynomial(K) eq DefiningPolynomial(RationalsAsNumberField());
    return K;
  else
    return NumberField(Polredabs_db(DefiningPolynomial(K)));
  end if;
end intrinsic;
