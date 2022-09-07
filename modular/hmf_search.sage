from lmfdb import db

print("Computing list of HMFs with needed dimension and degree that are not base changes")
hmfs0 = db.hmf_forms.search({'dimension':4, 'deg':2, 'is_base_change':'no'})
hmfs = []
print("Searching for HMFS whose level is Galois invariant, whose Hecke eigenvalue field is totally real and in which 2 is inert")
cnt = 0
for rec in hmfs0:
    print("Checking %s" % rec['label'])
    fld_rec = db.nf_fields.lookup(rec['field_label'])
    R.<x> = PolynomialRing(QQ)
    K.<w> = NumberField(R(fld_rec['coeffs']))
    OK = K.ring_of_integers()
    # level ideal == conjugate
    level_gens = eval(rec['level_ideal'])
    level = OK.ideal(level_gens[1:])
    gal = K.automorphisms()
    sigma = gal[1]
    assert sigma(w) != w
    if sigma(level) == level:
        # Hecke stuff
        hecke_rec = db.hmf_hecke.lookup(rec['label'])
        H.<e> = NumberField(R(hecke_rec['hecke_polynomial']))
        if H.is_totally_real():
            OH = H.ring_of_integers()
            if OH.ideal(2).is_prime():
                inv_bool = False
                # check if H has involution
                if H.is_galois():
                    G_H = H.galois_group()
                    if G_H.order() % 2 == 0:
                        inv_bool = True
                else:
                    for phi in H.automorphisms():
                        if (phi(e) != e) and (phi(phi(e)) == e):
                            inv_bool = True
                if inv_bool:
                    print("Found one!")
                    cnt += 1
                    hmfs.append(rec['label'])
                    with open("/scratch/home/sschiavo/github/17T7/modular/hmf_output.txt","a") as my_file:
                        my_file.write("%s\n" % rec["label"])
print("Found %s suitable HMFs" % cnt)
