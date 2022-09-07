from lmfdb import db

print("Computing list of HMFs with needed dimension and degree that are not base changes")
hmfs0 = db.hmf_forms.search({'dimension':4, 'deg':2, 'is_base_change':'no'})
hmfs = []
print("Searching for HMFS whose level is Galois invariant, whose Hecke eigenvalue field is totally real and in which 2 is inert")
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
                print("Found one!")
                hmfs.append(rec['label'])
