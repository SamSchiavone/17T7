from lmfdb import db

print("Computing list of HMFs with needed dimension and degree")
hmfs0 = db.hmf_forms.search({'dimension':4, 'deg':2})
hmfs = []
print("Searching for HMFS whose level is Galois invariant")
for rec in hmfs0:
    print("Checking %s" % rec['label'])
    fld_rec = db.nf_fields.lookup(rec['field_label'])
    R.<x> = PolynomialRing(QQ)
    K.<w> = NumberField(R(fld_rec['coeffs']))
    OK = K.ring_of_integers()
    level_gens = eval(rec['level_ideal'])
    level = OK.ideal(level_gens)
    gal = K.automorphisms()
    sigma = gal[1]
    assert sigma(w) != w
    # level ideal == conjugate
    if sigma(level) == level:
        print("Galois invariant level found!")
        hmfs.append(rec['label'])
