import sys
#sys.path.append(os.path.join(os.environ["HOME"], "lmfdb"))
sys.path.append("lmfdb/psycodict")
from psycodict import psycopg2
sys.path.append("lmfdb")
from lmfdb import db

@cached_function
def lmfdb_field(label):
    R = PolynomialRing(QQ, 'x')
    return NumberField(R(db.nf_fields.lookup(label, "coeffs")), "w")

@cached_function
def number_field_from_string(input):
    R = PolynomialRing(QQ, 'x')
    return NumberField(R(input), "e")

def primes_list(K, base_field_label):
	a = K.gen()
	OK = K.ring_of_integers()
	prime_strings = db.hmf_fields.lucky({'label' : base_field_label}, 'primes')
	return [OK.fractional_ideal(sage_eval(s, {"w" : a})) for s in prime_strings]

def ap_list(K, hecke_eigenvalues, primes):
	e = K.gen()
	L = dict()
	for i in [0..len(primes)-1]:
		L[primes[i]] = K(sage_eval(str(hecke_eigenvalues[i]), {'e' : e}))
	return L

def has_surjective_trace(level_norm, hecke_eigenvalues, prime=2, proof=False):
	assert not(proof)
	# primes = hecke_eigenvalues.keys()
	hecke_field = list(hecke_eigenvalues.values())[0].parent()
	bad_prime_product = prime*level_norm
	good_primes = [ p for p in hecke_eigenvalues.keys() if not((p.norm().radical()).divides(bad_prime_product)) ]
	#frp = hecke_field.fractional_ideal(prime)
	#Fq = frp.residue_field()
	OH = hecke_field.ring_of_integers()
	frp = OH.fractional_ideal(prime)
	Fq = OH.quo(frp, 'zbar') # prime is inert in H
	reduced_eigenvals = [Fq(hecke_field(hecke_eigenvalues[p])) for p in good_primes]
	#histogram = [reduced_eigenvals.count(x) > 0 for x in Fq]
	#return histogram == [True for x in Fq]
	return len(set(reduced_eigenvals)) == prime**hecke_field.degree()

def get_elements_of_order(H, d):
	e = H.gens()[0]
    #auts = H.galois_group()
	auts = H.automorphisms()
	assert d == 2
	return [f for f in auts if (not f.is_identity()) and (f(f(e)) == e)]
    #return [f for f in auts if f.order() == d]

def check_congruence_mod_frp(F, frp, iota, p, ap_dict, field_degree=2, prime=2):
	assert F.degree() == field_degree
	#sigma = [f for f in F.automorphisms() if not f.is_identity()][0]
	sigma = get_elements_of_order(F, field_degree)[0]
	ap = ap_dict[p]
	apbar = ap_dict[sigma(p)]
	return apbar - iota(ap) in frp


def search(ab_var_dimension=4, field_degree=2, prime=2, field_labels=None):
    res = []
    res_failed = []

    R = PolynomialRing(QQ, 'x')
    query  = {'dimension':ab_var_dimension, 'deg':field_degree, 'is_base_change':'no'}
    log_file = open(f"output_SL2_F_{prime**ab_var_dimension}.txt", "w")
    if field_labels is None:
        field_labels = db.hmf_forms.distinct('field_label', query)
        # print(field_labels)

    for field_label in field_labels:
        F = lmfdb_field(field_label)
        assert F.degree() == field_degree
        w = F.gen()
        OK = F.ring_of_integers()
        gal = F.galois_group()
        assert gal.is_cyclic()
        sigma = gal[1]
        assert sigma(w) != w
        query['field_label'] = field_label
        records = list(db.hmf_forms.search(query, ['field_label', 'level_ideal', 'label', 'level_norm']))
        levels = set(elt['level_ideal'] for elt in records)
        ps = primes_list(F, field_label)
        invariant_levels = [] 
        for elt in levels:
            # level ideal == conjugate
            level_gens = eval(elt)
            level = OK.fractional_ideal(level_gens[1:])
            if sigma(level) == level:
                invariant_levels.append(elt)
        invariant_levels = set(invariant_levels)
        level_records = [elt for elt in records if elt['level_ideal'] in invariant_levels]

        labels = [elt['label'] for elt in level_records]
        #hecke_field_dict = dict((rec["label"], rec["hecke_polynomial"]) for rec in db.hmf_hecke.search({"label": {"$like": f'{field_label}-%'}}, ["label", "hecke_polynomial"]))
        hecke_field_dict = dict((rec["label"], [rec["hecke_polynomial"], rec["hecke_eigenvalues"]] ) for rec in db.hmf_hecke.search({"label": {"$in": labels}}, ["label", "hecke_polynomial", "hecke_eigenvalues"]))
        print(field_label, len(level_records))
        for rec in level_records:
            # Hecke stuff
            # hecke_rec = db.hmf_hecke.lookup(rec['label'], ["hecke_polynomial"])
            # H.<e> = NumberField(R(hecke_rec["hecke_polynomial"]))
            H = number_field_from_string(hecke_field_dict[rec['label']][0])
            if not H.is_totally_real():
                continue

            OH = H.ring_of_integers()
            frp = OH.fractional_ideal(prime)
            if not frp.is_prime():
                # res_failed.append(("p is not inert in Hecke field", rec))
                continue

            # check if Gal(F|Q) injects into Gal(H|Q)
            if H.is_galois() and H.galois_group().order() % field_degree != 0:
                # res_failed.append(("no involution by Galois", rec))
                continue

            sigmas = get_elements_of_order(H, field_degree)
            if not sigmas:
                # res_failed.append(("no elements of order {field_degree}", rec))
                continue

            aps = ap_list(H, hecke_field_dict[rec['label']][1], ps)
            if not has_surjective_trace(rec['level_norm'], aps, prime=prime):
                res_failed.append(("not surjective trace", rec))
                continue

            for sigma in sigmas:
                for p in ps:
                    if not check_congruence_mod_frp(F, frp, sigma, p, aps, prime=prime, field_degree=field_degree):
                        break
                else:
                    log_file.write(rec['label'] + "\n")
                    log_file.flush()
                    res.append(rec)
                    break
            else:
                res_failed.append(("not mod {prime}", rec))
    log_file.close()
    return [res, res_failed]

foo2 = search(field_labels=["2.2.8.1", "2.2.12.1", "2.2.5.1", "2.2.24.1"])
len(foo2[0])