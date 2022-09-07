from lmfdb import db

R.<x> = QQ[]
diff_list = []
with open("hmf_output.txt", "r") as hmf_output:
	for label in hmf_output:
		hecke_entry = db.hmf_hecke.lucky({'label' : label[:-1]})
		hecke_eigenvalues = hecke_entry['hecke_eigenvalues']
		f = sage_eval(hecke_entry['hecke_polynomial'], {'x' : x})
		K.<e> = NumberField(f)
		F16 = K.residue_field(2)
		reduced_eigenvals = [F16(sage_eval(str(x), {'e' : e})) for x in hecke_eigenvalues[1:]]
		histogram = [reduced_eigenvals.count(x) / len(reduced_eigenvals) for x in F16]
		R16.<y> = F16[]
		a = (y^4 + y + 1).roots()[0][0]
		def Maarten(x):
			if x == 0:
				return 16/255
			if x in {a, a^2, a^3, a+1, a^3+a^2, a^2+1, a^3+a, a^3+a^2+a+1}:
				return 1/17
			if x in {a^2+a, a^3+a+1, a^2+a+1, a^3+a^2+a, a^3+a^2+1, a^3+1, 1}:
				return 1/15
		Maarten_histogram = [Maarten(x) for x in F16]
		diff_histogram = [abs(histogram[i] - Maarten_histogram[i]) for i in range(16)]
		diff_list.append( ( max(diff_histogram), label[:-1] ) )
		
diff_list = sorted(diff_list, key = lambda x: x[0])
