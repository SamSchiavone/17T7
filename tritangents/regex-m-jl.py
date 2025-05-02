import re

with open("sols_refined.out", "r") as file:
    pol_data = file.read()

# Replace every occurrence of " +/- number E- number2" with "Pnumber2"
# Regex: match a space, then + or -, then optional spaces, then a number, then "E-", then a number
pattern_replace = r'\s\+/-\s*\d+(?:\.\d+)?[eE]-?(\d+)'

# Replacement string: "P<exponent>"
converted_pol_simple = re.sub(pattern_replace, r'P\1', pol_data)

converted_pol_simple = converted_pol_simple.replace("AcbVector[", "CC<I> := ComplexField(18000);\nvec := ")
converted_pol_simple = converted_pol_simple.replace("] + [", " + ")

converted_pol_simple = converted_pol_simple.replace("]im, [", "*I, ")
converted_pol_simple = converted_pol_simple.replace("]im], [[", "*I], [")

converted_pol_simple = converted_pol_simple.replace("]im]]", "*I]];\n")


# Save the result to a new file
output_path_final = "./sols_refined.m";
with open(output_path_final, "w") as file:
    file.write(converted_pol_simple)
