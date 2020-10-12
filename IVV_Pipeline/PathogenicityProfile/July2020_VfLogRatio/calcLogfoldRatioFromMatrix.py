import sys
import re
import math

genome_metadata_file = "Bacterial_human_pathogens_vs_nonhuman_genomes_sampled_20_per_species.txt"
gene_count_matrix_file = sys.argv[1]

genome_cat = {}
cat_count = {}
F = open(genome_metadata_file)
for line in F:
	genome, category, name = line.rstrip().split("\t")
	genome_cat[genome] = category
	if category not in cat_count:
		cat_count[category] = 0
	cat_count[category] += 1

sys.stderr.write("Num categories = {}\n".format(len(cat_count)))
for cat in cat_count:
	sys.stderr.write("Cat: {}\t{}\n".format(cat, cat_count[cat]))
cats = sorted(cat_count.keys())
cat_rat = float(cat_count[cats[0]])/len(genome_cat)
sys.stderr.write("Category ratio ({0}/total) = {1:.4f}\n".format(cats[0], cat_rat))

F = open(gene_count_matrix_file)
genomes = F.readline().rstrip().split("\t")
genomes.pop(0)
sys.stdout.write("VF_ID\tlog10({}/{})\tVirulence_factor\n".format(cats[0], cats[1]))
for line in F:
	g_count = line.rstrip().split("\t")
	vf = g_count.pop(0)
	vfid = re.search("\((\w+\d+)\)", vf).group(1)
	vf_cat_count = {cats[0]: 0, cats[1]: 0}
	for i, genome in enumerate(genomes):
		cat = genome_cat[genome]
		vf_cat_count[cat] += float(g_count[i])
	ratio = ((vf_cat_count[cats[0]]*(1-cat_rat))+1)/((vf_cat_count[cats[1]]*cat_rat)+1)
	log_ratio = math.log10(ratio)
	sys.stdout.write("{}\t{:.5f}\t{}\n".format(vfid, log_ratio, vf))


"""
==> Bacterial_human_pathogens_vs_nonhuman_genomes_sampled_20_per_species.txt <==
1311018.3	human	Acinetobacter baumannii
470.8890	human	Acinetobacter baumannii
470.9074	human	Acinetobacter baumannii
470.5522	human	Acinetobacter baumannii
470.9876	human	Acinetobacter baumannii
470.576	human	Acinetobacter baumannii
1116234.8	human	Acinetobacter baumannii
470.2270	human	Acinetobacter baumannii
1224747.3	human	Acinetobacter baumannii
470.2996	human	Acinetobacter baumannii
VF	100053.4	100053.6	100053.8	1001588.3
Trehalose-recycling ABC transporter (CVF651)	41	52	40	37
LOS (CVF494)	51	54	29	70
PDIM (phthiocerol dimycocerosate) and PGL (phenolic glycolipid) biosynthesis and transport (CVF288)	86	67	77	87
Pyoverdine (CVF551)	13	13	2	0
MtrCDE (CVF759)	24	24	13	11
Copper exporter (CVF658)	20	20	37	38
Flagella (CVF521)	103	121	82	52
Alginate regulation (CVF523)	43	43	33	45
RegX3 (CVF667)	36	19	35	33
"""
