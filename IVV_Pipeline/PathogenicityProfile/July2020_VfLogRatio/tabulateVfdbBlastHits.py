import sys
import re

""" 
Read multiple blastx files for multiple genomes, count each VF in bracket-defined field of 6th tab-seqarated field)
Get genome ID from blastx file name.
Write out matrix of VF count per genome to STDOUT
"""
genome_vf_count = {}
genome_list = []
vf_count = {}
for blastx_file in sys.argv[1:]:
	genome = blastx_file
	m = re.search(r"(\d+\.\d+)", blastx_file)
	if m:
		genome = m.group(1)
	genome_list.append(genome)
	genome_vf_count[genome] = {}
	sys.stderr.write(genome+"\n")
	F = open(blastx_file)
	for line in F:
		fields = line.split("\t")
		if len(fields) < 6:
			continue
		m = re.search("\[([^\]]+\))\]", fields[5])
		if m:
			vf = m.group(1)
			if vf not in genome_vf_count[genome]:
				genome_vf_count[genome][vf] = 0
				if vf not in vf_count:
					vf_count[vf] = 0
			genome_vf_count[genome][vf] += 1
			vf_count[vf] += 1
		else:
			sys.stderr.write("cannot parse {}\n".format(fields[5]))

sys.stdout.write("VF\t"+"\t".join(genome_list)+"\n")
for vf in sorted(vf_count, key = vf_count.get, reverse=True):
	sys.stdout.write(vf)
	for genome in genome_list:
		count = 0
		if vf in genome_vf_count[genome]:
			count = genome_vf_count[genome][vf]
		sys.stdout.write("\t{:d}".format(count))
	sys.stdout.write("\n")



"""
==> ffvf_means_logodds.txt <==
Foe	Friend	logodds
irtB	4.5210389021209	5.81954371878978	-0.252475111565807
rtxB	4.93716683679256	6.06501517673553	-0.205741626424912
irp6C	0.684019507769082	1.1099578967982	-0.484034830902171

==> ref_rep_vs_vfdb/312309.11_vfdb.blastx <==
NC_006841	VFG007277(gi:59713515)	711	0.0e+00	1387.9	VFG007277(gi:59713515) (hutR) hemin receptor [Heme receptors (CVF278)] [Vibrio fischeri ES114]
NC_006841	VFG006872(gi:59714051)	497	9.6e-289	999.2	VFG006872(gi:59714051) (tcpT) TcpT, toxin co-regulated pilus biosynthesis protein [Toxin-coregulated pilus (type IVB pilus) (CVF256)] [Vibrio fischeri ES114]
NC_006841	VFG006861(gi:59714055)	484	4.2e-260	904.0	VFG006861(gi:59714055) (tcpC) TcpC, toxin co-regulated pilus biosynthesis outer membrane protein [Toxin-coregulated pilus (type IVB pilus) (CVF256)] [Vibrio fischeri ES114]
NC_006841	VFG006855(gi:59714057)	445	3.3e-257	894.4	VFG006855(gi:59714057) (tcpB) TcpB, toxin co-regulated pilus biosynthesis protein [Toxin-coregulated pilus (type IVB pilus) (CVF256)] [Vibrio fischeri ES114]
==> 1144305.3_vfdb.blastx <==
AKKE01000069	VFG036966(gi:385342653)	508	3.9e-102	375.2	VFG036966(gi:385342653) (farB) antibacterial fatty acid resistance protein B [FarAB (CVF758)] [Neisseria meningitidis M01-240149]
AKKE01000069	VFG036964(gi:385327678)	508	6.7e-102	374.4	VFG036964(gi:385327678) (farB) multidrug resistance translocase [FarAB (CVF758)] [Neisseria meningitidis alpha710]
AKKE01000069	VFG036955(gi:218768903)	508	1.9e-101	372.9	VFG036955(gi:218768903) (farB) multidrug resistance translocase [FarAB (CVF758)] [Neisseria meningitidis Z2491]
AKKE01000069	VFG036961(gi:385338790)	508	1.9e-101	372.9	VFG036961(gi:385338790) (farB) multidrug resistance protein B [FarAB (CVF758)] [Neisseria meningitidis WUE 2594]
AKKE01000069	VFG036963(gi:254805645)	508	1.9e-101	372.9	VFG036963(gi:254805645) (farB) putative efflux pump protein [FarAB (CVF758)] [Neisseria meningitidis alpha14]
AKKE01000069	VFG036972(gi:313667722)	509	1.9e-101	372.9	VFG036972(gi:313667722) (farB) multidrug resistance translocase [FarAB (CVF758)] [Neisseria lactamica 020-06]
AKKE01000069	VFG036958(gi:121635543)	508	4.3e-101	371.7	VFG036958(gi:121635543) (farB) multidrug resistance translocase [FarAB (CVF758)] [Neisseria meningitidis FAM18]
AKKE01000069	VFG036969(gi:385854493)	508	4.3e-101	371.7	VFG036969(gi:385854493) (farB) antibacterial fatty acid resistance protein B [FarAB (CVF758)] [Neisseria meningitidis M01-240355]
AKKE01000069	VFG036956(gb|NP_273368)	508	5.6e-101	371.3	VFG036956(gb|NP_273368) (farB) fatty acid efflux system protein FarB [FarAB (VF0450)] [Neisseria meningitidis MC58]
AKKE01000069	VFG036962(gi:385323459)	508	5.6e-101	371.3	VFG036962(gi:385323459) (farB) multidrug resistance protein B [FarAB (CVF758)] [Neisseria meningitidis 8013]
"""
