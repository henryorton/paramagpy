s = """
2
Q
-2.35
0.32
-1.49
0.27
41
Q
2.56
0.33
0.56
0.21
3
I
-4.05
0.38
0.3
0.26
42
R
2.04
0.22


4
F
-3.58
0.42
1.45
0.30
43
L
2.41
0.41
0.33
0.28
5
V
-0.87
0.38
3.77
0.27
44
I
0.58
0.32
0.84
0.27
6
K
-1.34
0.43
1.25
0.38
45
F
1.3
0.28
3.82
0.26
7
T
0.89
0.25
3.62
0.20
47
G
-0.95
0.38
1.64
0.43
8
L
0.07
1.18
0.57
0.93
48
K
1.94
0.26
1.76
0.35
10
G
0.77
0.73
1.12
0.48
49
Q
0.34
0.23


11
K
0.11
0.72
-0.02
0.68
50
L
-0.17
0.33
-1.22
0.29
12
T
0.99
1.04
2.09
0.80
51
E
-3.86
0.29
-4.55
0.35
13
I
0.07
0.37
4.18
0.25
52
D
1.3
0.35
-2.63
0.20
14
T
-3.92
0.49
1.22
0.26
54
R
1.83
0.55


15
L
-4.58
0.42
1.39
0.25
55
T
-3.21
0.24


16
E
-4.09
0.57
-1.07
0.19
56
L
-0.9
0.67


17
V
-2.32
0.48
0.15
0.28
57
S
-1.19
0.27


18
E
-1.26
0.72
-0.92
0.49
58
D
-4.29
0.19


20
S
2.06
0.34


59
Y
-2.05
0.27


21
D
-0.34
0.79


60
N
2.31
0.21


25
N


-2.63
0.33
61
I
-0.07
0.31


26
V


1.89
0.34
62
Q
2.32
0.23
0.01
0.60
27
K


0.01
0.44
63
K
1.08
0.40


28
A


-3.52
0.21
64
E
-3.07
0.28
2.36
0.29
29
K


-0.63
0.34
65
S
2.99
0.28
-0.82
0.26
30
I


1.15
0.23
66
T
-3.62
0.35
2.23
0.29
31
Q


-2.16
0.19
67
L
-0.73
0.39
3.96
0.30
32
D


-3.82
0.16
68
H
-0.32
0.36
1.06
0.29
33
K


0.61
0.29
69
L
2.08
0.32


34
E
0.66
0.57
-1.34
0.22
70
V
2.63
0.27
-2.14
0.28
35
G
2.93
0.54
-4.32
0.25
71
L
2.77
0.14
-0.39
0.18
36
I
-0.84
0.37
-2.42
0.22
72
R
1.92
0.35
1.88
0.26
39
D


2.01
0.20
73
L
-0.47
0.75


40
Q
1.71
0.67
4.09
0.20
76
G
0.14
0.30
0.20
0.20
"""


ppp = d = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

ppp = {i:j for j,i in ppp.items()}

d1 = {}
d2 = {}


ss = s.split('\n')
for i in range(len(ss)):
	if ss[i].isalpha():
		key = ss[i-1], ss[i]
		v1a, v1b, v2a, v2b = ss[i+1:i+5]
		if v1a and v1b:
			d1[key] = v1a, v1b
		if v2a and v2b:
			d2[key] = v2a, v2b



with open('ubiquitin_a28c_c1_Tb_HN.rdc', 'w') as o:
	for key in sorted(list(d1), key=lambda x: int(x[0])):
		value, error = map(float,d1[key])
		s = int(key[0])
		line = "{0:<3d} N  {1:<3d} H  {2:5.2f}  {3:5.2f}\n".format(s,s,value,error)
		o.write(line)


with open('ubiquitin_s57c_c1_Tb_HN.rdc', 'w') as o:
	for key in sorted(list(d2), key=lambda x: int(x[0])):
		value, error = map(float,d2[key])
		s = int(key[0])
		line = "{0:<3d} N  {1:<3d} H  {2:5.2f}  {3:5.2f}\n".format(s,s,value,error)
		o.write(line)


def write_pales(fileName, data):
	hdr = "VARS   RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J D      DD    W\n"
	fmt = "%5d     %6s       %6s        %5d     %6s       %6s    %9.3f   %9.3f %.2f\n"

	with open(fileName, 'w') as o:
		o.write(hdr)
		o.write("FORMAT "+fmt)

		for key in sorted(list(data), key=lambda x: int(x[0])):
			seq, res = key
			seq = int(seq)
			res = ppp[res]
			value, error = map(float,data[key])
			line = fmt % (seq,res,'N',seq,res,'H',value,error,1.0)
			o.write(line)



write_pales('ubiquitin_a28c_c1_Tb_HN.tab', d1)
write_pales('ubiquitin_s57c_c1_Tb_HN.tab', d2)