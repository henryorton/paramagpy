def data_to_pdb(fileName, saveName):
	data = []
	with open(fileName) as f:
		for line in f:
			atom, *coords = line.split()
			coords = list(map(float, coords))
			data.append((atom, coords))

	s = ""
	for i, (atom, c) in enumerate(data):
		if 'P' in atom:
			element = 'P'
			continue
		elif 'O' in atom:
			element = 'O'
			continue
		elif 'N' in atom:
			element = 'N'
			continue
		elif 'C' in atom:
			element = 'C'
			continue
		else:
			element = 'H'

		fmt = "ATOM   {0:3d}  {1:<6s} RES A 1   {2:7.3f}{3:7.3f}{4:7.3f}  1.00  0.00  {5:1s}\n".format(
			i+1, atom, c[0], c[1], c[2], element)
		s += fmt
		print(fmt)

	# with open(saveName, 'w') as f:
		# f.write(s)


# ATOM     80  CA  GLY A   8      12.393  17.630  -3.793  1.00 14.06           C  
# ATOM     81  C   GLY A   8      12.178  18.142  -2.374  1.00 13.90           C  
# ATOM     82  O   GLY A   8      12.705  19.254  -1.990  1.00 14.66           O  
# ATOM     83  H   GLY A   8      11.773  15.623  -4.290  1.00 12.75           H  
# ATOM     84  HA1 GLY A   8      11.554  17.944  -4.414  1.00 14.06           H  
# ATOM     85  HA2 GLY A   8      13.329  18.045  -4.166  1.00 14.06           H  
# ATOM     86  N   ILE A   9      11.523  17.376  -1.523  1.00 11.84           N  
# ATOM     87  CA  ILE A   9      11.307  17.826  -0.127  1.00 11.32           C  

# data_to_pdb('Dy_R1.csv', 'Dy.pdb')
data_to_pdb('Tb_R1.csv', 'Tb.pdb')