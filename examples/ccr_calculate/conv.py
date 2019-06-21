d1 = []
d2 = []
with open("data_from_paper.txt") as o:
	for line in o:
		a, b, c, d = line.split(",")
		seq = int(a)
		unpack = lambda x: map(lambda y: float(y.strip()), x.split('±'))

		try:
			dv, de = unpack(b)
			pv, pe = unpack(c)
			d1.append((seq, pv-dv,(de**2+pe**2)**0.5))
		except:
			pass

		try:
			dv, de = unpack(b)
			pv, pe = unpack(d)
			d2.append((seq, pv-dv,(de**2+pe**2)**0.5))
		except:
			pass


with open("myo_cn.ccr", 'w') as o:
	for s, v, e in d1:
		line = "{0:3d} H{1:4d} N {2:4.1f} {3:3.1f}\n".format(s, s, v, e)
		o.write(line)

with open("myo_f.ccr", 'w') as o:
	for s, v, e in d2:
		line = "{0:3d} H{1:4d} N {2:5.1f} {3:3.1f}\n".format(s, s, v, e)
		o.write(line)