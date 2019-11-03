import sys
sys.path.append('../..')

from paramagpy import protein, fit, dataparse, metal

import numpy as np
import quaternion

def unique_q(q):
	axis = quaternion.as_rotation_vector(q)
	mag = np.linalg.norm(axis)
	i, j, k = axis
	if i>=0 and j>=0 and k>=0:
		i = i
		j = j
		k = k
	elif i<0 and j>=0 and k>=0:
		i = np.pi+i
		j = j
		k = k


	return quaternion.from_rotation_vector([i,j,k])	



def normalise_axis(axis):
	mag = np.linalg.norm(axis)
	if mag>2*np.pi:
		axis = (axis / mag) * (mag % (2*np.pi))
	return axis



axrh = (10.,2.)
m0 = metal.Metal(axrh=axrh)
m1 = metal.Metal(axrh=axrh)

bools = []
for i in range(10):
	# m0.eulers = np.random.uniform(0, 2*np.pi, size=3)
	# q = quaternion.from_rotation_matrix(m0.rotationMatrix)
	
	i = np.random.uniform(-2*np.pi,       0)
	j = np.random.uniform(0       , 2*np.pi)
	k = np.random.uniform(0       , 2*np.pi)

	axis = normalise_axis([i,j,k])
	q = quaternion.from_rotation_vector(axis)
	uq = unique_q(q)
	m0.rotationMatrix = quaternion.as_rotation_matrix(q)
	m1.rotationMatrix = quaternion.as_rotation_matrix(uq)
	print(np.allclose(m0.tensor, m1.tensor))
	# print(m0.tensor - m1.tensor)








# # Load the PDB file
# prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# # Load the fitted tensor
# met = metal.load_tensor('../data_files/calbindin_Er_HN_PCS_tensor.txt')


# q = quaternion.from_rotation_matrix(met.rotationMatrix)

# print(q)

# for i in range(100):
# 	print(i)
# 	eulers = np.random.uniform(-2*np.pi, 2*np.pi, 3)
# 	unique_eulers = metal.unique_eulers(eulers)
# 	R = metal.euler_to_matrix(unique_eulers)

# 	q = quaternion.from_euler_angles(eulers)
# 	# angle = np.arccos(q.real)*2
# 	# if angle>np.pi:
# 	# 	new_angle = np.pi - angle
# 	# 	qlnew = np.quaternion()
# 	# 	ql = np.log(q)
# 	# 	qlnew.imag = ql.imag * 0.5*(new_angle/angle)

# 	# 	q = np.exp(qlnew)

# 	q_expected = quaternion.from_rotation_matrix(R)

# 	axis = np.log(q_expected).imag
# 	angle = 2*np.linalg.norm(axis)
	
# 	print(angle/np.pi)

# 	# axis2 = np.log(q).imag

# 	# print(axis, np.linalg.norm(axis)*2)
# 	# print(axis2, np.linalg.norm(axis2)*2)
	
	





# angle = np.pi/2.

# q = quaternion.from_euler_angles(0.,0.,angle)

# ql = np.log(q)

# print(np.arccos(q.real)/angle)

# print(abs(ql)/angle)

