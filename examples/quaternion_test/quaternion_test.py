from paramagpy import protein, fit, dataparse, metal

import numpy as np
import quaternion

# # Load the PDB file
# prot = protein.load_pdb('../data_files/4icbH_mut.pdb')

# # Load the fitted tensor
# met = metal.load_tensor('../data_files/calbindin_Er_HN_PCS_tensor.txt')


# q = quaternion.from_rotation_matrix(met.rotationMatrix)

# print(q)

for i in range(100):
	print(i)
	eulers = np.random.uniform(-2*np.pi, 2*np.pi, 3)
	unique_eulers = metal.unique_eulers(eulers)
	R = metal.euler_to_matrix(unique_eulers)

	q = quaternion.from_euler_angles(eulers)
	# angle = np.arccos(q.real)*2
	# if angle>np.pi:
	# 	new_angle = np.pi - angle
	# 	qlnew = np.quaternion()
	# 	ql = np.log(q)
	# 	qlnew.imag = ql.imag * 0.5*(new_angle/angle)

	# 	q = np.exp(qlnew)

	q_expected = quaternion.from_rotation_matrix(R)

	axis = np.log(q_expected).imag
	angle = 2*np.linalg.norm(axis)
	
	print(angle/np.pi)

	# axis2 = np.log(q).imag

	# print(axis, np.linalg.norm(axis)*2)
	# print(axis2, np.linalg.norm(axis2)*2)
	
	





# angle = np.pi/2.

# q = quaternion.from_euler_angles(0.,0.,angle)

# ql = np.log(q)

# print(np.arccos(q.real)/angle)

# print(abs(ql)/angle)

