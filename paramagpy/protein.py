import heapq
import warnings
from operator import attrgetter
from random import randint

import numpy as np
import quaternion as quat
from Bio.PDB import PDBParser
from Bio.PDB.Atom import Atom, DisorderedAtom
from Bio.PDB.PDBExceptions import PDBConstructionWarning, PDBConstructionException
from Bio.PDB.Residue import Residue, DisorderedResidue
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder


def rotation_matrix(axis, theta):
    """Return the rotation matrix associated with counterclockwise
    rotation about the given axis by theta radians.

    Parameters
    ----------
    axis : array of floats
        the [x,y,z] axis for rotation.
    theta : angle of rotation about the axis

    Returns
    -------
    matrix : numpy 3x3 matrix object
        the rotation matrix
    """
    axis = np.array(axis)
    axis /= np.linalg.norm(axis)
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])


class CustomAtom(Atom):
    MU0 = 4 * np.pi * 1E-7
    HBAR = 1.0546E-34

    gyro_lib = {
        'H': 2 * np.pi * 42.576E6,
        'N': 2 * np.pi * -4.316E6,
        'C': 2 * np.pi * 10.705E6}  # rad/s/T

    csa_lib = {
        'H': (np.array([-5.8, 0.0, 5.8]) * 1E-6, 8. * (np.pi / 180.)),
        'N': (np.array([-62.8, -45.7, 108.5]) * 1E-6, 19. * (np.pi / 180.)),
        'C': (np.array([-86.5, 11.8, 74.7]) * 1E-6, 38. * (np.pi / 180.))}

    valency_lib = {
        'H': 1,
        'C': 4,
        'N': 3,
        'O': 2,
        'S': 2,
        'CA': 2
    }

    # Priority order for choosing the atom to measure the dihedral angle from
    atomic_number_lib = {
        'H': 1,
        'C': 6,
        'N': 7,
        'O': 8,
        'S': 16
    }

    """docstring for CustomAtom"""

    def __init__(self, *arg, **kwargs):
        super().__init__(*arg, **kwargs)
        self.coord = np.asarray(self.coord, dtype=np.float64)
        self.gamma = self.gyro_lib.get(self.element, 0.0)
        self._csa = None
        self.valency = self.valency_lib.get(self.element)
        self.bonded_atoms = None
        self.atomic_number = CustomAtom.atomic_number_lib.get(self.element)

    def __repr__(self):
        return "<Atom {0:3d}-{1:}>".format(self.parent.id[1], self.name)

    def top(self):
        return self.parent.parent.parent.parent

    @property
    def position(self):
        return self.coord * 1E-10

    @position.setter
    def position(self, value):
        self.coord = value * 1E10

    @property
    def csa(self):
        """
        Get the CSA tensor at the nuclear position
        This uses the geometry of neighbouring atoms
        and a standard library from Bax J. Am. Chem. Soc. 2000

        Returns
        -------
        matrix : 3x3 array
            the CSA tensor in the PDB frame
            if appropriate nuclear positions are not
            available <None> is returned.
        """

        if self._csa is not None:
            return self._csa

        def norm(x):
            return x / np.linalg.norm(x)

        res = self.parent
        resid = res.id
        respid = resid[0], resid[1] - 1, resid[2]
        resnid = resid[0], resid[1] + 1, resid[2]
        resp = res.parent.child_dict.get(respid)
        resn = res.parent.child_dict.get(resnid)

        pas, beta = self.csa_lib.get(self.name, (None, None))
        if resp:
            Hcond = self.element == 'H', 'N' in res, 'C' in resp, beta
            Ncond = self.element == 'N', 'H' in res, 'C' in resp, beta
        else:
            Hcond = (None,)
            Ncond = (None,)
        if resn:
            Ccond = self.element == 'C', 'H' in resn, 'N' in resn, beta
        else:
            Ccond = (None,)

        if all(Hcond):
            NC_vec = resp['C'].coord - res['N'].coord
            NH_vec = res['H'].coord - res['N'].coord
            z = norm(np.cross(NC_vec, NH_vec))
            R = rotation_matrix(-z, beta)
            x = norm(R.dot(NH_vec))
            y = norm(np.cross(z, x))

        elif all(Ncond):
            NC_vec = resp['C'].coord - res['N'].coord
            NH_vec = res['H'].coord - res['N'].coord
            y = norm(np.cross(NC_vec, NH_vec))
            R = rotation_matrix(-y, beta)
            z = norm(R.dot(NH_vec))
            x = norm(np.cross(y, z))

        elif all(Ccond):
            CN_vec = resn['N'].coord - res['C'].coord
            NH_vec = resn['H'].coord - resn['N'].coord
            x = norm(np.cross(NH_vec, CN_vec))
            R = rotation_matrix(x, beta)
            z = norm(R.dot(CN_vec))
            y = norm(np.cross(z, x))

        else:
            return np.zeros(9).reshape(3, 3)
        transform = np.vstack([x, y, z]).T
        tensor = transform.dot(np.diag(pas)).dot(transform.T)
        return tensor

    @csa.setter
    def csa(self, newTensor):
        if newTensor is None:
            self._csa = None
            return
        try:
            assert newTensor.shape == (3, 3)
        except (AttributeError, AssertionError):
            print("The specified CSA tensor does not have the correct format")
            raise
        self._csa = newTensor

    def dipole_shift_tensor(self, position):
        """
        Calculate the magnetic field shielding tensor at the given postition
        due to the nuclear dipole

        Assumes nuclear spin 1/2

        Parameters
        ----------
        position : array floats
            the position (x, y, z) in meters

        Returns
        -------
        dipole_shielding_tensor : 3x3 array
            the tensor describing magnetic shielding at the given position
        """
        pos = np.array(position, dtype=float) - self.position
        distance = np.linalg.norm(pos)
        preFactor = (self.MU0 * self.gamma * self.HBAR * 0.5) / (4. * np.pi)
        p1 = (1. / distance ** 5) * np.kron(pos, pos).reshape(3, 3)
        p2 = (1. / distance ** 3) * np.identity(3)
        return (preFactor * (3. * p1 - p2))

    def bonded_to(self, valency=None, recompute=False):
        """
        Compute the list of atoms which could be bonded to this atom

        Note: The function doesn't guarantee that the list returned is always correct. Use with caution.

        Parameters
        ----------
        valency : int
            The valency of the atom this function is being called from. If none, it will try to use the
            default value if it exists. If the atom is charged, please add the charge to valency to get
            accurate results.

        recompute : bool
            The list of atoms is only computed the first time and stored for later use. Calls succeeding
            the first, return the stored value instead. Set this argument to True if you would like to
            recompute the list.

        Returns
        -------
        atoms_bonded_to : list of CustomAtom
            The list of atoms which could be bonded to this atom
        """
        return self.parent.bonded_to(self, self.valency if valency is None else valency, recompute)


class CustomStructure(Structure):
    """This is an overload hack of the BioPython Structure object"""

    def __init__(self, *arg, **kwargs):
        super().__init__(*arg, **kwargs)

    def parse(self, dataValues, models=None):
        used = set([])
        data = []

        if type(models) == int:
            chains = self[models].get_chains()
        elif type(models) in (list, tuple):
            chains = []
            for m in models:
                chains += self[m].get_chains()
        else:
            chains = self.get_chains()

        if dataValues.dtype in ('PCS', 'PRE'):
            for chain in chains:
                for key in dataValues:
                    seq, name = key
                    if seq in chain:
                        resi = chain[seq]
                        if name in resi:
                            atom = resi[name]
                            data.append((atom, *dataValues[key]))
                            used.add(key)

        elif dataValues.dtype in ('RDC', 'CCR'):
            for chain in chains:
                for key in dataValues:
                    (seq1, name1), (seq2, name2) = key
                    if seq1 in chain and seq2 in chain:
                        resi1 = chain[seq1]
                        resi2 = chain[seq2]
                        if name1 in resi1 and name2 in resi2:
                            atom1 = resi1[name1]
                            atom2 = resi2[name2]
                            data.append((atom1, atom2, *dataValues[key]))
                            used.add(key)

        unused = set(dataValues) - used
        if unused:
            message = "WARNING: Some values were not parsed to {}:"
            print(message.format(self.id))
            print(list(unused))
        return data


class CustomResidue(Residue):
    """Paramagpy wrapper for BioPython's Residue entity"""

    # Atoms along the side-chain about which the atoms are to be rotated to generate different rotamers
    # Prepend ['CA'] for all molecules
    side_chain_lib = {
        'GLY': [],
        'ALA': [],
        'SER': ['CB'],
        'THR': ['CB'],
        'CYS': ['CB'],
        'VAL': ['CB'],
        'LEU': ['CB', 'CG'],
        'ILE': ['CB', 'CG1'],
        'MET': ['CB', 'CG', 'SD'],
        'PRO': [],
        'PHE': ['CB', 'CG'],
        'TYR': ['CB', 'CG'],
        'TRP': ['CB', 'CG'],
        'ASP': ['CB', 'CG'],
        'GLU': ['CB', 'CG', 'CD'],
        'ASN': ['CB', 'CG'],
        'GLN': ['CB', 'CG', 'CD'],
        'HIS': ['CB', 'CG'],
        'ARG': ['CB', 'CG', 'CD', 'NE'],
        'LYS': ['CB', 'CG', 'CD', 'CE']
    }

    # Residue backbone
    back_bone = {'N', 'C', 'O', 'CA'}

    def __init__(self, *arg, **kwargs):
        super().__init__(*arg, **kwargs)
        # These variables are set by fit.FitPCSToRotamer as and when required
        self._noise = None
        self._dihedral_full = None
        self._metal = None
        self.pcs_data = None
        self._min_pcs = None
        self._rot_path = None

    def bonded_to(self, source_atom, valency=None, recompute=False):
        """
        Compute the list of atoms which could be bonded to the source_atom in this residue

        Note: The function doesn't guarantee that the list returned is always correct. Use with caution.

        Parameters
        ----------
        source_atom: CustomAtom
            An atom part of this residue, the atoms bonded to which are returned

        valency : int
            The valency of the atom this function is being called from. If none, it will try to use the
            default value if it exists. If the atom is charged, please add the charge to valency to get
            accurate results.

        recompute : bool
            The list of atoms is only computed the first time and stored for later use. Calls succeeding
            the first, return the stored value instead. Set this argument to True if you would like to
            recompute the list.

        Returns
        -------
        atoms_bonded_to : list of CustomAtom
            The list of atoms which could be bonded to this atom
        """

        def dist(from_atom):
            return np.abs(from_atom - source_atom)

        def quick_select(start, end, val):
            if start >= end:
                return None

            mid = partition(start, end, randint(start, end))
            bucket_a_len = mid - start + 1

            if valency < bucket_a_len:
                quick_select(start, mid - 1, val)
            elif valency > bucket_a_len:
                quick_select(mid + 1, end, val - bucket_a_len)

        def partition(start, end, rnd):
            atoms[start], atoms[rnd] = atoms[rnd], atoms[start]
            _start = start
            pivot = dist(atoms[start])

            start += 1
            while True:
                while start < end and dist(atoms[start]) < pivot:
                    start += 1
                while start <= end and dist(atoms[end]) >= pivot:
                    end -= 1
                if start >= end:
                    break
                atoms[start], atoms[end] = atoms[end], atoms[start]

            atoms[_start], atoms[end] = atoms[end], atoms[_start]
            return end

        # If computed already, just return the list
        if source_atom.bonded_atoms is not None and not recompute:
            return source_atom.bonded_atoms

        # If id is supplied instead of CustomAtom, convert accordingly
        if type(source_atom) is str:
            source_atom = self[source_atom] if self.has_id(source_atom) else None

        # Error handling

        # If source_atom is not part of the residue, raise an exception
        if source_atom not in self.child_list and not source_atom.disordered_flag:
            raise ValueError("source_atom is not part of this residue, call this function with a valid source")

        # If valency isn't provided, use the CustomAtom instance's valency
        # If still not available, throw an exception (or possibly return atoms at distance < 1.8E-10)
        valency = source_atom.valency if valency is None else valency
        if valency is None:
            raise ValueError(
                "Couldn't obtain valency for this source atom, call this function with the valency parameter")

        # Only look for the nearest atoms at distance < 1.8E-10
        # dist > 0 ensures the source_atom itself isn't considered
        atoms = [atom for atom in self.get_atoms() if 0 < dist(atom) < 1.95]
        quick_select(0, len(atoms) - 1, valency)
        source_atom.bonded_atoms = atoms[:valency]
        return source_atom.bonded_atoms

    def set_dihedral(self, theta_vector):
        """
        Set the given dihedral angles in the residue. The dihedral angle has to be in the range (-pi, pi].
        Uses the right-hand rule as looking down from the backbone of the residue. The convention used for
        finding the dihedral angle about a bond is to find the atom with the highest atomic number attached
        both the atoms (the bond between whom the dihedral angle is being discussed) and use this as the four
        atoms for the dihedral calculations. If there is a tie in the highest atomic number, the atom with the
        highest lexicographical name in the PDB file among these tied atoms is chosen.

        Parameters
        ----------
        theta_vector : numpy.ndarray
         nd array of shape (#bonds on side chain which are rotated)
            A numpy array corresponding to the dihedral angles which have to be set in the residue

        """

        rot_path = self.__get_rot_path()

        # Error Handling - Exit raising an AttributeError if the theta vector supplied is of invalid dimension
        if len(theta_vector) != len(rot_path) - 1:
            raise AttributeError(
                f"The delta theta vector supplied is of invalid length. Expected: {len(rot_path) - 1}, "
                f"Actual: {len(theta_vector)}")

        current_dihedral_vector = self.get_dihedral_full()
        delta_theta_vector = theta_vector - current_dihedral_vector
        self.set_delta_dihedral(delta_theta_vector)

    @staticmethod
    def get_dihedral(atom_a, atom_b):
        """
        Get the given dihedral angles in the residue. The dihedral angle has to be in the range (-pi, pi].
        Uses the right-hand rule as looking down from the backbone of the residue. The convention used for
        finding the dihedral angle about the bond atom_a and atom_b is to find the atom with the highest
        atomic number attached both the atoms and use this as the four atoms for the dihedral calculations.
        If there is a tie in the highest atomic number, the atom with the highest lexicographical name in
        the PDB file among these tied atoms is chosen.

        Parameters
        ----------
        atom_a : CustomAtom
            Atom from which the dihedral angle is being measured.
        atom_b : CustomAtom
            Atom from which the dihedral angle is being measured.

        Returns
        -------
        dihedral_angle : float
            The dihedral angle in radians

        """
        atom_a_bonded_to = atom_a.bonded_to().copy()
        atom_b_bonded_to = atom_b.bonded_to().copy()

        # Error Handling
        if atom_a not in atom_b_bonded_to or atom_b not in atom_a_bonded_to:
            raise AttributeError("Atom A and atom B have to bonded to each other for this to work")

        atom_a_bonded_to.remove(atom_b)
        atom_b_bonded_to.remove(atom_a)

        list.sort(atom_a_bonded_to, key=attrgetter('atomic_number', 'name'), reverse=True)
        list.sort(atom_b_bonded_to, key=attrgetter('atomic_number', 'name'), reverse=True)

        if not atom_a_bonded_to:
            raise Exception(f"Could not find the neighbors of {atom_a} to calculate the dihedral angle")

        if not atom_b_bonded_to:
            raise Exception(f"Could not find the neighbors of {atom_b} to calculate the dihedral angle")

        return CustomResidue.__calc_dihedral(atom_a_bonded_to[0].coord, atom_a.coord, atom_b.coord,
                                             atom_b_bonded_to[0].coord)

    def get_dihedral_full(self):
        """
        Get the list of all the dihedral angle (in radians) in the residue

        Returns
        -------
        dihedral_angle : float
            List of all the dihedral angle (in radians) in the residue

        """
        rot_path = self.__get_rot_path()
        return np.array([CustomResidue.get_dihedral(rot_path[i], rot_path[i + 1]) for i in range(len(rot_path) - 1)])

    def set_delta_dihedral(self, delta_theta_vector):
        """
        Increase each dihedral angle by the corresponding delta_theta value. Increase in angle is according
        to the right-hand rule looking from the back-bone of the residue.

        Parameters
        ----------
        delta_theta_vector : numpy.ndarray
            ndarray of shape (# bonds on side chain which are rotated)
            A numpy array corresponding to the dihedral angles (in radians) which have to be set
            in the residue

        """
        # Error Handling - HETATM seems to be causing a problem, raising an exception until we know for sure
        # what is to be done
        if self.get_resname() not in CustomResidue.side_chain_lib:
            raise Exception("Residue is not an amino acid, could be a hetero atom")

        # If there is an amine-H, don't rotate it and make it a part of the backbone
        back_bone_atoms = set(map(lambda x: self[x], self.back_bone))
        iter_bb_atoms = set(back_bone_atoms)

        for atom in iter_bb_atoms:
            for nbr in atom.bonded_to():
                back_bone_atoms.add(nbr)

        rot_path = self.__get_rot_path()

        atoms_to_rotate = set(atom for atom in self.get_atoms() if atom not in back_bone_atoms)
        atoms_position = {x: 999 for x in atoms_to_rotate}
        for i in range(len(rot_path) - 1):
            for atom in rot_path[i + 1].bonded_to():
                if atom in atoms_to_rotate:
                    atoms_position[atom] = i
                atoms_to_rotate.discard(atom)

        # Error Handling - Exit raising an AttributeError if the theta vector supplied is of invalid dimension
        if len(delta_theta_vector) != len(rot_path) - 1:
            raise AttributeError(
                f"The delta theta vector supplied is of invalid length. Expected: {len(rot_path) - 1}, "
                f"Actual: {len(delta_theta_vector)}")

        for i in range(len(delta_theta_vector)):
            if delta_theta_vector[i] == 0:
                continue

            # The rotation axis is "along the previous atom in rot_path and the current" (aka bond b/w the two)
            rot_vector = rot_path[i + 1].coord - rot_path[i].coord
            q = CustomResidue.__rot_quat(delta_theta_vector[i], rot_vector)

            self.__rotate_atoms(atoms_position, i, q, rot_path[i].coord)

    def grid_search_rotamer(self, rotation_param, fit_pcs=False, top_n=1):
        """
        Perform a grid search generating rotamers using the parameters given by rotation_param

        Parameters
        ----------
        rotation_param : numpy.ndarray
            ndarray of shape (# bonds on side chain which are rotated, 3)
            For each bond, 3 parameters are required. The starting dihedral angle, the ending dihedral angle
            and the number of steps between start and finish (including the start and end angle)

        fit_pcs : bool
            If True, the function tries to calculate the rotamer with PCS values closest to the experimental
            PCS values (see :py:meth: `paramagpy.protein.CustomStructure.parse`)

        top_n : int
            Returns the top_n rotamers which fits best to the experimental PCS values

        """
        # Error Handling - HETATM seems to be causing a problem, raising an exception until we know for sure
        # what is to be done
        if self.get_resname() not in CustomResidue.side_chain_lib:
            raise Exception("Residue is not an amino acid, could be a hetero atom")

        for _p in rotation_param:
            if not -np.pi < _p[0] <= np.pi or not -np.pi < _p[1] <= np.pi:
                raise AttributeError("The dihedral angles should be in range (-pi, pi]")
            if _p[2] < 0:
                raise AttributeError(
                    "The steps for rotation cannot be negative. Note: 0 steps denotes no rotation about that bond")

        # If there is an amine-H, don't rotate it and make it a part of the backbone
        back_bone_atoms = set(map(lambda x: self[x], self.back_bone))
        iter_bb_atoms = set(back_bone_atoms)  # Copy to a new set to iterate over this while changing

        for atom in iter_bb_atoms:
            for nbr in atom.bonded_to():
                back_bone_atoms.add(nbr)

        rot_path = self.__get_rot_path()

        if len(rot_path) == 1:
            return None

        atoms_to_rotate = set(atom for atom in self.get_atoms() if atom not in back_bone_atoms)
        atoms_position = {x: 999 for x in atoms_to_rotate}
        for i in range(len(rot_path) - 1):
            for atom in rot_path[i + 1].bonded_to():
                if atom in atoms_to_rotate:
                    atoms_position[atom] = i
                atoms_to_rotate.discard(atom)

        self.set_dihedral(rotation_param[:, 0])
        self._min_pcs = []
        self._dihedral_full = np.array(rotation_param[:, 0])

        _param1 = rotation_param[0]
        q1 = CustomResidue.__rot_quat(CustomResidue.__get_angle_increment(_param1[0], _param1[1], _param1[2]),
                                      rot_path[1].coord - rot_path[0].coord)
        self.__grid_search_rotamer_helper(q1, 0, rotation_param, rot_path, atoms_position, fit_pcs, top_n)

        return self._min_pcs if fit_pcs and len(self._min_pcs) > 0 else None

    def __grid_search_rotamer_helper(self, q, i, rotation_param, rot_path, atoms_position, fit_pcs=False, top_n=1):
        _param = rotation_param[i]
        for j in range(1, int(_param[2])):
            if i + 1 < len(rot_path) - 1:
                _param1 = rotation_param[i + 1]
                q_next = CustomResidue.__rot_quat(
                    CustomResidue.__get_angle_increment(_param1[0], _param1[1], _param1[2]),
                    rot_path[i + 2].coord - rot_path[i + 1].coord)
                self.__grid_search_rotamer_helper(q_next, i + 1, rotation_param, rot_path,
                                                  atoms_position, fit_pcs, top_n)
            self._dihedral_full[i] += CustomResidue.__get_angle_increment(_param[0], _param[1], _param[2])
            self._dihedral_full[i] -= 2 * np.pi if self._dihedral_full[i] > np.pi else 0
            self.__rotate_atoms(atoms_position, i, q, rot_path[i].coord, fit_pcs=i + 1 == len(rot_path) - 1,
                                top_n=top_n)

        if int(_param[2]) > 0:
            if i + 1 < len(rot_path) - 1:
                _param1 = rotation_param[i + 1]
                q_next = CustomResidue.__rot_quat(
                    CustomResidue.__get_angle_increment(_param1[0], _param1[1], _param1[2]),
                    rot_path[i + 2].coord - rot_path[i + 1].coord)
                self.__grid_search_rotamer_helper(q_next, i + 1, rotation_param, rot_path,
                                                  atoms_position, fit_pcs, top_n)
            q_reset = CustomResidue.__rot_quat(2 * np.pi - (_param[1] - _param[0]),
                                               rot_path[i + 1].coord - rot_path[i].coord)
            self._dihedral_full[i] += 2 * np.pi - (_param[1] - _param[0])
            self._dihedral_full[i] -= 2 * np.pi if self._dihedral_full[i] > np.pi else 0
            self.__rotate_atoms(atoms_position, i, q_reset, rot_path[i].coord, fit_pcs=i + 1 == len(rot_path) - 1,
                                top_n=top_n)

    def __rotate_atoms(self, atoms_position, i, q, origin, fit_pcs=False, top_n=1):
        _q = quat.as_quat_array(np.empty(4))
        _q.real = 0
        for atom in atoms_position:
            if i > atoms_position[atom]:
                continue
            _q.imag = atom.coord - origin
            atom_coord_shifted = (q * _q * q.conj())
            atom.set_coord(atom_coord_shifted.imag + origin)

        if fit_pcs:
            coord_matrix = np.empty((len(self.pcs_data[0]), 3))

            for idx, tup in enumerate(self.pcs_data[0]):
                coord_matrix[idx] = tup.position

            pcs_calc = self._metal.fast_pcs(coord_matrix)
            pcs_calc_binned = np.bincount(self.pcs_data[3], weights=pcs_calc)

            pcs_exp_binned = np.bincount(self.pcs_data[3], weights=self.pcs_data[1])
            bin_count = np.bincount(self.pcs_data[3])
            bin_count[bin_count == 0] = 1

            _x = (pcs_exp_binned - pcs_calc_binned) / bin_count
            if not self._noise:
                if not np.any(np.abs(self.pcs_data[2]) < 1E-6):
                    pcs_err_binned = np.bincount(self.pcs_data[3], weights=self.pcs_data[2])
                    _x = (pcs_exp_binned - pcs_calc_binned) / pcs_err_binned
                else:
                    warnings.warn(
                        "Cannot calculate Mahalanobis distance because at least one of the PCS errors is zero. "
                        "Using RMSD instead.")

            pcs_dist = np.linalg.norm(_x)/len(pcs_calc_binned)

            try:
                if len(self._min_pcs) < top_n or top_n == -1:
                    heapq.heappush(self._min_pcs, (-pcs_dist, np.array(self._dihedral_full), pcs_calc))
                elif len(self._min_pcs) == top_n and -self._min_pcs[0][0] > pcs_dist:
                    heapq.heappop(self._min_pcs)
                    heapq.heappush(self._min_pcs, (-pcs_dist, np.array(self._dihedral_full), pcs_calc))
            except ValueError:
                pass

    def __get_rot_path(self):
        if self._rot_path is None:
            self._rot_path = [self['CA']] + list(
                map(lambda x: self[x], CustomResidue.side_chain_lib[self.get_resname()]))
        return self._rot_path

    @staticmethod
    def __get_angle_increment(start, stop, steps):
        if start == stop or int(steps) <= 1:
            return 0
        diff = (stop - start)
        if diff < 0:
            diff += 2 * np.pi
        return diff / int(steps - 1)

    @staticmethod
    def __rot_quat(theta, vector):
        v = vector / np.linalg.norm(vector)
        s = np.sin(theta / 2)
        c = np.cos(theta / 2)
        _q = quat.as_quat_array(np.empty(4))
        _q.real = c
        _q.imag = v * s
        return _q

    @staticmethod
    def __calc_dihedral(v1, v2, v3, v4):
        b = np.diff(np.array([v1, v2, v3, v4]), axis=0)
        n1, n2 = np.cross(b[0], b[1]), np.cross(b[1], b[2])

        # Our new co-ordinate frame (x, y, z) := (n1, n2, m1)
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)
        b[1] /= np.linalg.norm(b[1])
        m1 = np.cross(n1, b[1])

        return -np.arctan2(m1.dot(n2), n1.dot(n2))


class CustomStructureBuilder(StructureBuilder):
    """This is an overload hack of BioPython's CustomStructureBuilder"""

    def __init__(self, *arg, **kwargs):
        super().__init__()
        self.structure = None

    def init_structure(self, structure_id):
        self.structure = CustomStructure(structure_id)

    def init_residue(self, resname, field, resseq, icode):
        """Create a new Residue object.

        Arguments:
         - resname - string, e.g. "ASN"
         - field - hetero flag, "W" for waters, "H" for
           hetero residues, otherwise blank.
         - resseq - int, sequence identifier
         - icode - string, insertion code

        """
        if field != " ":
            if field == "H":
                # The hetero field consists of H_ + the residue name (e.g. H_FUC)
                field = "H_" + resname
        res_id = (field, resseq, icode)
        if field == " ":
            if self.chain.has_id(res_id):
                # There already is a residue with the id (field, resseq, icode).
                # This only makes sense in the case of a point mutation.
                warnings.warn("WARNING: Residue ('%s', %i, '%s') "
                              "redefined at line %i."
                              % (field, resseq, icode, self.line_counter),
                              PDBConstructionWarning)
                duplicate_residue = self.chain[res_id]
                if duplicate_residue.is_disordered() == 2:
                    # The residue in the chain is a DisorderedResidue object.
                    # So just add the last Residue object.
                    if duplicate_residue.disordered_has_id(resname):
                        # The residue was already made
                        self.residue = duplicate_residue
                        duplicate_residue.disordered_select(resname)
                    else:
                        # Make a new residue and add it to the already
                        # present DisorderedResidue
                        new_residue = CustomResidue(res_id, resname, self.segid)
                        duplicate_residue.disordered_add(new_residue)
                        self.residue = duplicate_residue
                        return
                else:
                    if resname == duplicate_residue.resname:
                        warnings.warn("WARNING: Residue ('%s', %i, '%s','%s')"
                                      " already defined with the same name "
                                      "at line  %i."
                                      % (field, resseq, icode, resname,
                                         self.line_counter),
                                      PDBConstructionWarning)
                        self.residue = duplicate_residue
                        return
                    # Make a new DisorderedResidue object and put all
                    # the Residue objects with the id (field, resseq, icode) in it.
                    # These residues each should have non-blank altlocs for all their atoms.
                    # If not, the PDB file probably contains an error.
                    if not self._is_completely_disordered(duplicate_residue):
                        # if this exception is ignored, a residue will be missing
                        self.residue = None
                        raise PDBConstructionException(
                            "Blank altlocs in duplicate residue %s ('%s', %i, '%s')"
                            % (resname, field, resseq, icode))
                    self.chain.detach_child(res_id)
                    new_residue = CustomResidue(res_id, resname, self.segid)
                    disordered_residue = DisorderedResidue(res_id)
                    self.chain.add(disordered_residue)
                    disordered_residue.disordered_add(duplicate_residue)
                    disordered_residue.disordered_add(new_residue)
                    self.residue = disordered_residue
                    return
        self.residue = CustomResidue(res_id, resname, self.segid)
        self.chain.add(self.residue)

    def init_atom(self, name, coord, b_factor, occupancy, altloc, fullname,
                  serial_number=None, element=None):
        """Create a new Atom object.
        Arguments:
         - name - string, atom name, e.g. CA, spaces should be stripped
         - coord - Numeric array (Float0, size 3), atomic coordinates
         - b_factor - float, B factor
         - occupancy - float
         - altloc - string, alternative location specifier
         - fullname - string, atom name including spaces, e.g. " CA "
         - element - string, upper case, e.g. "HG" for mercury
        """
        residue = self.residue
        # if residue is None, an exception was generated during
        # the construction of the residue
        if residue is None:
            return
        # First check if this atom is already present in the residue.
        # If it is, it might be due to the fact that the two atoms have atom
        # names that differ only in spaces (e.g. "CA.." and ".CA.",
        # where the dots are spaces). If that is so, use all spaces
        # in the atom name of the current atom.
        if residue.has_id(name):
            duplicate_atom = residue[name]
            # atom name with spaces of duplicate atom
            duplicate_fullname = duplicate_atom.get_fullname()
            if duplicate_fullname != fullname:
                # name of current atom now includes spaces
                name = fullname
                warnings.warn("Atom names %r and %r differ "
                              "only in spaces at line %i."
                              % (duplicate_fullname, fullname,
                                 self.line_counter),
                              PDBConstructionWarning)
        self.atom = CustomAtom(name, coord, b_factor, occupancy, altloc,
                               fullname, serial_number, element)
        if altloc != " ":
            # The atom is disordered
            if residue.has_id(name):
                # Residue already contains this atom
                duplicate_atom = residue[name]
                if duplicate_atom.is_disordered() == 2:
                    duplicate_atom.disordered_add(self.atom)
                else:
                    # This is an error in the PDB file:
                    # a disordered atom is found with a blank altloc
                    # Detach the duplicate atom, and put it in a
                    # DisorderedAtom object together with the current
                    # atom.
                    residue.detach_child(name)
                    disordered_atom = DisorderedAtom(name)
                    residue.add(disordered_atom)
                    disordered_atom.disordered_add(self.atom)
                    disordered_atom.disordered_add(duplicate_atom)
                    residue.flag_disordered()
                    warnings.warn("WARNING: disordered atom found "
                                  "with blank altloc before line %i.\n"
                                  % self.line_counter,
                                  PDBConstructionWarning)
            else:
                # The residue does not contain this disordered atom
                # so we create a new one.
                disordered_atom = DisorderedAtom(name)
                residue.add(disordered_atom)
                # Add the real atom to the disordered atom, and the
                # disordered atom to the residue
                disordered_atom.disordered_add(self.atom)
                residue.flag_disordered()
        else:
            # The atom is not disordered
            residue.add(self.atom)


def load_pdb(fileName, ident=None):
    """
    Read PDB from file into biopython structure object

    Parameters
    ----------
    fileName : str
        the path to the file
    ident : str (optional)
        the desired identity of the structure object

    Returns
    -------
    values : :class:`paramagpy.protein.CustomStructure`
        a structure object containing the atomic coordinates
    """
    if not ident:
        ident = fileName
    parser = PDBParser(structure_builder=CustomStructureBuilder())
    return parser.get_structure(ident, fileName)
