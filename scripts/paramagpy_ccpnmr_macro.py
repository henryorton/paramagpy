import re, os, subprocess

# Input the location of the paramagpy PCS fitting script here:
# Defaults to home directory
PARAMAGPY_SCRIPT_PATH = os.path.expanduser("~/") + "paramagpy_fit_pcs.py"


def paramagpyMACRO(argServer):
	experimentalPCSFileName, diamagPeaks, paramagPeaks = calculatePCS(argServer)
	fitted_tensor = fitTensor(argServer, experimentalPCSFileName)
	calculatedPCSFileName = 'pcscal.npc'
	fittedTensorFileName = 'pcs_tensor.txt'
	readPCS(argServer, diamagPeaks, paramagPeaks, calculatedPCSFileName)



def calculatePCS(argServer):

	def get_peak_dims(peak):
		return tuple(peak.sortedPeakDims())

	def get_peak_values(peak):
		return tuple(i.value for i in get_peak_dims(peak))

	def get_peak_names(peak):
		return tuple(i.annotation for i in get_peak_dims(peak))

	def get_peak_atoms(peak):
		out = []
		for peakDim in get_peak_dims(peak):
			contribs = list(peakDim.peakDimContribs)
			if len(contribs)>0:
				out.append(contribs[0].resonance.assignNames[0])
			else:
				out.append("") 
		return tuple(out)

	def extract_peaks(peakList):
		out = {}
		for peak in peakList.sortedPeaks():
			name = get_peak_names(peak)[0]
			if name:
				match = re.match(r"^(\d+)", name)
				if match:
					seq = int(match.group())
				else:
					continue
			else:
				continue
			out[seq] = peak
		return out

	def subtract_common_peaks(diaPeakList, paraPeakList):
		dia = extract_peaks(diaPeakList)
		para = extract_peaks(paraPeakList)
		common = set(dia) & set(para)
		out = {}
		for seq in common:
			diaPeak = dia[seq]
			paraPeak = para[seq]
			diaValues = get_peak_values(diaPeak)
			paraValues = get_peak_values(paraPeak)
			atoms = get_peak_atoms(diaPeak)
			deltas = [j-i for i, j in zip(diaValues, paraValues)]
			for atom, delta in zip(atoms, deltas):
				out[(seq, atom)] = delta
		return out


	def write_npc(pcsValues, fileName='pcsexp.npc'):
		with open(fileName, 'w') as o:
			for key in sorted(pcsValues):
				seq, atom = key
				value = pcsValues[key]
				line = "{0:<4d}{1:<2}{2:15.8f}  0.0\n".format(seq, atom, value)
				o.write(line)
			o.close()
		return fileName


	project = argServer.getProject()
	diamagPeaks = argServer.getPeakList()
	paramagPeaks = argServer.getPeakList()

	delta_peaks = subtract_common_peaks(diamagPeaks, paramagPeaks)
	fileName = os.path.abspath(write_npc(delta_peaks))
	argServer.showInfo("Wrote PCS to file:\n{}".format(fileName))
	return fileName, diamagPeaks, paramagPeaks





def fitTensor(argServer, experimentalPCSFileName):

	if not os.path.isfile(PARAMAGPY_SCRIPT_PATH):
		argServer.showWarning("Could not find paramagpy script. Please ensure the following script exists:\n{}".format(PARAMAGPY_SCRIPT_PATH))
		return

	if not os.access(PARAMAGPY_SCRIPT_PATH, os.X_OK):
		argServer.showWarning("The paramagpy script does not have executable permissions. Please run in the terminal:\nchmod +x {}".format(PARAMAGPY_SCRIPT_PATH))
		return

	subprocess.call([PARAMAGPY_SCRIPT_PATH, experimentalPCSFileName])




def readPCS(argServer, diamagPeaks, paramagPeaks, calculatedPCSFileName):

	def get_spectral_widths(spectrum):
		specDims = spectrum.sortedDataDims()
		out = []
		for dataDim in specDims:
			dr = dataDim.sortedDataDimRefs()[0]
			diff = dr.spectralWidthOrig - dr.spectralWidth
			low = dr.refValue - dr.spectralWidth + dr.spectralWidthOrig/2.
			high = dr.refValue + dr.spectralWidthOrig/2.
			out.append((low, high))
		return tuple(out)

	def get_peak_dims(peak):
		return tuple(peak.sortedPeakDims())

	def get_peak_names(peak):
		return tuple(i.annotation for i in get_peak_dims(peak))

	def get_peak_values(peak):
		return tuple(i.value for i in get_peak_dims(peak))

	def get_peak_atoms(peak):
		out = []
		for peakDim in get_peak_dims(peak):
			contribs = list(peakDim.peakDimContribs)
			if len(contribs)>0:
				out.append(contribs[0].resonance.assignNames[0])
			else:
				out.append("") 
		return tuple(out)

	def extract_peaks(peakList):
		out = {}
		for peak in peakList.sortedPeaks():
			name = get_peak_names(peak)[0]
			if name:
				match = re.match(r"^(\d+)", name)
				if match:
					seq = int(match.group())
				else:
					continue
			else:
				continue
			out[seq] = peak
		return out

	pcs_values = {}
	with open(calculatedPCSFileName) as o:
		for line in o:
			sequence, atomName, pcsValue, error = line.split()
			key = int(sequence), atomName
			pcs_values[key] = float(pcsValue)

	paraSpec = paramagPeaks.parent
	Hrange, Xrange = get_spectral_widths(paraSpec)
	newPeaksList = paraSpec.newPeakList(isSimulated=True)

	diaPeaks = extract_peaks(diamagPeaks)

	for seq in diaPeaks:
		diaPeak = diaPeaks[seq]

		newPeakName = get_peak_names(diaPeak)[0] + '_calc'

		Hnuc, Xnuc = get_peak_atoms(diaPeak)
		Hkey = (seq, Hnuc)
		Xkey = (seq, Xnuc)
		diaH, diaX = get_peak_values(diaPeak)

		if Hkey in pcs_values and Xkey in pcs_values:
			paraH = diaH + pcs_values[Hkey]
			paraX = diaX + pcs_values[Xkey]

			if Hrange[0]<paraH<Hrange[1] and Xrange[0]<paraX<Xrange[1]:
				peak = newPeaksList.newPeak()
				peak.setAnnotation(newPeakName)
				H, X = get_peak_dims(peak)
				H.setAnnotation(str(seq))
				H.setValue(paraH)
				X.setAnnotation(str(seq))
				X.setValue(paraX)
