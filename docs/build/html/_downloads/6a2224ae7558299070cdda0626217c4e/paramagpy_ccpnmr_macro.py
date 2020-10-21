import os, subprocess

from memops.gui.Entry import Entry
from memops.gui.FileSelect import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.Label import Label
from memops.gui.PulldownList import PulldownList
from memops.gui.Button import Button
from memops.gui.LabelFrame import LabelFrame
from memops.gui.MessageReporter import showWarning, showOkCancel, showYesNo, showInfo

from ccpnmr.analysis.popups.BasePopup import BasePopup
from ccpnmr.analysis.core.AssignmentBasic import getResonanceResidue

from ccpnmr.analysis.core import ExperimentBasic


def paramagpyMACRO(argServer):
	popup = PCSPopup(argServer.parent)
	popup.open()

class PCSPopup(BasePopup):

	def __init__(self, parent, *args, **kw):
		BasePopup.__init__(self, parent, *args, **kw)

	def body(self, guiFrame):
	
		spectraFrame = LabelFrame(guiFrame, text='Spectra', grid=(0,0))
		
		Label(spectraFrame, text='Diamagnetic Spectrum:', grid=(0,0))
		self.diaSpecPulldown = PulldownList(
			spectraFrame, callback=self.changeDiaSpec, grid=(0,1))

		Label(spectraFrame, text='Paramagnetic Spectrum:', grid=(1,0))
		self.paraSpecPulldown = PulldownList(
			spectraFrame, callback=self.changeParaSpec, grid=(1,1))

		pcsIOFrame = LabelFrame(guiFrame, text='PCS', grid=(1,0))
		
		Label(pcsIOFrame, text='Write Experimental PCS:', grid=(0,0))
		self.expPCSEntry = Entry(pcsIOFrame, text='pcsexp.npc', width=32, grid=(0,1))
		Button(pcsIOFrame, command=self.writePCSFileChange,
			grid=(0,2), text='Choose File')

		Label(pcsIOFrame, text='Read Calculated PCS:', grid=(1,0))
		self.calPCSEntry = Entry(pcsIOFrame, text='pcscal.npc', width=32, grid=(1,1))
		Button(pcsIOFrame, command=self.readPCSFileChange,
			grid=(1,2), text='Choose File')

		scriptFrame = LabelFrame(guiFrame, text='Paramagpy Script', grid=(2,0))
		
		Label(scriptFrame, text='PCS fitting script:', grid=(0,0))
		self.scriptEntry = Entry(scriptFrame, text='paramagpy_fit_pcs.py', width=32, grid=(0,1))
		Button(scriptFrame, command=self.writeScriptFileChange,
			grid=(0,2), text='Choose File')

		commandFrame = LabelFrame(guiFrame, text='Commands', grid=(3,0))

		Button(commandFrame, command=self.writePCS, grid=(0,0), text='Write PCS')
		Button(commandFrame, command=self.fitTensor, grid=(0,1), text='Fit Tensor')
		Button(commandFrame, command=self.readPCS, grid=(0,2), text='Read PCS')
		
		self.updateSpecPulldown()
	

	def updateSpecPulldown(self):
		spectra = ExperimentBasic.getSpectra(self.project)
		names = []
		for spectrum in spectra:
			names.append("{}:{}".format(spectrum.experiment.name, spectrum.name))
		
		self.diaSpecPulldown.setup(names, spectra, 0)
		self.paraSpecPulldown.setup(names, spectra, 0)

		if len(spectra)>1:
			self.diaSpec = spectra[0]
			self.paraSpec = spectra[1]

	def changeDiaSpec(self, spectrum):
		self.diaSpec = spectrum

	def changeParaSpec(self, spectrum):
		self.paraSpec = spectrum


	def writePCSFileChange(self):
		fileTypes = [  FileType('PCS', ['*.npc']), FileType('All', ['*'])]
		fileSelectPopup = FileSelectPopup(self,
			file=os.path.basename(self.expPCSEntry.get()),
			directory=os.path.dirname(self.expPCSEntry.get()),
			file_types=fileTypes,
			title='Save PCS values to file', dismiss_text='Cancel',
			selected_file_must_exist=False, multiSelect=False,)

		self.expPCSEntry.set(fileSelectPopup.getFile())

	def readPCSFileChange(self):
		fileTypes = [  FileType('PCS', ['*.npc']), FileType('All', ['*'])]
		fileSelectPopup = FileSelectPopup(self,
			file=os.path.basename(self.calPCSEntry.get()),
			directory=os.path.dirname(self.calPCSEntry.get()),
			file_types=fileTypes,
			title='Save PCS values to file', dismiss_text='Cancel',
			selected_file_must_exist=False, multiSelect=False,)

		self.calPCSEntry.set(fileSelectPopup.getFile())

	def writeScriptFileChange(self):
		fileTypes = [  FileType('python', ['*.py']), FileType('All', ['*'])]
		fileSelectPopup = FileSelectPopup(self,
			file=os.path.basename(self.scriptEntry.get()),
			directory=os.path.dirname(self.scriptEntry.get()),
			file_types=fileTypes,
			title='Paramagpy script', dismiss_text='Cancel',
			selected_file_must_exist=False, multiSelect=False,)

		self.scriptEntry.set(fileSelectPopup.getFile())


	def writePCS(self):
		diaPeaks = self.diaSpec.getActivePeakList()
		paraPeaks = self.paraSpec.getActivePeakList()
		deltaPeaks = self.subtract_peak_lists(diaPeaks, paraPeaks)
		self.write_npc(deltaPeaks)

	def readPCS(self):
		print(self.diaSpec.experiment.name)
		print(self.paraSpec.experiment.name)

		calPCSName = os.path.abspath(self.calPCSEntry.get())
		if not os.path.isfile(calPCSName):
			showWarning("PCS file not found",
				"Cound not find backcalculated PCS values. Please ensure the following file exists:\n{}".format(calPCSName))
			return

		pcs_values = {}
		with open(calPCSName) as o:
			for line in o:
				try:
					sequence, atomName, pcsValue, error = line.split()
					key = int(sequence), atomName
					pcs_values[key] = float(pcsValue)
				except ValueError:
					print("Line ignored in calculated PCS file:\n{}".format(line))

		peakListDia = self.diaSpec.getActivePeakList()
		newPeakList = self.paraSpec.newPeakList(isSimulated=True)
		newPeakList.analysisPeakList.setSymbolColor("#FF0000")
		newPeakList.setDetails("Back-calculated PCS")

		for diaPeak in peakListDia.sortedPeaks():
			paraPeak = newPeakList.newPeak()
			paraPeak.setAnnotation(diaPeak.getAnnotation())
			for ddim, pdim in zip(diaPeak.sortedPeakDims(), paraPeak.sortedPeakDims()):
				pdim.setAnnotation(ddim.getAnnotation())
				if not ddim.dataDimRef:
					paraPeak.delete()
					break

				contrib = ddim.findFirstPeakDimContrib()
				if not contrib:
					paraPeak.delete()
					break
				resonance = contrib.getResonance()
				if not resonance:
					paraPeak.delete()
					break
				if not resonance.resonanceSet:
					paraPeak.delete()
					break

				residue = getResonanceResidue(resonance)
				sequence = residue.seqCode
				atomName = resonance.getAssignNames()[0]
				key = sequence, atomName

				if key in pcs_values:
					newShift = ddim.getValue() + pcs_values[key]
					if self.within_axis(pdim, newShift):
						pdim.setValue(newShift)
					else:
						paraPeak.delete()
						break
				else:
					paraPeak.delete()
					break


	def fitTensor(self):
		scriptName = os.path.abspath(self.scriptEntry.get())
		expPCSName = os.path.abspath(self.expPCSEntry.get())
		if not os.path.isfile(scriptName):
			showWarning("Script not found",
				"Could not find paramagpy script. Please ensure the following script exists:\n{}".format(scriptName))
			return

		if not os.access(scriptName, os.X_OK):
			showWarning("Invalid permissions",
				"The paramagpy script does not have executable permissions. Please run in the terminal:\nchmod +x {}".format(scriptName))
			return

		if not os.path.isfile(expPCSName):
			showWarning("PCS file not found",
				"Cound not find experimental PCS file.\nEnsure the following file exists:{}".format(expPCSName))
			return

		subprocess.call([scriptName, expPCSName])


	def write_npc(self, pcsValues):
		fileName = os.path.abspath(self.expPCSEntry.get())
		with open(fileName, 'w') as o:
			for key in sorted(pcsValues):
				seq, atom = key
				value = pcsValues[key]
				line = "{0:<4d}{1:<2}{2:15.8f}  0.0\n".format(seq, atom, value)
				o.write(line)
			o.close()
		showInfo("Wrote PCS", "Wrote experimental PCS to: {}".format(fileName))

	@staticmethod
	def extract_peaks(peakList):
		out = []
		for peak in peakList.sortedPeaks():
				tmp = {}
				for peakDim in peak.sortedPeakDims():
					if not peakDim.dataDimRef:
						continue

					contrib = peakDim.findFirstPeakDimContrib()
					if not contrib:
						continue
					resonance = contrib.getResonance()
					if not resonance:
						continue
					if not resonance.resonanceSet:
						continue

					residue = getResonanceResidue(resonance)
					sequence = residue.seqCode
					atomName = resonance.getAssignNames()[0]
					tmp[(sequence, atomName)] = peakDim
				out.append((peak, tmp))
		return out


	@staticmethod
	def within_axis(peakDim, value):
		dataDim = peakDim.getDataDim()
		dr = dataDim.findFirstDataDimRef()
		diff = dr.spectralWidthOrig - dr.spectralWidth
		low = dr.refValue - dr.spectralWidth + dr.spectralWidthOrig/2.
		high = dr.refValue + dr.spectralWidthOrig/2.
		return low < value < high


	def subtract_peak_lists(self, diamagPeaks, paramagPeaks):

		dpeaks = self.extract_peaks(diamagPeaks)
		ppeaks = self.extract_peaks(paramagPeaks)

		dvals = {}
		pvals = {}

		for peak, dic in dpeaks:
			dvals.update(dic)

		for peak, dic in ppeaks:
			pvals.update(dic)

		out = {}
		for key in set(dvals) & set(pvals):
			out[key] = pvals[key].value - dvals[key].value

		return out

