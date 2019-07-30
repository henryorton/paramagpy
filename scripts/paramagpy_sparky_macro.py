import re
import os
import subprocess
import sputil
import tkutil
import Tkinter
import sparky

BACK_CALC_PEAK_COLOUR = 'green'


def read_write_pcs_files(session):
	sputil.the_dialog(ReadWritePCSFiles, session).show_window(1)


class ReadWritePCSFiles(tkutil.Dialog):
	def __init__(self, session):
		self.session = session
		self.grid_widgets = []

		tkutil.Dialog.__init__(self, session.tk, 'Fit Delta-Chi Tensor to PCS')

		self.dia_menu = sputil.spectrum_menu(session, self.top, 'Diamagnetic Spectrum:  ')
		self.dia_menu.frame.pack(side = 'top', anchor = 'w')

		self.para_menu = sputil.spectrum_menu(session, self.top, 'Paramagnetic Spectrum: ')
		self.para_menu.frame.pack(side = 'top', anchor = 'w')

		self.npc_write_path = tkutil.file_field(self.top, 'PCS save file: ', '.npc', save=1)
		self.npc_write_path.frame.pack(side = 'top', anchor = 'e')
		self.npc_write_path.set("pcsexp.npc")

		self.npc_read_path = tkutil.file_field(self.top, 'PCS read file: ', '.npc')
		self.npc_read_path.frame.pack(side = 'top', anchor = 'e')
		self.npc_read_path.set("pcscalc.npc")

		self.message = Tkinter.Label(self.top)
		self.message.pack(side = 'top', anchor = 'w')
		self.message['text'] = "\n"

		br = tkutil.button_row(self.top,
			('Write PCS', self.write_npc_cb),
			('Read PCS', self.read_npc_cb),
			('Clear PCS', self.clear_peaks_cb),
			('Close', self.close_cb))
		br.frame.pack(side = 'top', anchor = 'w')


	@staticmethod
	def unpack_peaks(peakList):
		peaks = {}
		for peak in peakList:
			if peak.note == "PCS_backcalc":
				continue
			peaks[peak.assignment] = peak.frequency
		return peaks

	def write_npc_cb(self):
		fileName = os.path.abspath(self.npc_write_path.get())
		diaSpec = self.dia_menu.spectrum()
		paraSpec = self.para_menu.spectrum()

		if not all([diaSpec, paraSpec, fileName]):
			raise IOError("You need to select two spectra and a valid file name")

		dia = self.unpack_peaks(diaSpec.peak_list())
		para = self.unpack_peaks(paraSpec.peak_list())
		common = set(dia) & set(para)

		pcs = {}
		for assig in common:
			new_freq = []
			for d, p in zip(dia[assig], para[assig]):
				new_freq.append(p-d)
			pcs[assig] = tuple(new_freq)

		out = {}
		for assig in pcs:
			for (tmp, atom), value in zip(sputil.parse_assignment(assig), pcs[assig]):
				seq = int(re.findall(r"\d+", tmp)[0])
				out[(seq, atom)] = value

		with open(fileName, 'w') as f:
			for key in sorted(out):
				seq, atom = key
				line = "{0:<3d} {1:<3s} {2:6.3f} 0.000\n".format(
					seq, atom, out[key])
				f.write(line)

		self.message['text'] = "PCS values written to file: \n{}".format(fileName)

	def read_npc_cb(self):
		fileName = os.path.abspath(self.npc_read_path.get())
		diaSpec = self.dia_menu.spectrum()
		paraSpec = self.para_menu.spectrum()

		if not all([diaSpec, paraSpec, fileName]):
			raise IOError("You need to select two spectra and a valid file name")

		pcs = {}
		with open(self.npc_read_path.get()) as f:
			for line in f:
				try:
					seq, atom, value, error = line.split()
					pcs[(seq, atom)] = float(value)
				except ValueError:
					print("Line ignored while reading npc file:\n{}".format(line))

		for peak in diaSpec.peak_list():
			new_freq = []
			parsed_assig = sputil.parse_assignment(peak.assignment)
			for (tmp, atom), value in zip(parsed_assig, peak.frequency):
				seq = re.findall(r"\d+", tmp)[0]
				key = seq, atom
				if key in pcs:
					new_freq.append(value + pcs[key])

			if len(new_freq)!=paraSpec.dimension:
				continue

			for freq, (low, high) in zip(new_freq, zip(*paraSpec.region)):
				if not low < freq < high:
					break

			peak = paraSpec.place_peak(tuple(new_freq))
			for i, assig in enumerate(parsed_assig):
				peak.assign(i, assig[0], assig[1])

			peak.note = "PCS_backcalc"
			peak.show_assignment_label()
			peak.label.color = BACK_CALC_PEAK_COLOUR
			peak.color = BACK_CALC_PEAK_COLOUR

		self.message['text'] = "PCS values read from: \n{}".format(fileName)



	def clear_peaks_cb(self):
		paraSpec = self.para_menu.spectrum()

		if not paraSpec:
			raise IOError("You need to select a paramagnetic spectrum")

		selected_ornaments = self.session.selected_ornaments()
		self.session.unselect_all_ornaments()

		for peak in paraSpec.peak_list():
			if peak.note == "PCS_backcalc":
				peak.selected = 1
		self.session.command_characters("")

		for ornament in selected_ornaments:
			if sparky.object_exists(ornament):
				ornament.selected = 1

		self.message['text'] = "PCS peaks deleted from spectrum:\n{}".format(paraSpec.name)

		



