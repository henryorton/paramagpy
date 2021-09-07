#!/home/u5376227/.virtualenvs/paramagpyenv/bin/python
import tkinter as tk
from tkinter import filedialog, messagebox, ttk, font
from tkinter.colorchooser import askcolor
import numpy as np
import os, sys, subprocess, shutil
from pprint import pprint
from collections import OrderedDict

import paramagpy
from paramagpy import metal, fit, protein, dataparse

def platform_settings():
	if sys.platform == 'darwin':
		settings = {
			'listboxFontSize': 14
		}
	elif sys.platform == 'linux':
		settings = {
			'listboxFontSize': 10
		}
	else:
		settings = {}
	return settings


def format_path(path, length):
	if len(path) > length:
		fmt = "...{{:<{}}}".format(length-3)
		return fmt.format(path[len(path)-length+3:])
	else:
		fmt = "{{:<{}}}".format(length)
		return fmt.format(path)


def unpack_ranges(string):
	out = []
	for i in string.split(','):
		if '-' in i:
			j, k = map(int, i.split('-'))
			if j>=k:
				raise ValueError
			out += range(j,k+1)
		else:
			out += [int(i)]
	return set(out)


class Tooltip:
	"""Notifications that contain information 
	and appear on mouse hover"""
	def __init__(self, widget, text, bg='#FFFFEA', pad=(5, 3, 5, 3),
		waittime=1000, wraplength=250):
		self.waittime = waittime
		self.wraplength = wraplength
		self.widget = widget
		self.text = text
		self.widget.bind("<Enter>", self.onEnter)
		self.widget.bind("<Leave>", self.onLeave)
		self.widget.bind("<ButtonPress>", self.onLeave)
		self.bg = bg
		self.pad = pad
		self.id = None
		self.tw = None

	def onEnter(self, event=None):
		self.schedule()

	def onLeave(self, event=None):
		self.unschedule()
		self.hide()

	def schedule(self):
		self.unschedule()
		self.id = self.widget.after(self.waittime, self.show)

	def unschedule(self):
		id_ = self.id
		self.id = None
		if id_:
			self.widget.after_cancel(id_)

	def show(self):
		def tip_pos_calculator(widget, label,
			tip_delta=(10, 5), pad=(5, 3, 5, 3)):
			w = widget
			s_width, s_height = w.winfo_screenwidth(), w.winfo_screenheight()
			width, height = (pad[0] + label.winfo_reqwidth() + pad[2],
							 pad[1] + label.winfo_reqheight() + pad[3])
			mouse_x, mouse_y = w.winfo_pointerxy()
			x1, y1 = mouse_x + tip_delta[0], mouse_y + tip_delta[1]
			x2, y2 = x1 + width, y1 + height
			x_delta = x2 - s_width
			if x_delta < 0:
				x_delta = 0
			y_delta = y2 - s_height
			if y_delta < 0:
				y_delta = 0
			offscreen = (x_delta, y_delta) != (0, 0)
			if offscreen:
				if x_delta:
					x1 = mouse_x - tip_delta[0] - width

				if y_delta:
					y1 = mouse_y - tip_delta[1] - height
			offscreen_again = y1 < 0
			if offscreen_again:
				y1 = 0
			return x1, y1

		bg = self.bg
		pad = self.pad
		widget = self.widget

		self.tw = tk.Toplevel(widget)
		self.tw.wm_overrideredirect(True)
		win = tk.Frame(self.tw,
					   background=bg,
					   borderwidth=0)
		label = tk.Label(win,
						  text=self.text,
						  justify=tk.LEFT,
						  background=bg,
						  relief=tk.SOLID,
						  borderwidth=0,
						  wraplength=self.wraplength)
		label.grid(padx=(pad[0], pad[2]),
				   pady=(pad[1], pad[3]),
				   sticky=tk.NSEW)
		win.grid()
		x, y = tip_pos_calculator(widget, label)
		self.tw.wm_geometry("+%d+%d" % (x, y))

	def hide(self):
		tw = self.tw
		if tw:
			tw.destroy()
		self.tw = None


class CustomTextDefaultEntry(tk.Entry):
	"""Entry widget that has greyed out default text"""
	def __init__(self, parent, default_text, 
		returnKey=None, mode='models', width=25):
		super().__init__(parent, width=width)
		self.default_text = default_text
		self.mode = mode
		self.insert(0, default_text)
		self.returnKey = returnKey
		if self.mode=='models':
			self.defaultColour = 'grey'
			self.config(fg = self.defaultColour)
		elif self.mode=='files':
			self.defaultColour = 'black'
			self.config(fg = self.defaultColour)
		elif self.mode=='atoms':
			self.defaultColour = 'black'
			self.config(fg = self.defaultColour)
		self.bind('<FocusIn>', self.on_entry_click)
		self.bind('<FocusOut>', self.on_focusout)
		if returnKey:
			self.bind('<Return>', lambda x: returnKey())

	def on_entry_click(self, event):
		s = self.get()
		if self.mode=='models':
			if s == self.default_text:
				self.clear()
				self.insert(0, '')
				self.config(fg = 'black')
		elif self.mode=='files':
			f, e = os.path.splitext(s)
			self.icursor(len(f))
			self.select_range(0, len(f))
		elif self.mode=='atoms':
			self.select_range(0, len(s))

	def on_focusout(self, event):
		if self.get() == '':
			self.insert(0, self.default_text)
			self.config(fg = self.defaultColour)
		if self.mode=='atoms':
			self.returnKey()

	def clear(self):
		self.delete(0, "end")


class NumericEntry(tk.Entry):
	"""Includes float validation and variable storage in tk.StringVar"""
	def __init__(self, root, parse, display, state='normal', 
		label=None, formatter="{:.3f}", onlyPositive=False, dtype=float):

		self.dtype = dtype
		self.label = label
		self.formatter = lambda x: formatter.format(x)
		self.textVar = tk.StringVar()
		self.floatVar = tk.DoubleVar()
		self.onlyPositive = onlyPositive

		if callable(parse):
			self.parser = parse
			self.display = display
		else:
			self.dummy_store(display)
			self.parser = lambda x: self.dummy_store(x)
			self.display = self.dummy_get

		vcmd = (root.register(self.validate),'%P')
		super().__init__(root, textvariable=self.textVar, validate='key', 
			validatecommand=vcmd, width=6, state=state)
		
		self.bind('<FocusIn>', lambda x: self.highlight())
		self.bind('<FocusOut>', lambda x: self.update())
		self.bind('<Return>', lambda x: [self.update(), self.highlight()])
		self.bind('<KP_Enter>', lambda x: [self.update(), self.highlight()])
		self.update()

	def update(self):
		var = self.formatter(self.display())
		self.textVar.set(var)

	def validate(self, value_if_allowed):
		if self.onlyPositive:
			initchars = ['+','']
			if self.dtype==float:
				initchars += ['.','+.']
		else:
			initchars = ['+','-','']
			if self.dtype==float:
				initchars += ['.','+.','-.']

		if value_if_allowed in initchars:
			return True
		
		try:
			val = self.dtype(value_if_allowed)
			if self.onlyPositive and val<0.0:
				raise ValueError
			self.parser(val)
			return True
		except ValueError:
			return False

		return False

	def dummy_store(self, value):
		self.floatVar.set(self.dtype(value))

	def dummy_get(self):
		return self.dtype(self.floatVar.get())

	def highlight(self):
		self.icursor(len(self.get()))
		self.select_range(0, 'end')

	def enable(self):
		self.config(state='normal')
		if self.label:
			self.label.config(fg='black')

	def disable(self):
		self.config(state='disabled')
		if self.label:
			self.label.config(fg='grey')


class Popup(tk.Toplevel):
	def __init__(self, parent, title=""):
		super().__init__(parent)
		self.parent = parent
		x = self.parent.winfo_rootx()+self.parent.winfo_width()//2
		y = self.parent.winfo_rooty()+self.parent.winfo_height()//2
		self.geometry("+{0:d}+{1:d}".format(x,y))
		self.transient(parent)
		self.title(title)
		self.protocol("WM_DELETE_WINDOW", self.death)
		self.grab_set()

	def last_words(self, *args):
		pass

	def death(self, *args):
		self.last_words()
		self.parent.focus_set()
		self.destroy()


class ProgressPopup(Popup):
	def __init__(self, parent, variable, text, num_iter=None, 
		mode='determinate', auto_close=True, title="Fitting progress"):
		super().__init__(parent, title)
		self.auto_close = auto_close
		if num_iter:
			self.points = num_iter
		else:
			self.points = 50
		self.last = 0.0
		self.incr = 1.0/self.points
		self.var = variable
		self.traceID = self.var.trace("w", self.change)
		self.label = tk.Label(self, text=text)
		self.label.pack()
		self.prog = ttk.Progressbar(self, maximum=self.points, 
			length=200, mode=mode)
		self.prog.pack()
		self.update()

	def change(self, *args):
		val = self.var.get()
		if val >= 1.0:
			if self.auto_close:
				self.death()
			else:
				self.update()
		elif abs(val-self.last) > self.incr:
			self.last = val
			self.prog['value'] = int(val*self.points)
			self.prog.update()

	def set_label(self, newtext):
		self.label.config(text=newtext)

	def last_words(self, *args):
		self.var.trace_vdelete("w", self.traceID)


class PlotTensorPopup(Popup):
	def __init__(self, parent):
		super().__init__(parent, "Plot tensor isosurface")
		self.pathWidth = 40
		self.params = {}

		self.curdir = self.parent.parent.loadData.currentDir()
		self.default_fileName = "isosurface"

		tk.Label(self, text='Location:').grid(row=0,column=0,sticky='E')
		self.lbl_dir = tk.Label(self, 
			text=format_path(self.curdir, self.pathWidth)+'/')
		self.lbl_dir.grid(row=0,column=1)
		ttk.Button(self, text="...", 
			command=self.change_dir).grid(row=0,column=2)

		tk.Label(self, text='Save Directory Name:').grid(
			row=1,column=0,sticky='E')
		self.saveName = tk.StringVar(value='isosurface')
		self.saveName.trace('w',self.saveNameChange)
		ttk.Entry(self, textvariable=self.saveName).grid(
			row=1,column=1,sticky='ew')

		tk.Label(self, text='PyMol Script File:').grid(
			row=2,column=0,sticky='E')
		self.pymol_file_label = tk.Label(self, text="")
		self.pymol_file_label.grid(row=2,column=1, sticky='W')

		tk.Label(self, text='Mesh File:').grid(row=3,column=0,sticky='E')
		self.mesh_file_label = tk.Label(self, text="")
		self.mesh_file_label.grid(row=3,column=1, sticky='W')

		tk.Label(self, text="Points per \u00c5:").grid(
			row=4,column=0,sticky='E')
		self.dens_input = NumericEntry(self, None, 2, formatter="{:d}", onlyPositive=True, dtype=int)
		self.dens_input.grid(row=4,column=1, sticky='W')
		self.params['dens'] = self.dens_input.floatVar

		tk.Label(self, text="Grid size \u00c5:").grid(
			row=5,column=0,sticky='E')
		self.size_input = NumericEntry(self, None, 40.0, formatter="{:.1f}",
			onlyPositive=True)
		self.size_input.grid(row=5,column=1, sticky='W')
		self.params['size'] = self.size_input.floatVar

		if self.parent.parent.frm_pdb.prot:
			state = 'normal'
			value = 1
		else:
			state = 'disabled'
			value = 0

		tk.Label(self, text="Include PDB structure:").grid(
			row=6,column=0,sticky='E')
		self.params['pdb'] = tk.IntVar(value=value)
		chk = ttk.Checkbutton(self, variable=self.params['pdb'], 
			state=state).grid(row=6,column=1, sticky='W')

		if self.parent.dtype=='PCS':
			clevel = 1.0
			clevelUnit = 'ppm'
		elif self.parent.dtype=='PRE':
			clevel = 50.0
			clevelUnit = '/s'
		tk.Label(self, text="Contour Level "+clevelUnit).grid(
			row=7,column=0,sticky='E')
		self.cont_input = NumericEntry(self, None, clevel, formatter="{:.2f}",
			onlyPositive=True)
		self.cont_input.grid(row=7,column=1, sticky='W')
		self.params['cont'] = self.cont_input.floatVar

		ttk.Button(self, text='Save', command=self.save).grid(row=8,column=0)
		ttk.Button(self, text='Save + Pymol', 
			command=self.save_and_call_pymol).grid(row=8,column=1)

		self.saveNameChange()
		self.update()

	def change_dir(self):
		d = filedialog.askdirectory(initialdir=self.curdir)
		if d:
			self.curdir = d
			self.lbl_dir.config(text=format_path(d, self.pathWidth))

	def saveNameChange(self, *args):
		pymolLabel = os.path.join(self.saveName.get(), self.get_pymol_file())
		meshLabel = os.path.join(self.saveName.get(), self.get_mesh_file())
		self.pymol_file_label.config(text=pymolLabel)
		self.mesh_file_label.config(text=meshLabel)

	def get_pymol_file(self):
		return "{}.pml".format(self.saveName.get())

	def get_mesh_file(self):
		return "{}.pml.ccp4".format(self.saveName.get())

	def get_pymol_file_path(self):
		return os.path.join(self.curdir, 
			self.saveName.get(), self.get_pymol_file())

	def get_mesh_file_path(self):
		return os.path.join(self.curdir, 
			self.saveName.get(), self.get_mesh_file())
		 
	def save(self):
		saveDir = os.path.join(self.curdir, self.saveName.get())
		if not os.path.exists(saveDir):
			os.makedirs(saveDir)
		pymolFile = self.get_pymol_file_path()
		meshFile = self.get_mesh_file_path()
		tensor = self.parent.tensor

		if self.params['pdb'].get():
			pdbPath = self.parent.parent.frm_pdb.prot.id
			pdbFile = os.path.basename(pdbPath)
			shutil.copyfile(pdbPath, os.path.join(saveDir, pdbFile))
		else:
			pdbFile = None

		mesh, bounds = tensor.make_mesh(
			density=int(self.params['dens'].get()), 
			size=self.params['size'].get())
		if self.parent.dtype=='PCS':
			value_mesh = tensor.pcs_mesh(mesh)
		elif self.parent.dtype=='PRE':
			rtype = self.parent.parent.loadData.rtype()
			usesbm=self.parent.parent.fopts.params['sbm'].get()
			usedsa=self.parent.parent.fopts.params['dsa'].get()
			value_mesh = tensor.pre_mesh(mesh, rtype=rtype, 
				dsa=usedsa, sbm=usesbm)
		tensor.write_isomap(value_mesh, bounds, fileName=meshFile)
		tensor.write_pymol_script(isoval=self.params['cont'].get(), 
			surfaceName=self.saveName.get(),scriptName=pymolFile, 
			meshName=os.path.basename(meshFile), 
			pdbFile=pdbFile)
		self.death()

	def save_and_call_pymol(self):
		self.save()
		subprocess.Popen(["pymol", "{}".format(self.get_pymol_file_path())])


class PlotCorrelationPopup(Popup):
	def __init__(self, parent):
		try:
			from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
			from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
			from matplotlib.figure import Figure
			from matplotlib import rcParams
		except ImportError:
			print("You need to install matplotlib module for this function")
			raise
		title = "Plot {} correlation".format(parent.dtype)
		super().__init__(parent.parent, title)
		self.parent = parent
		self.dtype = parent.dtype
		self.curdir = self.parent.loadData.currentDir()
		rcParams["savefig.directory"] = self.curdir
		# rcParams['savefig.dpi'] = 300
		self.params = {}

		self.frm_opt = tk.Frame(self)
		self.frm_opt.pack()
		self.frm_plt = tk.Frame(self)
		self.frm_plt.pack()

		self.frm_mff = MultipleFitFrame(self.parent.parent.parent, 
			state='plotting', root=self.frm_opt)
		self.frm_mff.grid(row=0,column=0, rowspan=4, sticky='EW')

		self.params['leg'] = tk.IntVar(value=0)
		chk = ttk.Checkbutton(self.frm_opt, text="Include Legend",
			variable=self.params['leg']).grid(row=0,column=2, sticky='W')

		self.params['inf'] = tk.IntVar(value=0)
		chk = ttk.Checkbutton(self.frm_opt, text="Include Tensor info.",
			variable=self.params['inf']).grid(row=1,column=2, sticky='W')

		self.params['ann'] = tk.IntVar(value=0)
		chk = ttk.Checkbutton(self.frm_opt, text="Include Atom info.",
			variable=self.params['ann']).grid(row=2,column=2, sticky='W')

		self.params['avg'] = tk.IntVar(value=0)
		chk = ttk.Checkbutton(self.frm_opt, text="Average Models",
			variable=self.params['avg']).grid(row=3,column=2, sticky='W')

		ttk.Button(self.frm_opt, text='\u2193   Re-plot   \u2193', 
			command=self.plot).grid(row=4,column=0,columnspan=4,sticky='EW')

		self.fig = Figure(figsize=(5, 5), dpi=100, tight_layout=True)
		self.axes = self.fig.add_subplot(111)
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.frm_plt)
		self.canvas.draw()
		self.canvas.get_tk_widget().pack()
		self.toolbar = NavigationToolbar2Tk(self.canvas, self.frm_plt)
		self.toolbar.update()

		self.update()
		self.frm_mff.update()
		self.plot()


	def last_words(self, *args):
		self.frm_mff.remove_traces()

	@staticmethod
	def format_coord(x, y):
		return "Exp: {0:5.2f}, Cal:{1:5.2f}".format(x, y)

	def plot(self):
		self.axes.clear()
		self.axes.set_xlabel("Experiment")
		self.axes.set_ylabel("Calculated")
		self.axes.format_coord = self.format_coord

		minig, maxig = None, None

		for tab in self.frm_mff.get_chosen_tabs():
			if tab.data is None:
				continue
			if self.params['avg'].get():
				data = tab.get_model_average()
			else:
				data = tab.data
			d = data[data['use'] & ~np.isnan(data['exp'])]
			mini = min([np.min(d[i]) for i in ['exp','cal']])
			maxi = max([np.max(d[i]) for i in ['exp','cal']])
			if minig is None:
				minig = mini
			if maxig is None:
				maxig = maxi
			if mini < minig:
				minig = mini
			if maxi < maxig:
				maxig = maxi
				
			self.axes.errorbar(d['exp'],d['cal'],xerr=d['err'],marker='o', 
				lw=0, elinewidth=1, ms=3, label=tab.name, color=tab.colour.colour)

			if self.params['ann'].get():
				for atom, exp, cal in d[['atm','exp','cal']]:
					_, mdl, chn, (_, seq, _), (atm, _) = atom.get_full_id()
					s = "{}{}".format(atm,seq)
					self.axes.annotate(s, xy=(exp, cal), fontsize=8)

			if self.params['leg'].get():
				self.axes.legend(bbox_to_anchor=(0.95,0.05), loc="lower right")

			if self.params['inf'].get():
				s = tab.tensorFit.tensor.info(comment=False)[:-1]
				self.axes.text(0.05, 0.95, s, family='monospace',fontsize=6.5,
					bbox={'facecolor':'white', 'pad':2},
					transform=self.axes.transAxes, horizontalalignment='left', 
					verticalalignment='top')

			scale = 1.1
			self.axes.plot([minig*scale,maxig*scale], 
						   [minig*scale,maxig*scale], '-k', lw=0.5, zorder=0)
			self.axes.set_xlim(minig*scale, maxig*scale)
			self.axes.set_ylim(minig*scale, maxig*scale)
			self.axes.set_aspect(1.0)

		self.canvas.draw()


class SaveDataPopup(Popup):
	def __init__(self, parent):
		title = "Save {} {}".format(parent.dtype, parent.parent.name)
		super().__init__(parent, title)
		self.params = {}

		self.curdir = self.parent.parent.loadData.currentDir()
		self.default_fileName = self.parent.parent.name.replace(" ","_").lower()

		self.varFileType = tk.StringVar(value=".npc")
		self.varAverage = tk.StringVar(value='average')

		tk.Label(self, text='File type:').grid(
			row=0,column=0,sticky='E')
		ttk.Radiobutton(self, text='.npc', variable=self.varFileType, 
			value='.npc').grid(row=0,column=1, sticky='W')
		ttk.Radiobutton(self, text='.csv', variable=self.varFileType, 
			value='.csv').grid(row=0,column=2, sticky='W')

		ttk.Separator(self, orient='horizontal').grid(
			row=1,column=0,columnspan=3,sticky='ew')

		tk.Label(self, text='Models:').grid(
			row=2,column=0,sticky='E')
		ttk.Radiobutton(self, text='Average', variable=self.varAverage, 
			value='average').grid(row=2,column=1, sticky='W')
		ttk.Radiobutton(self, text='Current', variable=self.varAverage, 
			value='current').grid(row=2,column=2, sticky='W')

		ttk.Button(self, text='Cancel', command=self.death).grid(row=3,column=0)
		ttk.Button(self, text='Save', command=self.save).grid(row=3,column=2)

		self.update()

	def save(self):
		kwargs = {'defaultextension':self.varFileType.get(), 
			'initialdir':self.curdir, 
			'initialfile':self.default_fileName}
		fileName = filedialog.asksaveasfilename(**kwargs)
		if not fileName:
			self.death()
			return

		if self.varAverage.get()=='current':
			tmp = self.parent.parent.data
			data = tmp[tmp['mdl']==self.parent.currentModel.get()]
		elif self.varAverage.get()=='average':
			data = self.parent.parent.get_model_average()

		outData = []

		if self.varFileType.get()=='.csv':
			header = ",".join([i[0] for i in self.parent.fields.values()])
			outData.append(header+'\n')
			for row in data:
				d = {}
				_, _, d['chn'], (_, d['seq'], _), (d['atm'], _) = row['atm'].get_full_id()
				d['res'] = row['atm'].parent.resname
				d['exp'] = row['exp']
				d['cal'] = row['cal']
				d['err'] = row['err']
				d['dev'] = abs(row['exp'] - row['cal'])
				if row['use']:
					d['use']= 'x'
				else:
					d['use'] = ' '
				line = self.parent.printRowFmt.format(**d)
				outData.append(line+'\n')

		elif self.varFileType.get()=='.npc':
			for row in data:
				_, _, chn, (_, seq, _), (atm, _) = row['atm'].get_full_id()
				line = "{0:<4d}{1:<5}{2:11.5f}{3:11.5f}\n".format(
					seq, atm, row['cal'], 0.0)
				outData.append(line)

		with open(fileName, 'w') as o:
			o.writelines(outData)

		self.death()


class SelectionPopup(Popup):
	"""A popup allowing atom and residue selection for experimental data"""
	def __init__(self, parent):
		title = "Atoms and residues for fitting"
		super().__init__(parent, title)
		tk.Label(self, text="Atom:").grid(row=0,column=0,sticky='W')
		self.atoms = {}
		for i, atom in enumerate(self.parent.default_atom_selection):
			if atom in parent.atom_selection:
				value = True
			else:
				value = False
			self.atoms[atom] = tk.BooleanVar(value=value)
			ttk.Checkbutton(self, text=atom,variable=self.atoms[atom]
				).grid(row=i+1,column=0, sticky='W')

		ttk.Separator(self, orient='vertical').grid(
			row=1,column=1,rowspan=len(self.atoms)+1,sticky='NS')

		tk.Label(self, text="Residue:").grid(row=0,column=2,columnspan=2,sticky='W')
		self.residues = {}
		for i, residue in enumerate(protein.standard_aa_names):
			if residue in parent.resi_selection:
				value = True
			else:
				value = False
			self.residues[residue] = tk.BooleanVar(value=value)
			ttk.Checkbutton(self, text=residue,variable=self.residues[residue]
				).grid(row=i%4+1,column=2+i//4, sticky='W',padx=4)

		ttk.Button(self, text='Cancel', command=self.death).grid(row=5,column=2)
		ttk.Button(self, text='Save', command=self.save).grid(row=5,column=5)
		self.update()

	def save(self):
		self.parent.atom_selection = set(
			[i for i,j in self.atoms.items() if j.get()])
		self.parent.resi_selection = set(
			[i for i,j in self.residues.items() if j.get()])
		self.parent.update_views(2)
		self.death()


class ErrorSimulationPopup(Popup):
	"""A popup allowing atom and residue selection for experimental data"""
	def __init__(self, parent):
		try:
			from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
			from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk
			from matplotlib.figure import Figure
			from matplotlib import rcParams
		except ImportError:
			print("You need to install matplotlib module for this function")
			raise
		title = "Uncertianty Analysis"
		super().__init__(parent, title)

		self.frm_left = tk.Frame(self)
		self.frm_opt = tk.Frame(self.frm_left)
		self.frm_disp = tk.Frame(self.frm_left)
		self.frm_left.pack(side='left')
		self.frm_opt.pack()
		self.frm_disp.pack()

		self.startTensor = self.parent.tensor
		self.errorTensor = metal.Metal()
		self.errorTensor.temperature = 0.0
		self.errorTensor.B0 = 0.0

		if self.parent.dtype in ['PCS','PRE','CCR']:
			self.varMulti = tk.BooleanVar(value=False)
			tk.Label(self.frm_opt, text='Fitting:').grid(
				row=0,column=0,sticky='W')
			ttk.Radiobutton(self.frm_opt, text='Use Single Fit', variable=self.varMulti, 
				value=False).grid(row=1,column=0, sticky='W')
			ttk.Radiobutton(self.frm_opt, text='Use Multiple Fit', variable=self.varMulti, 
				value=True).grid(row=2,column=0, sticky='W')

		ttk.Separator(self.frm_opt, orient='vertical').grid(
			row=0,column=1,rowspan=4,sticky='ns')

		tk.Label(self.frm_opt, text='Noise Source:').grid(
			row=0,column=2,sticky='W')
		
		self.var_error_method = tk.StringVar(value='mc')
		self.var_error_method.trace("w", self.error_method_changed)

		but_mdl = ttk.Radiobutton(self.frm_opt, text='Structure Models', 
			variable=self.var_error_method, 
			value='md')
		but_mdl.grid(row=1,column=2, sticky='W')
		ttk.Radiobutton(self.frm_opt, text='Experimental Uncertainty', 
			variable=self.var_error_method, 
			value='mc').grid(row=2,column=2, sticky='W')
		ttk.Radiobutton(self.frm_opt, text='Data Fraction', 
			variable=self.var_error_method, 
			value='bs').grid(row=3,column=2, sticky='W')

		ttk.Separator(self.frm_opt, orient='vertical').grid(
			row=0,column=4,rowspan=4,sticky='ns')

		tk.Label(self.frm_opt, text="Parameters:").grid(
			row=0,column=5,columnspan=2,sticky='W')
		lab_iter = tk.Label(self.frm_opt, text="Iterations:")
		lab_iter.grid(row=1,column=5,sticky='W')
		self.iterations = NumericEntry(self.frm_opt, None, 200, 
			formatter="{:d}", dtype=int, onlyPositive=True, label=lab_iter)
		self.iterations.grid(row=1,column=6, sticky='W')
		lab_frac = tk.Label(self.frm_opt, text="Sample Fraction:")
		lab_frac.grid(row=2,column=5,sticky='W')
		self.bootfrac = NumericEntry(self.frm_opt, None, 0.9, 
			formatter="{:.3f}", dtype=float, onlyPositive=True, label=lab_frac)
		self.bootfrac.grid(row=2,column=6, sticky='W')
		
		ttk.Separator(self.frm_opt, orient='horizontal').grid(
			row=4,column=0,columnspan=8,sticky='ew')

		ttk.Button(self.frm_opt, text='Run', 
			command=self.run).grid(row=5,column=0,columnspan=8,sticky='EW')

		self.textbox = tk.Text(self, state="disabled", 
			width=29, height=32)
		self.textbox.pack(side='right')

		self.fig = Figure(figsize=(5, 3), dpi=100)#, tight_layout=True)
		self.axes = self.fig.add_subplot(111, projection='hammer')
		self.canvas = FigureCanvasTkAgg(self.fig, master=self.frm_disp)
		self.canvas.draw()
		self.canvas.get_tk_widget().pack()
		self.toolbar = NavigationToolbar2Tk(self.canvas, self.frm_disp)
		self.toolbar.update()

		if not self.parent.fopts.frm_pdb.has_models():
			but_mdl.config(state='disabled')

		self.update()
		self.error_method_changed()
		self.plot()
		self.type_tensors()

	def error_method_changed(self, *args):
		method = self.var_error_method.get()
		if method=='md':
			self.iterations.disable()
			self.bootfrac.disable()
		elif method=='mc':
			self.iterations.enable()
			self.bootfrac.disable()
		elif method=='bs':
			self.iterations.enable()
			self.bootfrac.enable()

	def type_tensors(self, *args):
		s = "Fitted Tensor:\n"
		s += self.startTensor.info(comment=False)
		s += "\nError Tensor:\n"
		s += self.errorTensor.info(comment=False)[:-1]
		self.textbox.config(state='normal')
		self.textbox.delete(1.0,tk.END)
		self.textbox.insert(tk.END, s)
		self.textbox.config(state='disabled')

	def run(self, *args):
		if self.parent.dtype in ['PCS','PRE','CCR']:
			if self.varMulti.get():
				dataTabs = self.parent.parent.parent.parent.frm_fmul.get_chosen_tabs()
			else:
				dataTabs = [self.parent.parent]
		else:
			dataTabs = [self.parent.parent]

		method = self.var_error_method.get()
		iters = int(self.iterations.floatVar.get())
		bsfrac = self.bootfrac.floatVar.get()

		if self.parent.dtype=='PCS':
			all_metals, std_metal = self.parent.fopts.error_pcs(dataTabs, method, iters, 
				bsfrac=bsfrac)
		elif self.parent.dtype=='PRE':
			all_metals, std_metal = self.parent.fopts.error_pre(dataTabs, method, iters, 
				bsfrac=bsfrac)
		elif self.parent.dtype=='RDC':
			all_metals, std_metal = self.parent.fopts.error_rdc(dataTabs, method, iters, 
				bsfrac=bsfrac)
		elif self.parent.dtype=='CCR':
			all_metals, std_metal = self.parent.fopts.error_ccr(dataTabs, method, iters, 
				bsfrac=bsfrac)

		self.errorTensor = std_metal
		def transform(vector):
			x, y, z = vector
			theta = np.arctan2(y, x)
			phi = -np.arccos(z) + np.pi/2.
			return theta, phi

		spcoords = []
		for metal in all_metals:
			x, y, z = metal.rotationMatrix.T
			spcoords.append(tuple(map(transform, [x,y,z])))

		self.type_tensors()
		self.plot(zip(*spcoords))

	def plot(self, points=None):
		self.axes.clear()
		self.axes.set_xlabel("theta")
		self.axes.set_ylabel("phi")
		self.axes.grid()
		if points is not None:
			for data, col, label in zip(points, ['r','g','b'], ['x','y','z']):
				theta, phi = zip(*data)
				self.axes.scatter(theta, phi, s=0.4, c=col, label=label, zorder=10)
			self.axes.legend()
		self.canvas.draw()


class MorePopup(Popup):
	"""A popup containing lots of information about the tensor"""
	def __init__(self, parent):
		title = "More tensor information"
		tensor = parent.tensor
		super().__init__(parent, title)
		self.text = tk.Text(self, width=55, height=30, font="Courier 10")
		self.text.tag_configure("bold", font="Courier 10 bold")
		self.entry("Basic info.:", tensor.info(comment=False)[:-1])
		self.entry("Eigenvalues [dx,dy,dz]:", tensor.eigenvalues)
		self.entry("Delta Chi Tensor:", tensor.tensor_traceless)
		self.entry("Chi Tensor:", tensor.tensor)
		self.entry("Isotropy:", tensor.isotropy)
		self.entry("Rotation matrix:", tensor.rotationMatrix)
		self.entry("Saupe tensor:", tensor.tensor_saupe)

		self.text.pack(expand=True, fill='x')
		self.update()

	def entry(self, title, text):
		self.text.insert(tk.END, title+'\n', "bold")
		self.text.insert(tk.END, str(text)+'\n\n')


class DataLoad(tk.LabelFrame):
	"""
	For loading PCS, RDC and PRE data and parsing to <data>
	data : dictionary of {(sequence, atomName) : value}
	"""
	GAMMAH = 2*np.pi*42.576E6

	def __init__(self, parent):
		super().__init__(parent, text='Experimental Data')
		self.parent = parent # DataTab
		self.dtype = parent.dtype
		self.pathWidth = 20
		self.data = dataparse.DataContainer()
		self.fields = {}
		self.currentFile = None

		ttk.Button(self, text='Read {} Data'.format(self.dtype),
			command=self.load_data).grid(row=0,column=0)

		self.lbl_dataFile = tk.Label(self, text=" "*self.pathWidth)
		self.lbl_dataFile.grid(row=0,column=1,columnspan=3)

		def parse(value):
			mag = (2*np.pi * 1E6 * value)/self.GAMMAH
			self.parent.tensorStart.tensor.B0 = mag
			self.parent.tensorFit.tensor.B0 = mag
		def disp():
			mag = self.parent.tensorFit.tensor.B0
			return (mag*self.GAMMAH)/(2*np.pi * 1E6)
		self.fields['mag'] = NumericEntry(self, parse, 
			disp, onlyPositive=True, formatter="{:.1f}")
		self.fields['mag'].label = tk.Label(self, text='B0/MHz')
		self.fields['mag'].label.grid(row=1,column=0)
		self.fields['mag'].grid(row=1,column=1)

		def parse(value):
			self.parent.tensorStart.tensor.temperature = value
			self.parent.tensorFit.tensor.temperature = value
		def disp():
			return self.parent.tensorFit.tensor.temperature
		self.fields['tem'] = NumericEntry(self, parse, 
			disp, onlyPositive=True, formatter="{:.2f}")
		self.fields['tem'].label = tk.Label(self, text='Temp./K')
		self.fields['tem'].label.grid(row=1,column=2)
		self.fields['tem'].grid(row=1,column=3)

		if self.dtype=='PRE':
			tk.Label(self, text='Relaxation Type:').grid(row=2,column=0)
			self.rtype_var = tk.StringVar(value='r2')
			ttk.Radiobutton(self, text='R1', variable=self.rtype_var, 
				value='r1').grid(row=2,column=1, sticky='W')
			ttk.Radiobutton(self, text='R2', variable=self.rtype_var, 
				value='r2').grid(row=2,column=2, sticky='W')

		elif self.dtype in ['RDC','CCR']:
			self.show_pair_var = tk.BooleanVar(value=False)
			ttk.Checkbutton(self, text='Show Custom Atoms', variable=self.show_pair_var, 
				command=self.show_pair_change,
				).grid(row=2,column=0, columnspan=2, sticky='W')
			self.pair_entry1 = CustomTextDefaultEntry(self, "H",
				returnKey=self.show_pair_change, mode='atoms', width=4)
			self.pair_entry1.grid(row=2,column=2, columnspan=1, sticky='W')
			self.pair_entry2 = CustomTextDefaultEntry(self, "N",
				returnKey=self.show_pair_change, mode='atoms', width=4)
			self.pair_entry2.grid(row=2,column=3, columnspan=1, sticky='W')

	def rtype(self):
		if self.dtype=='PRE':
			return self.rtype_var.get()
		else:
			return None

	def show_pair_change(self):
		self.parent.update(2)

	def currentDir(self):
		if self.currentFile:
			return os.path.dirname(self.currentFile)
		else:
			return os.getcwd()

	def load_data(self):
		filedict = {'PCS':'.npc','RDC':'.rdc','PRE':'.pre','CCR':'.ccr'}
		fileName = filedialog.askopenfilename(
			title="Choose {} file".format(self.dtype),
			defaultextension=filedict[self.dtype],
			filetypes=[('{} file'.format(filedict[self.dtype]),
				filedict[self.dtype]),('All files','*.*')])
		if not fileName:
			return
		if fileName and self.dtype=='PCS':
			self.data = dataparse.read_pcs(fileName)
		elif fileName and self.dtype=='RDC':
			self.data = dataparse.read_rdc(fileName)
		elif fileName and self.dtype=='PRE':
			self.data = dataparse.read_pre(fileName)
		elif fileName and self.dtype=='CCR':
			self.data = dataparse.read_ccr(fileName)
		else:
			self.data = dataparse.DataContainer()
		self.lbl_dataFile.config(text=format_path(fileName, self.pathWidth))
		self.currentFile = fileName
		self.parent.update(2)


class Treeviewer(ttk.Treeview):
	def __init__(self, root, columns):
		self.fields = columns
		self.columns = list(self.fields.keys())
		self.data_reference = {}
		self.length = 0

		super().__init__(root, columns=self.columns, 
			show='headings', height=20)

		for col, (name, formatter, length) in columns.items():
			self.heading(col, text=name,
				command=lambda c=col: self.sort(c, True))
			pixels = font.Font().measure('X'*length)
			self.column(col, minwidth=0, width=pixels, stretch='NO')

		self.last_sort = None

	def clear(self):
		self.delete(*self.get_children())

	def set_data(self, data):
		data_ref = []
		for i, fd in enumerate(data):
			row = [self.fields[n][1].format(v) for n, v in fd.items()]
			data_ref.append(tuple(fd.values()))
			self.insert('', 'end', values=row, iid=str(i))
		self.data_reference = {}
		for column, values in zip(self.columns, zip(*data_ref)):
			self.data_reference[column] = np.array(values)
			self.length = len(values)

	def sort(self, column, reverse):
		idxs = np.argsort(self.data_reference[column])
		if reverse:
			idxs = idxs[::-1]

		for i, iid in enumerate(idxs):
			self.move(str(iid), '', i)

		self.last_sort = lambda c=column: self.sort(c, reverse)
		next_sort = lambda c=column: self.sort(c, not reverse)
		self.heading(column, command=next_sort)

	def sort_previous(self):
		self.last_sort()


class DataView(tk.LabelFrame):
	"""
	A listbox to view pdb and experimental data in rows
	"""
	def __init__(self, parent):
		super().__init__(parent, text='View Data')
		self.parent = parent # DataTab
		self.dtype = parent.dtype

		self.fmt_qfac = "Q factor: {:.4f}"
		self.fields = OrderedDict([
			['use',('?'    ,'{:^1}'  , 2)],
			['chn',('Chn.' ,'{:^3}'  , 3)],
			['seq',('Seq.' ,'{:^4}'  , 5)],
			['res',('Res.' ,'{:^5}'  , 5)],
			['atm',('Atom' ,'{:^4}'  , 5)],
			['cal',('Calc.','{:^.3f}', 7)],
			['exp',('Exp.' ,'{:^.3f}', 7)],
			['err',('Err.' ,'{:^.3f}', 7)],
			['dev',('Dev.' ,'{:^.3f}', 7)]
		])

		self.printRowFmt = """{use:1s},{chn:1s},{seq:<3d},{res:<3s},{atm:6s},{cal:7.3f},{exp:7.3f},{err:7.3f},{dev:7.3f}"""

		self.rowReference = []
		self.currentIndex = tk.IntVar(value=0)
		self.currentModel = tk.IntVar(value=0)
		self.currentIndex.trace("w",self.idx_change)

		self.lbl_qfac = tk.Label(self)
		self.set_qfac(np.nan)
		self.lbl_qfac.grid(row=0,column=0)

		ttk.Button(self, text="Plot", 
			command=self.plot).grid(row=0,column=1)

		ttk.Button(self, text="Save", 
			command=self.save).grid(row=0,column=2)

		self.hide_nan = tk.BooleanVar(value=False)
		self.hide_nan.trace('w', self.hide_nan_change)
		ttk.Checkbutton(self, text="Hide 'nan'",
			variable=self.hide_nan).grid(row=0,column=3)

		self.listBox = Treeviewer(self, self.fields)

		self.listBox.grid(row=1,column=0,columnspan=4,pady=10,padx=5)
		self.listBox.bind("<Double-Button-1>", self.set_position)
		self.listBox.bind("x", self.toggle_use)
		self.scrollBar = ttk.Scrollbar(self, orient=tk.VERTICAL)
		self.scrollBar.grid(row=1,column=3,sticky='NSE')
		self.scrollBar.config(command=self.listBox.yview)
		self.listBox.config(yscrollcommand=self.scrollBar.set)

		self.buttonFrame = tk.Frame(self)
		self.buttonFrame.grid(row=2,column=0,columnspan=3)
		tk.Label(self.buttonFrame, text='Change models:').grid(row=0,column=0)
		self.buttonR = ttk.Button(self.buttonFrame, text='<', 
			command=self.prev_model, state='disabled')
		self.buttonF = ttk.Button(self.buttonFrame, text='>', 
			command=self.next_model, state='disabled')
		self.lblModel = tk.Label(self.buttonFrame, 
			textvariable=self.currentModel)
		self.buttonR.grid(row=0,column=1)
		self.lblModel.grid(row=0,column=2,padx=5)
		self.buttonF.grid(row=0,column=3)

	@property
	def models(self):
		models = self.parent.parent.parent.parent.parent.frm_coords.frm_pdb.models
		return sorted(list(models))

	def idx_change(self, *args):
		idx = self.currentIndex.get()
		if idx < 0:
			self.currentModel.set(0)
		else:
			self.currentModel.set(self.models[idx])

	def set_current_model(self, model):
		idx = self.models.index(model)
		self.currentIndex.set(idx)
		self.update()

	def draw_view(self):
		self.rowReference = []
		if self.parent.data is None:
			return

		fd = self.fields.copy()
		curModel = self.currentModel.get()
		data = self.parent.data[self.parent.data['mdl']==curModel]

		dispData = []
		if self.dtype in ['PCS','PRE']:
			for row in data:
				if row['use']:
					fd['use'] = 'x'
				else:
					fd['use'] = ' '
				atom = row['atm']
				_, mdl,fd['chn'],(_,seq,_),(fd['atm'],_) = atom.get_full_id()
				fd['seq'] = int(seq)
				fd['res'] = atom.parent.resname
				fd['exp'] = row['exp']
				fd['cal'] = row['cal']
				fd['err'] = row['err']
				fd['dev'] = abs(row['exp'] - row['cal'])
				# line = [self.fields[n][1].format(v) for n, v in fd.items()]
				if self.hide_nan.get():
					if np.isnan(fd['exp']):
						continue
				dispData.append(fd.copy())
				self.rowReference.append(row)

		elif self.dtype in ['RDC','CCR']:
			for row in data:
				if row['use']:
					fd['use'] = 'x'
				else:
					fd['use'] = ' '
				atom1 = row['atm']
				atom2 = row['atx']
				_, mdl,chn1,(_,seq1,_),(atm1,_) = atom1.get_full_id()
				_, mdl,chn2,(_,seq2,_),(atm2,_) = atom2.get_full_id()
				res1 = atom1.parent.resname
				res2 = atom2.parent.resname
				fd['seq'] = "{}-{}".format(seq1, seq2)
				fd['res'] = "{}-{}".format(res1, res2)
				fd['atm'] = "{}-{}".format(atm1, atm2)
				fd['chn'] = "{}-{}".format(chn1, chn2)
				fd['exp'] = row['exp']
				fd['cal'] = row['cal']
				fd['err'] = row['err']
				fd['dev'] = abs(row['exp'] - row['cal'])
				dispData.append(fd.copy())
				self.rowReference.append(row)

		self.listBox.set_data(dispData)

	def set_qfac(self, value):
		self.lbl_qfac.config(text=self.fmt_qfac.format(value))

	def plot(self):
		PlotCorrelationPopup(self.parent)

	def save(self):
		SaveDataPopup(self)

	def hide_nan_change(self, *args):
		self.update(reselect=False, rescroll=False)

	def update(self, reselect=True, rescroll=True):
		scrollPos = self.listBox.yview()
		curFocus = self.listBox.focus()
		curSelect = self.listBox.selection()
		self.listBox.clear()
		numModels = len(self.models)

		if self.parent.data is None:
			self.currentIndex.set(-1)
		elif numModels <=1:
			self.currentIndex.set(0)
			self.buttonR.config(state='disabled')
			self.buttonF.config(state='disabled')
		elif self.currentIndex.get() <= 0:
			self.currentIndex.set(0)
			self.buttonR.config(state='disabled')
			self.buttonF.config(state='active')
		elif self.currentIndex.get() >= numModels-1:
			self.currentIndex.set(numModels-1)
			self.buttonR.config(state='active')
			self.buttonF.config(state='disabled')
		else:
			self.buttonR.config(state='active')
			self.buttonF.config(state='active')
		self.draw_view()
		if rescroll:
			self.listBox.yview_moveto(scrollPos[0])
		if curSelect and reselect:
			self.listBox.focus(curFocus)
			self.listBox.selection_set(curSelect)
		if self.listBox.last_sort:
			self.listBox.last_sort()

	def next_model(self):
		self.currentIndex.set(self.currentIndex.get()+1)
		self.update()

	def prev_model(self):
		self.currentIndex.set(self.currentIndex.get()-1)
		self.update()

	def set_position(self, *args):
		foc = self.listBox.focus()
		if not foc:
			return
		idx = int(foc)
		row = self.rowReference[idx]
		atom = row['atm']
		self.parent.tensorStart.tensor.position = atom.coord*1E-10
		self.parent.tensorStart.update()

	def toggle_use(self, *args):
		idxs = map(int, self.listBox.selection())
		for idx in idxs:
			row = self.rowReference[idx]
			if not np.isnan(row['exp']):
				data = self.parent.data
				data['use'][np.where(data['idx']==row['idx'])] = not row['use']
		self.update()


class TensorFrame(tk.LabelFrame):
	"""
	The initial tensor fields 
	"""
	tensorVariables = OrderedDict([
		('x','x/\u00c5'),
		('y','y/\u00c5'),
		('z','z/\u00c5'),
		('a','\u03b1/\u00b0'),
		('b','\u03b2/\u00b0'),
		('g','\u03b3/\u00b0'),
		('ax','\u0394\u03c7ax/10\u207b\u00b3\u00b2'),
		('rh','\u0394\u03c7rh/10\u207b\u00b3\u00b2'),
		('iso','\u03c7iso/10\u207b\u00b3\u00b2'),
		('ref','ref./ppm'),
		('taur','\u03c4r/ns'),
		('t1e','T1e/ps')])

	def __init__(self, parent, text):
		super().__init__(parent, text=text)
		self.parent = parent # DataTab
		self.dtype = parent.dtype
		self.fields = OrderedDict()
		self.tensor = metal.Metal()

		if 'Initial' in text:
			tk.Label(self, text="Lanthanide Template:").grid(row=0,
				column=0,columnspan=3)
			self.templateVar = tk.StringVar(value="Zero")
			self.templateVar.trace("w", self.template_change)
			lanths = list(metal.Metal.lanth_lib.keys())
			ttk.Combobox(self, textvariable=self.templateVar,
				values=lanths,width=4, 
				state="readonly").grid(row=0,column=3,columnspan=2)

		for i, name in enumerate(['x','y','z']):
			def wrap(axis):
				ax = axis
				def parse(value):
					self.tensor.position[ax] = value*1E-10
				def disp():
					return self.tensor.position[ax]*1E10
				return parse, disp
			parse, disp = wrap(i)
			self.fields[name] = NumericEntry(self, parse, disp)

		def parse(value):
			self.tensor.taur = value*1E-9
		def disp():
			return self.tensor.taur*1E9
		self.fields['taur'] = NumericEntry(self, parse, disp, onlyPositive=True)

		for i, name in enumerate(['a','b','g']):
			def wrap(axis):
				ax = axis
				def parse(value):
					self.tensor.eulers[ax] = value*(np.pi/180.)
				def disp():
					return self.tensor.eulers[ax]*(180./np.pi)
				return parse, disp
			parse, disp = wrap(i)
			self.fields[name] = NumericEntry(self, parse, disp)

		def parse(value):
			self.tensor.t1e = value*1E-13
		def disp():
			return self.tensor.t1e*1E13
		self.fields['t1e'] = NumericEntry(self, parse, disp, onlyPositive=True)

		for i, name in enumerate(['ax','rh']):
			def wrap(axis):
				ax = axis
				def parse(value):
					self.tensor.axrh[ax] = value*1E-32
				def disp():
					return self.tensor.axrh[ax]*1E32
				return parse, disp
			parse, disp = wrap(i)
			self.fields[name] = NumericEntry(self, parse, disp)

		def parse(value):
			self.tensor.isotropy = value*1E-32
		def disp():
			return self.tensor.isotropy*1E32
		self.fields['iso'] = NumericEntry(self, parse, disp, onlyPositive=True)

		def parse(value):
			self.tensor.shift = value
		def disp():
			return self.tensor.shift
		self.fields['ref'] = NumericEntry(self, parse, disp)

		i=0
		ofs = 1
		for name, field in self.fields.items():
			label = tk.Label(self, text=self.tensorVariables[name])
			label.grid(row=ofs+i%4, column=(i//4)*2)
			field.grid(row=ofs+i%4, column=(i//4)*2+1)
			i+=1
			field.label = label
			self.fields[name] = field

		ttk.Button(self, text='Copy', command=self.copy).grid(row=4+ofs, 
			column=0, columnspan=2, sticky='EW')
		ttk.Button(self, text='Paste', command=self.paste).grid(row=4+ofs, 
			column=2, columnspan=2, sticky='EW')
		ttk.Button(self, text='Set UTR', command=self.set_utr).grid(row=4+ofs, 
			column=4, columnspan=2, sticky='EW')

		if 'Fitted' in text:
			ttk.Button(self, text='More', command=self.more).grid(row=5+ofs, 
				column=0, columnspan=2, sticky='EW')
			ttk.Button(self, text='Error Sim.', command=self.error_sim).grid(row=5+ofs, 
				column=2, columnspan=2, sticky='EW')
			if self.dtype in ['PCS', 'PRE']:
				ttk.Button(self, text='Plot', command=self.plot).grid(row=5+ofs, 
					column=4, columnspan=2, sticky='EW')


		self.update()

	@property
	def fopts(self):
		return self.parent.parent.parent.frm_fopts

	def template_change(self, *args):
		lanth = self.templateVar.get()
		if self.dtype=='PRE':
			self.tensor.set_lanthanide(lanth, set_dchi=False)
		else:
			self.tensor.set_lanthanide(lanth, set_dchi=True)
		self.update()

	def update(self):
		for field in self.fields.values():
			field.update()

		pcs_fields = set(['x','y','z','a','b','g','ax','rh','ref'])
		rdc_fields = set(['a','b','g','ax','rh'])
		pre_fields = set(['x','y','z','a','b','g','ax','rh','iso','taur','t1e'])
		ccr_fields = set(['x','y','z','a','b','g','ax','rh','iso','taur'])

		fpars = self.fopts.params
		if self.dtype=='PCS':
			for var in self.fields:
				field = self.fields[var]
				if var in pcs_fields:
					field.enable()
				else:
					field.disable()

			if fpars['ref'].get():
				self.fields['ref'].label.config(fg='black')
			else:
				self.fields['ref'].label.config(fg='red')

			if fpars['pos'].get():
				for i in ['x','y','z']:
					self.fields[i].label.config(fg='black')
			else:
				for i in ['x','y','z']:
					self.fields[i].label.config(fg='red')

		elif self.dtype=='RDC':
			for var in self.fields:
				field = self.fields[var]
				if var in rdc_fields:
					field.enable()
				else:
					field.disable()

		elif self.dtype=='PRE':
			for var in self.fields:
				field = self.fields[var]
				if var in pre_fields:
					field.enable()
				else:
					field.disable()

			if fpars['pos'].get():
				for i in ['x','y','z']:
					self.fields[i].label.config(fg='black')
			else:
				for i in ['x','y','z']:
					self.fields[i].label.config(fg='red')

			if fpars['dchi'].get():
				for i in ['ax','rh','a','b','g']:
					self.fields[i].label.config(fg='black')
			else:
				for i in ['ax','rh','a','b','g']:
					self.fields[i].label.config(fg='red')

			if fpars['taur'].get():
				self.fields['taur'].label.config(fg='black')
			else:
				self.fields['taur'].label.config(fg='red')

			if fpars['iso'].get():
				self.fields['iso'].label.config(fg='black')
			else:
				self.fields['iso'].label.config(fg='red')

			if fpars['taue'].get():
				self.fields['t1e'].label.config(fg='black')
			else:
				self.fields['t1e'].label.config(fg='red')

		elif self.dtype=='CCR':
			for var in self.fields:
				field = self.fields[var]
				if var in ccr_fields:
					field.enable()
				else:
					field.disable()

			if fpars['pos'].get():
				for i in ['x','y','z']:
					self.fields[i].label.config(fg='black')
			else:
				for i in ['x','y','z']:
					self.fields[i].label.config(fg='red')

			if fpars['dchi'].get():
				for i in ['ax','rh','a','b','g']:
					self.fields[i].label.config(fg='black')
			else:
				for i in ['ax','rh','a','b','g']:
					self.fields[i].label.config(fg='red')

			if fpars['taur'].get():
				self.fields['taur'].label.config(fg='black')
			else:
				self.fields['taur'].label.config(fg='red')

			if fpars['iso'].get():
				self.fields['iso'].label.config(fg='black')
			else:
				self.fields['iso'].label.config(fg='red')

	def copy(self):
		self.clipboard_clear()
		self.clipboard_append(self.tensor.info(comment=False))
		self.parent.parent.parent.parent.copied_tensor = self.tensor.copy()

	def paste(self):
		copied_tensor = self.parent.parent.parent.parent.copied_tensor
		if copied_tensor:
			self.tensor = copied_tensor.copy()
		self.update()

	def plot(self):
		PlotTensorPopup(self)

	def error_sim(self):
		ErrorSimulationPopup(self)

	def set_utr(self):
		self.tensor.set_utr()
		self.update()

	def more(self):
		MorePopup(self)


class DataTab(tk.Frame):
	"""
	A tab of the DataNotebook
	contains data loading, viewing and tensor info
	dataView : spreadsheet of data
	"""
	def __init__(self, parent, name):
		super().__init__(parent)
		self.parent = parent # DataNotebook
		self.dtype = parent.dtype
		self.name = name
		self.data = None

		self.tensorStart = TensorFrame(self,'Initial Tensor')
		self.tensorFit = TensorFrame(self,'Fitted Tensor')
		self.loadData = DataLoad(self)
		self.viewData = DataView(self)

		self.viewData.grid(row=0,column=0,rowspan=10,sticky='NS')
		self.loadData.grid(row=0,column=1,sticky='EW')
		self.tensorStart.grid(row=1,column=1)
		self.tensorFit.grid(row=3,column=1)

		ttk.Button(self, text="\u2193  Fit Tensor  \u2193", 
			command= lambda: self.fopts.single_calc(self)
			).grid(row=2,column=1,sticky='EW')

		ttk.Button(self, text='Back-calculate {}'.format(self.dtype),
			command= self.back_calc).grid(row=4, column=1, sticky='EW')

		# A check button for keeping track of available data:
		self.hasData = tk.IntVar(value=0)
		
	def is_current(self):
		if self==self.parent.current_tab:
			return True
		else:
			return False

	def rtype(self):
		return self.loadData.rtype()

	def update(self, strength=0):
		if self.loadData.data:
			self.hasData.set(1)
		else:
			self.hasData.set(0)
		self.tensorStart.update()
		self.tensorFit.update()
		if strength>=2:
			self.parse_data()
		if strength>=1 and self.is_current():
			self.viewData.update()
		
	@property
	def frm_pdb(self):
		return self.parent.parent.parent.parent.frm_coords.frm_pdb

	@property
	def fopts(self):
		return self.parent.parent.frm_fopts

	@property
	def fmul(self):
		return self.parent.parent.frm_fmul

	def parse_data(self):
		data = self.parent.parent.parent.templates[self.dtype]
		atom_selection = self.frm_pdb.atom_selection
		resi_selection = self.frm_pdb.resi_selection
		if data is not None:
			self.data = data.copy()
		else:
			return
		if self.dtype in ['PCS','PRE']:
			expdata = self.loadData.data
			usedKeys = set([])
			for row in self.data:
				atom = row['atm']
				_, mdl, chn, (_,seq,_), (atm, _) = atom.get_full_id()
				res = atom.parent.resname
				key = seq, atm
				exp, err = expdata.get(key, (None, None))
				if exp is not None:
					usedKeys.add(key)
					row['exp'] = exp
					row['err'] = err
					if atom.element in atom_selection and res in resi_selection:
						row['use'] = True

			unused = set(self.loadData.data) - usedKeys
			if unused:
				message = "WARNING: Some atoms were not found in the PDB:{}"
				print(message.format(unused))

		elif self.dtype in ['RDC','CCR']:
			dtype = self.parent.parent.parent.structdtype[self.dtype]
			expdata = self.loadData.data
			if not expdata:
				if self.loadData.show_pair_var.get():
					a1 = self.loadData.pair_entry1.get()
					a2 = self.loadData.pair_entry2.get()
					for row in self.data:
						if row['atm'].name==a1:
							_, mdl, chn, seq, (atm, _) = row['atm'].get_full_id()
							row['atx'] = self.frm_pdb.prot[mdl][chn][seq].child_dict.get(a2, None)
					self.data = self.data[self.data['atx']!=None]
				else:
					self.data = np.array([], dtype=dtype)
				return
			pdata = self.frm_pdb.prot.parse(expdata)
			df = []
			for r in pdata:
				
				if r['mdl'] not in self.frm_pdb.models:
					continue

				res1 = r['atm'].parent.resname
				res2 = r['atx'].parent.resname
				row = (r['mdl'], True, r['atm'], r['atx'], np.nan, r['exp'], r['err'], r['idx'])

				if r['atm'].element in atom_selection and res1 in resi_selection:
					if r['atx'].element in atom_selection and res2 in resi_selection:
						df.append(row)
			self.data = np.array(df, dtype=dtype)

	def get_fit_data(self):
		return self.data[self.data['use']]

	def back_calc(self):
		if self.dtype=='PCS':
			poss = np.array([atom.coord for atom in self.data['atm']])*1E-10
			pcss = self.tensorFit.tensor.fast_pcs(poss)
			if self.fopts.params['racs'].get():
				csas = np.array([atom.csa for atom in self.data['atm']])
				racs = self.tensorFit.tensor.fast_racs(csas)
				pcss += racs
			if self.fopts.params['rads'].get():
				rads = self.tensorFit.tensor.fast_rads(poss)
				pcss += rads
			self.data['cal'] = pcss

		elif self.dtype=='RDC':
			vecs = []
			gams = []
			for row in self.data:
				vecs.append(row['atx'].position - row['atm'].position)
				gams.append(row['atx'].gamma * row['atm'].gamma)
			vecs = np.array(vecs)
			gams = np.array(gams)
			self.data['cal'] = self.tensorFit.tensor.fast_rdc(vecs, gams)

		elif self.dtype=='PRE':
			poss = np.array([atom.coord for atom in self.data['atm']])*1E-10
			gamm = np.array([atom.gamma for atom in self.data['atm']])
			if self.fopts.params['csa'].get(): 
				csas = np.array([atom.csa for atom in self.data['atm']])
			else:
				csas = 0.0
			pres = self.tensorFit.tensor.fast_pre(poss, gamm, self.rtype(), 
				dsa=self.fopts.params['dsa'].get(), 
				sbm=self.fopts.params['sbm'].get(), csaarray=csas)

			self.data['cal'] = pres

		elif self.dtype=='CCR':
			poss = []
			gams = []
			dsts = []
			for row in self.data:
				poss.append(row['atm'].position)
				gams.append(row['atm'].gamma)
				dsts.append(row['atx'].dipole_shift_tensor(row['atm'].position))
			poss = np.array(poss)
			gams = np.array(gams)
			dsts = np.array(dsts)
			self.data['cal'] = self.tensorFit.tensor.fast_ccr(poss, gams, dsts)

		self.set_qfactor()
		self.update(1)

	def set_qfactor(self):
		filt = self.data[self.data['use']]
		qfac = fit.qfactor(filt, ensembleAverage=self.fopts.params['eav'].get())
		self.viewData.set_qfac(qfac)
		return qfac

	def get_model_average(self):
		model = list(self.frm_pdb.models)[0]
		new = self.data[self.data['mdl']==model].copy()
		for row in new:
			calcs = self.data['cal'][self.data['idx']==row['idx']]
			row['mdl'] = -1
			row['cal'] = np.sum(calcs)/len(calcs)
		return new


class ColourSampler(tk.Canvas):
	"""docstring for ColourSampler"""
	class_cols = ['#0000ff','#ff0000','#00ff00','#bbbb00','#00bbbb']
	class_iter = [0]

	def __init__(self, parent, defaultColour=None):
		super().__init__(parent, width=40, height=20)
		self.parent = parent
		if not defaultColour:
			defaultColour = self.class_cols[self.class_iter[0]]
			self.class_iter[0] = (self.class_iter[0]+1)%5
		self.colour = defaultColour
		self.bind("<Button-1>", self.change_colour)
		self.box = self.create_rectangle(0, 0, 40, 20, fill=self.colour)

	def change_colour(self, *args):
		rgb, hexi = askcolor(initialcolor=self.colour)
		if hexi:
			self.colour = hexi
		self.itemconfig(self.box, fill=hexi)


class MultipleFitFrame(tk.LabelFrame):
	"""Multiple fit PCS frame"""
	def __init__(self, parent, state='fitting', root=None):
		if not root:
			root = parent
		if state=='fitting':
			label = "Multiple Fit {}".format(parent.dtype)
		elif state=='plotting':
			label = "Plot multiple {} data".format(parent.dtype)

		super().__init__(root, text=label)
		self.parent = parent # MethodsTab
		self.dtype = parent.dtype
		self.params = []

		if state=='fitting':
			num_dtabs = len(self.parent.ntb_data.tabs)
			ttk.Button(self, text="Multiple Fit Tensor", 
				command=self.multiple_calc).grid(row=2,columnspan=num_dtabs)
		for i, tab in enumerate(self.parent.ntb_data.tabs):
			useMulti = tk.IntVar(value=0)
			if state=='plotting' and tab.is_current() and tab.hasData.get():
				useMulti.set(1)
			chk = ttk.Checkbutton(self, text=tab.name, 
				variable=useMulti, state='disabled')
			chk.grid(row=0,column=i)
			chk.traceID = tab.hasData.trace("w", self.has_data_change)
			if state=='plotting':
				tab.colour = ColourSampler(self)
				tab.colour.grid(row=1, column=i)

			self.params.append((tab, chk, useMulti))
		self.update()

	@property
	def fopts(self):
		return self.parent.frm_fopts

	def has_data_change(self, *args):
		self.update()
		for tab, chk, var in self.params:
			if tab.hasData.get():
				var.set(1)

	def update(self, *args):
		for tab, chk, var in self.params:
			if tab.hasData.get():
				chk.config(state='normal')
			else:
				chk.config(state='disabled')

	def get_chosen_tabs(self):
		tabs = []
		for tab, chk, var in self.params:
			if var.get():
				tabs.append(tab)
		return tabs

	def multiple_calc(self):
		tabs = self.get_chosen_tabs()
		self.fopts.multiple_calc(tabs)

	def remove_traces(self):
		for tab, chk, var in self.params:
			tab.hasData.trace_vdelete("w", chk.traceID)


class FittingOptionsFrame(tk.LabelFrame):
	"""
	Fitting options and logical combinations
	Affects starting tensor values
	"""
	varNames = {
	'svd':'SVD Gridsearch',
	'nlr':'NLR Gradient Descent',
	'ref':'Fit Offset',
	'pos':'Fit Position',
	'rad':'SVD Radius/\u00c5',
	'den':'SVD Grid Spacing/\u00c5',
	'eav':'Ensemble Average',
	'racs':'Use RACS',
	'rads':'Use RADS',
	'taur':'Fit \u03c4r',
	'taue':'Fit T1e',
	'dsa':'Use DSA',
	'sbm':'Use SBM',
	'csa':'Use CSA',
	'dchi':'Fit \u0394\u03c7 tensor',
	'iso':'Fit \u03c7iso',
	}

	def __init__(self, parent):
		super().__init__(parent, text='Fitting Options')
		self.parent = parent # MethodsTab
		self.dtype = parent.dtype
		self.params = {}
		self.fields = {}

		if self.dtype=='PCS':
			self.set_checkbox('ref', 0, 0, 1, 0)
			self.set_checkbox('pos', 1, 0, 1)
			self.set_checkbox('eav', 2, 0, 1, 0)
			ttk.Separator(self,orient='vertical').grid(
				row=0,column=1,rowspan=3,sticky='ns')
			self.set_checkbox('svd', 0, 2, 2)
			self.set_field('rad', 1, 2, 1, 10.0)
			self.set_field('den', 2, 2, 1, 1.0)
			ttk.Separator(self,orient='vertical').grid(
				row=0,column=4,rowspan=3,sticky='ns')
			self.set_checkbox('nlr', 0, 5, 1)
			self.set_checkbox('racs', 1, 5, 1, 0)
			self.set_checkbox('rads', 2, 5, 1, 0)

		elif self.dtype=='RDC':
			self.set_checkbox('eav', 0, 0, 1, 0)

		elif self.dtype=='PRE':
			self.set_checkbox('iso', 0, 0, 1, 0)
			self.set_checkbox('pos', 1, 0, 1)
			self.set_checkbox('eav', 2, 0, 1, 0)

			self.set_checkbox('dchi', 0, 1, 1, 0)
			self.set_checkbox('taur', 1, 1, 1, 0)
			self.set_checkbox('taue', 2, 1, 1, 0)

			self.set_checkbox('dsa', 0, 2, 1)
			self.set_checkbox('sbm', 1, 2, 1)
			self.set_checkbox('csa', 2, 2, 1, 0)

		elif self.dtype=='CCR':
			self.set_checkbox('iso', 0, 0, 1, 0)
			self.set_checkbox('pos', 1, 0, 1)
			self.set_checkbox('eav', 2, 0, 1, 0)

			self.set_checkbox('dchi', 0, 1, 1, 0)
			self.set_checkbox('taur', 1, 1, 1, 0)


	@property
	def frm_pdb(self):
		return self.parent.parent.parent.frm_coords.frm_pdb

	@property
	def tensorStart(self):
		return self.parent.ntb_data.current_tab.tensorStart

	def set_checkbox(self, variableName, row, col, colspan, state=1):
		var = tk.IntVar(value=state)
		chk = ttk.Checkbutton(self, text=self.varNames[variableName], 
			variable=var, command=self.update)
		Tooltip(chk, tt[variableName])
		chk.grid(row=row, column=col, columnspan=colspan, sticky='W')
		self.params[variableName] = var
		self.fields[variableName] = chk
		return chk

	def set_field(self, variableName, row, col, colspan, default):
		field = NumericEntry(self, None, default, onlyPositive=True)
		label = tk.Label(self, text=self.varNames[variableName])
		Tooltip(field, tt[variableName])
		Tooltip(label, tt[variableName])
		field.label = label
		self.fields[variableName] = field
		self.params[variableName] = field.floatVar
		label.grid(row=row,column=col,columnspan=colspan,sticky='W')
		field.grid(row=row,column=col+1,columnspan=colspan)

	def update(self, strength=0):
		if self.dtype=='PCS':
			if not self.params['svd'].get():
				self.fields['rad'].disable()
				self.fields['den'].disable()
			elif not self.params['pos'].get():
				self.fields['rad'].floatVar.set(0.0)
				self.fields['rad'].update()
				self.fields['rad'].disable()
				self.fields['den'].disable()
			else:
				self.fields['rad'].enable()
				self.fields['den'].enable()

		self.parent.update(0)

	def single_calc(self, dataTab):
		if self.dtype=='PCS':
			self.fit_pcs([dataTab])
		if self.dtype=='RDC':
			self.fit_rdc(dataTab)
		if self.dtype=='PRE':
			self.fit_pre([dataTab])
		if self.dtype=='CCR':
			self.fit_ccr([dataTab])

	def multiple_calc(self, dataTabs):
		if self.dtype=='PCS':
			self.fit_pcs(dataTabs)
		if self.dtype=='RDC':
			raise TypeError("RDC does not support multiple fit")
		if self.dtype=='PRE':
			self.fit_pre(dataTabs)
		if self.dtype=='CCR':
			self.fit_ccr(dataTabs)

	def get_data(self, dataTabs):
		return [tab.get_fit_data() for tab in dataTabs]

	def get_params(self):
		if self.dtype=='PCS':
			pars = ['ax','rh','a','b','g']
			if self.params['pos'].get():
				pars += ['x','y','z']
			if self.params['ref'].get():
				pars += ['shift']
			return pars
		elif self.dtype=='PRE':
			pars = []
			if self.params['pos'].get():
				pars += ['x','y','z']
			if self.params['iso'].get():
				pars += ['iso']
			if self.params['dchi'].get():
				pars += ['ax','rh','a','b','g']
			if self.params['taur'].get():
				pars += ['taur']
			if self.params['taue'].get():
				pars += ['t1e']
			return pars
		elif self.dtype=='CCR':
			pars = []
			if self.params['pos'].get():
				pars += ['x','y','z']
			if self.params['iso'].get():
				pars += ['iso']
			if self.params['dchi'].get():
				pars += ['ax','rh','a','b','g']
			if self.params['taur'].get():
				pars += ['taur']
			return pars


	def check_data(self, dataTabs):
		if self.frm_pdb.prot is None:
			messagebox.showerror("Error", "PDB coordinates missing")
			return False
		for tab in dataTabs:
			if tab.loadData.currentFile is None:
				messagebox.showerror("Error", "Experimental data missing")
				return False
		return True


	def fit_pcs(self, dataTabs):
		if not self.check_data(dataTabs):
			return

		metals = [tab.tensorStart.tensor for tab in dataTabs]
		datas = self.get_data(dataTabs)

		progVar = tk.DoubleVar(value=0.0)
		progbar = ProgressPopup(self, progVar, "", auto_close=False)

		svdline = "SVD Gridsearch . . ."
		nlrline = "NLR Fitting . . ."

		if self.params['svd'].get():
			progbar.set_label(svdline)
			progVar.set(0.0)
			radius = float(self.fields['rad'].get())
			points = int(radius/float(self.fields['den'].get()))
			if points<1:
				points = 1
			ref = self.params['ref'].get()
			metals, calcs = fit.svd_gridsearch_fit_metal_from_pcs(
				metals, datas, 
				radius=radius, 
				points=points, 
				offsetShift=ref, 
				ensembleAverage=self.params['eav'].get(),
				progress=progVar)

		if self.params['nlr'].get():
			progbar.set_label(nlrline)
			progVar.set(0.0)
			pars = self.get_params()
			metals, calcs = fit.nlr_fit_metal_from_pcs(
				metals, datas, pars,
				userads=self.params['rads'].get(), 
				useracs=self.params['racs'].get(),
				ensembleAverage=self.params['eav'].get(),
				progress=progVar)

		progbar.death()

		for tab, metal in zip(dataTabs, metals):
			tab.tensorFit.tensor = metal.copy()
			tab.update(0)
			tab.back_calc()

		if self.frm_pdb.has_models():
			qs = {mdl:0.0 for mdl in self.frm_pdb.models}
			for mdl in qs:
				for calc in calcs:
					filt = calc[calc['use']]
					filt = filt[filt['mdl']==mdl]
					q = fit.qfactor(filt, ensembleAverage=self.params['eav'].get())
					qs[mdl] += q / len(calcs)
			qsort = sorted(qs.items(), key=lambda x: x[1])
			line = "Model {0:} found with minimum Q-factor of {1:5.3f}"
			messagebox.showinfo("Model with best fit found", 
				line.format(qsort[0][0], qsort[0][1]))
			tab.viewData.set_current_model(qsort[0][0])

	def error_pcs(self, dataTabs, method, iters=None, bsfrac=None):
		if not self.check_data(dataTabs):
			return

		if method=='md':
			title = "Structure Sourced Uncertainty Analysis"
		if method=='mc':
			title = "Experimental Error Sourced Uncertainty Analysis"
		elif method=='bs':
			title = "Sample Fractions Sourced Uncertainty Analysis"

		progVar = tk.DoubleVar(value=0.0)
		progbar = ProgressPopup(self, progVar, 
			"Repeating NLR fit . . .",
			title=title, auto_close=False)

		kwargs = {
			'initMetals':[tab.tensorStart.tensor for tab in dataTabs],
			'dataArrays':self.get_data(dataTabs),
			'params':self.get_params(),
			'ensembleAverage':self.params['eav'].get(),
			'userads':self.params['rads'].get(), 
			'useracs':self.params['racs'].get(), 
			'progress':progVar
		}

		if method=='md':
			all_metals, std_metals = fit.fit_error_models(
				fit.nlr_fit_metal_from_pcs, **kwargs)

		elif method=='mc':
			all_metals, std_metals = fit.fit_error_monte_carlo(
				fit.nlr_fit_metal_from_pcs, iters, **kwargs)

		elif method=='bs':
			all_metals, std_metals = fit.fit_error_bootstrap(
				fit.nlr_fit_metal_from_pcs, iters, bsfrac, **kwargs)

		progbar.death()

		for tab, am, sm in zip(dataTabs, all_metals, std_metals):
			if tab.is_current():
				return am, sm


	def fit_pre(self, dataTabs):
		if not self.check_data(dataTabs):
			return

		metals = [tab.tensorStart.tensor for tab in dataTabs]
		datas = self.get_data(dataTabs)

		progVar = tk.DoubleVar(value=0.0)
		progbar = ProgressPopup(self, progVar, "", auto_close=False)

		nlrline = "NLR Fitting . . ."

		progbar.set_label(nlrline)
		progVar.set(0.0)
		metals, calcs = fit.nlr_fit_metal_from_pre(
			metals, datas,
			rtypes=[tab.rtype() for tab in dataTabs],
			params=self.get_params(),
			usesbm=self.params['sbm'].get(), 
			usedsa=self.params['dsa'].get(),
			usecsa=self.params['csa'].get(),
			ensembleAverage=self.params['eav'].get(),
			progress=progVar)

		progbar.death()

		for tab, metal in zip(dataTabs, metals):
			tab.tensorFit.tensor = metal.copy()
			tab.update(0)
			tab.back_calc()

		if self.frm_pdb.has_models():
			qs = {mdl:0.0 for mdl in self.frm_pdb.models}
			for mdl in qs:
				for calc in calcs:
					filt = calc[calc['use']]
					filt = filt[filt['mdl']==mdl]
					q = fit.qfactor(filt, ensembleAverage=self.params['eav'].get())
					qs[mdl] += q / len(calcs)
			qsort = sorted(qs.items(), key=lambda x: x[1])
			line = "Model {0:} found with minimum Q-factor of {1:5.3f}"
			messagebox.showinfo("Model with best fit found", 
				line.format(qsort[0][0], qsort[0][1]))
			tab.viewData.set_current_model(qsort[0][0])


	def error_pre(self, dataTabs, method, iters=None, bsfrac=None):
		if not self.check_data(dataTabs):
			return

		if method=='md':
			title = "Structure Sourced Uncertainty Analysis"
		if method=='mc':
			title = "Experimental Error Sourced Uncertainty Analysis"
		elif method=='bs':
			title = "Sample Fractions Sourced Uncertainty Analysis"

		progVar = tk.DoubleVar(value=0.0)
		progbar = ProgressPopup(self, progVar, 
			"Repeating NLR fit . . .",
			title=title, auto_close=False)

		kwargs = {
			'initMetals':[tab.tensorStart.tensor for tab in dataTabs],
			'dataArrays':self.get_data(dataTabs),
			'rtypes':[tab.rtype() for tab in dataTabs],
			'params':self.get_params(),
			'ensembleAverage':self.params['eav'].get(),
			'usesbm':self.params['sbm'].get(),
			'usedsa':self.params['dsa'].get(),
			'usecsa':self.params['csa'].get(),
			'progress':progVar
		}

		if method=='md':
			all_metals, std_metals = fit.fit_error_models(
				fit.nlr_fit_metal_from_pre, **kwargs)

		elif method=='mc':
			all_metals, std_metals = fit.fit_error_monte_carlo(
				fit.nlr_fit_metal_from_pre, iters, **kwargs)

		elif method=='bs':
			all_metals, std_metals = fit.fit_error_bootstrap(
				fit.nlr_fit_metal_from_pre, iters, bsfrac, **kwargs)

		progbar.death()

		for tab, am, sm in zip(dataTabs, all_metals, std_metals):
			if tab.is_current():
				return am, sm


	def fit_rdc(self, dataTab):
		if not self.check_data([dataTab]):
			return

		metal = dataTab.tensorStart.tensor.copy()
		[data] = self.get_data([dataTab])

		progVar = tk.DoubleVar(value=0.0)
		progbar = ProgressPopup(self, progVar, "", auto_close=False)

		nlrline = "SVD Fitting . . ."

		progbar.set_label(nlrline)
		progVar.set(0.0)
		[metal], [calc] = fit.svd_fit_metal_from_rdc(
			[metal], [data],
			ensembleAverage=self.params['eav'].get(),
			progress=progVar)


		dataTab.tensorFit.tensor = metal.copy()
		dataTab.update(0)
		dataTab.back_calc()
		progbar.death()

		if self.frm_pdb.has_models():
			qs = {}
			for mdl in self.frm_pdb.models:
				filt = calc[calc['use']]
				filt = filt[filt['mdl']==mdl]
				qs[mdl] = fit.qfactor(filt, ensembleAverage=self.params['eav'].get())
			qsort = sorted(qs.items(), key=lambda x: x[1])
			line = "Model {0:} found with minimum Q-factor of {1:5.3f}"
			messagebox.showinfo("Model with best fit found", 
				line.format(qsort[0][0], qsort[0][1]))
			dataTab.viewData.set_current_model(qsort[0][0])


	def error_rdc(self, dataTabs, method, iters=None, bsfrac=None):
		if not self.check_data(dataTabs):
			return

		if method=='md':
			title = "Structure Sourced Uncertainty Analysis"
		if method=='mc':
			title = "Experimental Error Sourced Uncertainty Analysis"
		elif method=='bs':
			title = "Sample Fractions Sourced Uncertainty Analysis"

		progVar = tk.DoubleVar(value=0.0)
		progbar = ProgressPopup(self, progVar, 
			"Repeating SVD calc . . .",
			title=title, auto_close=False)

		kwargs = {
			'initMetals':[tab.tensorStart.tensor for tab in dataTabs],
			'dataArrays':self.get_data(dataTabs),
			'ensembleAverage':self.params['eav'].get(),
			'progress':progVar
		}

		if method=='md':
			all_metals, std_metals = fit.fit_error_models(
				fit.svd_fit_metal_from_rdc, **kwargs)

		elif method=='mc':
			all_metals, std_metals = fit.fit_error_monte_carlo(
				fit.svd_fit_metal_from_rdc, iters, **kwargs)

		elif method=='bs':
			all_metals, std_metals = fit.fit_error_bootstrap(
				fit.svd_fit_metal_from_rdc, iters, bsfrac, **kwargs)

		progbar.death()

		for tab, am, sm in zip(dataTabs, all_metals, std_metals):
			if tab.is_current():
				return am, sm


	def fit_ccr(self, dataTabs):
		if not self.check_data(dataTabs):
			return

		metals = [tab.tensorStart.tensor for tab in dataTabs]
		datas = self.get_data(dataTabs)

		progVar = tk.DoubleVar(value=0.0)
		progbar = ProgressPopup(self, progVar, "", auto_close=False)

		nlrline = "NLR Fitting . . ."

		progbar.set_label(nlrline)
		progVar.set(0.0)
		metals, calcs = fit.nlr_fit_metal_from_ccr(
			metals, datas,
			params=self.get_params(),
			ensembleAverage=self.params['eav'].get(),
			progress=progVar)

		progbar.death()

		for tab, metal in zip(dataTabs, metals):
			tab.tensorFit.tensor = metal.copy()
			tab.update(0)
			tab.back_calc()

		if self.frm_pdb.has_models():
			qs = {mdl:0.0 for mdl in self.frm_pdb.models}
			for mdl in qs:
				for calc in calcs:
					filt = calc[calc['use']]
					filt = filt[filt['mdl']==mdl]
					q = fit.qfactor(filt, ensembleAverage=self.params['eav'].get())
					qs[mdl] += q / len(calcs)
			qsort = sorted(qs.items(), key=lambda x: x[1])
			line = "Model {0:} found with minimum Q-factor of {1:5.3f}"
			messagebox.showinfo("Model with best fit found", 
				line.format(qsort[0][0], qsort[0][1]))
			tab.viewData.set_current_model(qsort[0][0])


	def error_ccr(self, dataTabs, method, iters=None, bsfrac=None):
		if not self.check_data(dataTabs):
			return

		if method=='md':
			title = "Structure Sourced Uncertainty Analysis"
		if method=='mc':
			title = "Experimental Error Sourced Uncertainty Analysis"
		elif method=='bs':
			title = "Sample Fractions Sourced Uncertainty Analysis"

		progVar = tk.DoubleVar(value=0.0)
		progbar = ProgressPopup(self, progVar, 
			"Repeating NLR fit . . .",
			title=title, auto_close=False)

		kwargs = {
			'initMetals':[tab.tensorStart.tensor for tab in dataTabs],
			'dataArrays':self.get_data(dataTabs),
			'params':self.get_params(),
			'ensembleAverage':self.params['eav'].get(),
			'progress':progVar
		}

		if method=='md':
			all_metals, std_metals = fit.fit_error_models(
				fit.nlr_fit_metal_from_ccr, **kwargs)

		elif method=='mc':
			all_metals, std_metals = fit.fit_error_monte_carlo(
				fit.nlr_fit_metal_from_ccr, iters, **kwargs)

		elif method=='bs':
			all_metals, std_metals = fit.fit_error_bootstrap(
				fit.nlr_fit_metal_from_ccr, iters, bsfrac, **kwargs)

		progbar.death()

		for tab, am, sm in zip(dataTabs, all_metals, std_metals):
			if tab.is_current():
				return am, sm


class DataNotebook(ttk.Notebook):
	"""
	Notebook for datasets
	"""
	def __init__(self, parent, num_tabs=6):
		super().__init__(parent)
		self.parent = parent # MethodsTab
		self.dtype = parent.dtype
		self.tabs = []

		for i in range(1,num_tabs+1):
			name = "Data {}".format(i)
			tab = DataTab(self, name)
			self.add(tab, text=name)
			self.tabs.append(tab)

		self.bind("<<NotebookTabChanged>>", lambda x: self.update(1))

	def update(self, strength=0):
		for tab in self.tabs:
			tab.update(strength)

	@property
	def current_tab(self):
		return self.tabs[self.index('current')]


class MethodsTab(tk.Frame):
	"""
	A single tab for PCS, RDC, PRE or CCR
	frm_fitopts : fitting options frame
	frm_data : data notebook containing data display and tensors
	"""
	def __init__(self, parent, dtype):
		super().__init__(parent)
		self.parent = parent # MethodsNotebook
		self.dtype = dtype

		self.frm_fopts = FittingOptionsFrame(self)
		self.ntb_data = DataNotebook(self)

		self.frm_fopts.grid(row=0,column=0)
		self.ntb_data.grid(row=1,column=0,columnspan=2)

		if self.dtype in ['PCS','PRE','CCR']:
			self.frm_fmul = MultipleFitFrame(self)
			self.frm_fmul.grid(row=0,column=1,sticky='NS')

	def update(self, strength=0):
		if self.is_current():
			self.ntb_data.update(strength)

	def is_current(self):
		if self==self.parent.current_tab:
			return True
		else:
			return False


class MethodsNotebook(ttk.Notebook):
	"""
	Contains tabs of PCS, RDC and CSA
	Generates numpy structured array templates for data storage and
	are stored in dictionary <templates>
	"""
	structdtype = {
		'PCS':np.dtype([
					('mdl', int   ),
					('use', bool  ),
					('atm', object),
					('cal', float ),
					('exp', float ),
					('err', float ),
					('idx', int   )]),
		'RDC':np.dtype([
					('mdl',  int   ),
					('use',  bool  ),
					('atm',  object),
					('atx',  object),
					('cal',  float ),
					('exp',  float ),
					('err',  float ),
					('idx',  int   )]),
		'PRE':np.dtype([
					('mdl', int   ),
					('use', bool  ),
					('atm', object),
					('cal', float ),
					('exp', float ),
					('err', float ),
					('idx', int   )]),
		'CCR':np.dtype([
					('mdl',  int   ),
					('use',  bool  ),
					('atm',  object),
					('atx',  object),
					('cal',  float ),
					('exp',  float ),
					('err',  float ),
					('idx',  int   )])}

	def __init__(self, parent):
		super().__init__(parent)
		self.parent = parent # MainGui
		self.dtypes = ['PCS','RDC','PRE','CCR']
		self.tabs = []
		self.templates = {i:None for i in self.dtypes}
		self.copied_tensor = None
		self.atomSet = set([])
		for dtype in self.dtypes:
			tab = MethodsTab(self, dtype)
			self.add(tab, text=dtype)
			self.tabs.append(tab)

		self.bind("<<NotebookTabChanged>>", lambda x: self.update(2))

	def update(self, strength=0):
		"""
		Templates are only updated if a change 
		is made to the PDB file or selection
		"""
		if strength>=2:
			self.make_templates()
		for tab in self.tabs:
			tab.update(strength)

	@property
	def current_tab(self):
		return self.tabs[self.index('current')]

	@property
	def frm_pdb(self):
		return self.parent.frm_coords.frm_pdb

	def make_templates(self):
		self.atomSet = set([])
		templates = {i:[] for i in self.dtypes}
		if not self.frm_pdb.prot:
			return
		nan = np.nan
		for m, a in self.frm_pdb.get_models_atoms():
			templates['PCS'].append(
				(m, False, a, nan, nan, nan, a.serial_number))
			templates['RDC'].append(
				(m, False, a, None, nan, nan, nan, a.serial_number))
			self.atomSet.add(a.name)

		for dtype in self.dtypes:
			template = np.array(templates[dtype], dtype=self.structdtype[dtype])
			self.templates[dtype] = template
		self.templates['PRE'] = self.templates['PCS'].copy()
		self.templates['CCR'] = self.templates['RDC'].copy()

		old_atomNames = set(self.frm_pdb.default_atom_selection)
		new_atomNames = self.atomSet - old_atomNames
		self.frm_pdb.atom_selection |= new_atomNames
		self.frm_pdb.default_atom_selection = list(self.atomSet)


class CSAFrame(tk.LabelFrame):
	"""
	Frame containing CSA tensor inputs
	Note this modifies class attribues for protein.CustomAtom dict 
	allowing all instances of CustomAtom to be affected globally
	"""
	def __init__(self, parent):
		super().__init__(parent, text='CSA tensor parameters')
		self.parent = parent # CoordinatesFrame
		Tooltip(self, tt['csa_frame'])

		for i, lbl in enumerate(['xx/ppm','yy/ppm','zz/ppm','\u03b2/\u00b0']):
			tk.Label(self, text=lbl).grid(row=0,column=i+1)

		for i, atom in enumerate(['H','N','C']):
			tk.Label(self, text=atom).grid(row=i+1,column=0)
			pas, beta = protein.CustomAtom.csa_lib[atom]
			for j in range(3):
				def wrap(a, b):
					a, b = atom, j
					def parse(value):
						protein.CustomAtom.csa_lib[a][0][b] = value*1E-6
					def disp():
						return protein.CustomAtom.csa_lib[a][0][b]*1E6
					return parse, disp

				parse, disp = wrap(atom, j)
				NumericEntry(self, parse, disp).grid(row=i+1,column=j+1)

			def wrap(a):
				a = atom
				def parse(value):
					pas, beta = protein.CustomAtom.csa_lib[a] 
					protein.CustomAtom.csa_lib[a] = pas, value*(np.pi/180.)
				def disp():
					return protein.CustomAtom.csa_lib[a][1]*(180./np.pi)
				return parse, disp
			parse, disp = wrap(atom)
			NumericEntry(self, parse, disp).grid(row=i+1,column=4)


class CSAPopup(Popup):
	"""A popup allowing CSA parameters to be set"""
	def __init__(self, parent):
		title = "Set CSA parameters"
		super().__init__(parent, title)
		CSAFrame(self).pack()
		ttk.Button(self, text='Ok', command=self.death).pack()
		self.update()

	def save(self):
		self.death()


class PDBFrame(tk.LabelFrame):
	"""
	Frame for PDB input and model selection
	prot : the protein object
	models : the parsed models
	"""
	def __init__(self, parent):
		super().__init__(parent, text='PDB & Models')
		self.parent = parent # CoordinatesFrame
		self.num_chars = 30
		self.prot = None
		self.models = set([])

		self.default_atom_selection = ['H','N','C','O']
		self.default_resi_selection = list(protein.standard_aa_names)

		self.atom_selection = set(self.default_atom_selection)
		self.resi_selection = set(self.default_resi_selection)
		self.reset_selection()

		b = ttk.Button(self, text='Read PDB file',command=self.read_pdb)
		Tooltip(b, tt['read_pdb'])
		b.grid(row=0,column=0,sticky='EW')

		self.lbl_pdb_file = tk.Label(self,text=" "*self.num_chars)
		self.lbl_pdb_file.grid(row=0,column=1)

		# tk.Label(self, text="Models for pdb:").grid(row=1, column=0)
		b = ttk.Button(self, text='Set Models', 
			command=self.parse_models)
		Tooltip(b, tt['parse_models'])
		b.grid(row=1,column=0,sticky='EW')

		self.txt_models = CustomTextDefaultEntry(self,
			"Choose models. Default: all", returnKey=self.parse_models)
		Tooltip(self.txt_models, tt['select_models'])
		self.txt_models.grid(row=1, column=1)

		ttk.Separator(self, orient='vertical').grid(
			row=0,column=3,rowspan=3,sticky='NS',padx=3)

		self.selction = None
		b = ttk.Button(self, text='Atom Selection',command=self.change_selection)
		Tooltip(b, tt['selection'])
		b.grid(row=0,column=4,sticky='EW')

		self.frm_csa = CSAFrame(self)
		b = ttk.Button(self, text='Set CSA',command=self.set_csa)
		# Tooltip(b, tt['selection'])
		b.grid(row=1,column=4,sticky='EW')

	def has_models(self):
		if len(self.models)>1:
			return True
		else:
			return False

	def update_views(self, strength=2):
		"""
		When the PDB is changed, all experimental data is parsed
		and all views are updated
		"""
		self.parent.parent.update(strength)

	def reset_selection(self):
		self.atom_selection = set(self.default_atom_selection)
		self.resi_selection = set(self.default_resi_selection)

	def read_pdb(self):
		"""Open fild dialog to fetch the PDB file path, load and parse"""
		fileName = filedialog.askopenfilename(
			title="Choose PDB file",
			defaultextension='.pdb',
			filetypes=[('PDB file','.pdb'),('All files','.*')])

		if fileName:
			self.prot = protein.load_pdb(fileName)
			self.lbl_pdb_file.config(
				text=format_path(fileName, self.num_chars))
			self.parse_models()

	def parse_models(self):
		"""
		Printer pages type notation to specify the desired models
		If incorrect notation is entered, the models default to all
		"""
		if self.prot:
			prot_models = set([i.id for i in self.prot])
		else:
			self.models = set([])
			self.store()
			return

		s = self.txt_models.get()
		if s==self.txt_models.default_text:
			self.models = prot_models
			self.store()
			return

		try:
			models = unpack_ranges(s)
		except ValueError:
			messagebox.showerror("Error", 
				"Error parsing models: Default selection set to all")
			self.models = prot_models
			self.txt_models.clear()
			self.store()
			return None

		diff = models - prot_models
		if len(diff)>0:
			messagebox.showwarning("Warning", 
				"Additional models not found in PDB:\n{}".format(list(diff)))
		self.models = models & prot_models
		if len(self.models)==0:
			messagebox.showwarning("Warning", 
				"No models selected: Default selection set to all")
			self.models = prot_models
			self.txt_models.clear()
		self.store()

	def store(self):
		print("Parsed models:\n{}".format(self.models))
		self.update_views(2)

	def get_atoms(self):
		"""Generator for all atoms in the protein"""
		if self.prot:
			for m in self.models:
				for a in self.prot[m].get_atoms():
					yield a
		else:
			yield

	def get_models_atoms(self):
		if self.prot:
			for m in self.models:
				for a in self.prot[m].get_atoms():
					yield m, a
		else:
			yield

	def change_selection(self, *args):
		SelectionPopup(self)

	def set_csa(self, *args):
		CSAPopup(self)


class CoordinatesFrame(tk.Frame):
	"""The frame containing PDB input and CSA info"""
	def __init__(self, parent):
		super().__init__(parent)
		self.parent = parent # MainGui
		self.frm_pdb = PDBFrame(self)
		self.frm_pdb.grid(row=0,column=0,sticky="NESW")


class MainGui(tk.Frame):
	"""
	A Frame holding the entire application
	frm_coords : contains PDB and CSA frames
	ntb_methods : contains tabs of PCS, RDC and CSA
	"""
	def __init__(self, root):
		super().__init__(root)
		self.frm_coords = CoordinatesFrame(self)
		self.frm_coords.grid(row=0,column=0)
		self.ntb_methods = MethodsNotebook(self)
		self.ntb_methods.grid(row=1,column=0,sticky='EW')


	def update(self, strength=0):
		"""
		0 : update all fields, checkboxes and labels
		1 : update all data views
		2 : parse experimental data
		"""
		self.ntb_methods.update(strength)


tooltips_raw = """
read_pdb : Load a PDB file containing protein coordinates for all atoms
select_models : Select models to be used during fitting. This can be specified using printer-style formatting e.g. '1,3,7' or '1-5' or '1-5,7,9'. Be sure to click the "Parse models" button before continuing.
parse_models : Selects only the desired models from the PDB file
selection : Select the atoms and residues to be used during fitting
csa_frame : Set chemical shift anisotropy tensor parameters for backbone H, N and C nuclei. These are set via the three principle axes and single angle. Default values are taken from G. Cornilescu and A. Bax "Measurement of proton, nitrogen, and carbonyl chemical shielding anisotropies in a protein dissolved in a dilute liquid crystalline phase" J. Am. Chem. Soc. 2000, 122, 10143-10154.
svd : Single value decomposition is used to analytically solve the tensor anisotropy and orientation over a spherical grid of positions. The best position and tensor are then taken. Note that RACS and RADS cannot be included in this calculation.
nlr : Non-linear regression is used to solve a least-squares algorithm by gradient descent. This is generally used for refinement after an initial guess by the SVD grid search. Note that RACS and RADS can be taken into account by this algorithm.
ref : Allows fitting of an offset that shifts the entire PCS list by a given value. This may arise due to referencing errors between diamagnetic and paramagnetic peak lists when calculating the PCS. This should always be used when there are many PCS values available.
pos : Option to disable fitting of the tensor position. Note that the SVD grid search for PCS will collapse to a single point if the position is constrained.
rad : Radius of the SVD grid search sphere. This is taken with origin about the initial tensor position
den : The points per Angstrom to be taken in the SVD grid search sphere.
eav : When selected, calculations and fitting are performed with ensemble averaging between models of the PDB file.
racs : Residual anisotropic chemical shifts are included in the calculation. Fitting with RACS is only achieved with the NLR algorithm.
rads : Residual anisotropic dipolar shifts are included in the calculation. Fitting with RADS is only achieved with the NLR algorithm.
taur : The rotational correlation time is included as a parameter for fitting during the calculation.
taue : The electronic relaxation time is included as a parameter for fitting during the calculation.
dsa : Dipolar shielding anisotropy (DSA) commonly called the 'Curie Spin' relaxation theory is included for the calculation. Note that this method assumes a fast electronic relaxation time (T1e << TauR). Therefore calculations are independent of the T1e parameter. All calculations for DSA theory will include anisotropy as governed by the DeltaChi tensor of the paramagnetic centre.
sbm : Solomon-Bloembergen-Morgan (SBM) relaxation theory is included for the calculation. SBM theory uses an effective correlation time that is takes into account the rotational correlation time (TauR) and the electronic relaxation time (T1e), so both must be set when using this method.
csa : Chemical Shift Anisotropy (CSA) cross-correlated relaxation effects are included for the caluclation. For auto-relaxation of nuclear magnetisation, the CSA may only induce cross-correlated relaxation with the DSA mechanism. Therefore, the DSA relaxation method must be selected to observe CSA cross-correlation effects. The final output value is a sum of the DSA plus the DSAxCSA values, however the relaxation induced by only the CSA is not included. This means only the pure paramagnetic effect is reported.
dchi : The DeltaChi tensor is included for fitting, allowing the axial and rhombic parameters to be fit during the PRE calculation. Note that the DeltaChi tensor has only a very small effect on the PRE and therefore may lead to unwieldy results if included in the fitting.
iso : The isotropic Chi tensor is included for fitting during the caluclation. This defines the magnitude of the paramagnetic dipole.
"""

tt = {}
for line in tooltips_raw.split('\n'):
	if line:
		key, value = line.split(':',1)
		tt[key.strip()] = value.strip()

# Run the application
def run():
	settings = platform_settings()
	root = tk.Tk()
	root.title('Paramagpy GUI - v {}'.format(paramagpy.__version__))

	if hasattr(sys, '_MEIPASS'):
		dataDir = sys._MEIPASS
	else:
		dataDir = os.path.dirname(__file__)

	try:
		icon_path = os.path.join(dataDir, "icon.gif")
		icon = tk.PhotoImage(file=icon_path)
		root.tk.call('wm', 'iconphoto', root._w, icon)
	except tk.TclError:
		print("Note: icon not found")

	main = MainGui(root)
	main.pack(expand=True, fill='both')
	while True:
		try:
			root.mainloop()
			break
		except UnicodeDecodeError:
			pass

if __name__ == "__main__":
	run()
 
