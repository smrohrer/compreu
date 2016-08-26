from Tkinter import *
import ASEgraphene

#Requires 'hovering cursor' tooltip for all major components

root = Tk()
root.title("Graphene Nitrogenator")
Label(root, text="GUI for Nitrogenated Graphene Energy Calculations").pack()


class sheet_param:
	def __init__(self, master):
		global horizon_sheet_variable
		horizon_sheet_variable = StringVar()
		global vertical_sheet_variable
		vertical_sheet_variable = StringVar()
		param_frame = Frame(master)
		Label(param_frame, text="horizontal").grid(row=1, column=1)
		Label(param_frame, text="vertical").grid(row=1, column=3)
		Label(param_frame, text="Sheet Dimensions").grid(row=2, sticky=E)
		Label(param_frame, text="x").grid(row=2, column=2)
		Entry(param_frame, width=2, textvariable=horizon_sheet_variable).grid(row=2, column=1)
		Entry(param_frame, width=2, textvariable=vertical_sheet_variable).grid(row=2, column=3)
		param_frame.pack()


class checkbox:
	def __init__(self,master):
		self.cbox_frame = Frame(master)
		self.cbox_frame.pack()
		
	def symtrc_cbox(self):
		global symmetry_var
		symmetry_var = IntVar()
		Checkbutton(self.cbox_frame, text="Make molecule symmetric", variable=symmetry_var).pack(anchor=W)

	def unsaturate_cbox(self):
		global unsat_var
		unsat_var = IntVar()
		Checkbutton(self.cbox_frame, text="Leave molecule unsaturated", variable=unsat_var, onvalue=1, offvalue=0).pack(anchor=W)
		#unsat_var_string = str(unsat_var.get())

	def opt_geo_cbox(self):
		opt_geo_var = IntVar()
		Checkbutton(self.cbox_frame, text="Optimize Geometry", variable=opt_geo_var).pack(anchor=W)
		global selected_geometry_opt
		selected_geometry_opt = int(opt_geo_var.get())
#can't seem to get checkboxes to align


class drop_list:
	def __init__(self, master):
		self.drop_list_frame = Frame(master)
		self.drop_list_frame.pack()
		self.method_var = StringVar()
		self.calculator_var = StringVar()

	def calc_combobox(self):
		self.calculator_var.set("ORCA")
		OptionMenu(self.drop_list_frame, self.calculator_var, "ORCA", "NWChem", "Gaussian").grid(row=0, column=1)
		Label(self.drop_list_frame, text="Calculator").grid(row=0, sticky=E)
		global selected_calculator
		selected_calculator = str(self.calculator_var.get())

	def calc_method(self):
		self.method_var.set("am1")
		OptionMenu(self.drop_list_frame, self.method_var, "am1", "DFT").grid(row=1, column=1)
		Label(self.drop_list_frame, text="method").grid(row=1, sticky=E)
		global selected_calc_method
		selected_calc_method = str(self.method_var.get())



class button:
	def __init__(self, master):
		self.button_frame = Frame(master)
		self.button_frame.pack(expand="yes")

	def bottom_buttons(self):
		Button(self.button_frame, text="Calculate", command=gui_calculate).grid(row=0, column=3)
		Button(self.button_frame, text="View in Avogadro", command=gui_view).grid(row=0, column=0)

	def test_button(self):
		Button(self.button_frame, text="Print Variable", command=self.print_var).grid(row=1)

	def print_var(self):
		string = str(horizon_sheet_variable.get())
		print string


def build_param_frame(master):
	lbl_bld_frame = LabelFrame(master, text="Graphene Builder Parameters", padx=5, pady =5)
	lbl_bld_frame.pack(expand="yes")
	sheet_param(lbl_bld_frame)
	checkbox(lbl_bld_frame).symtrc_cbox()
	checkbox(lbl_bld_frame).unsaturate_cbox()

def calc_param_frame(master):
	lbl_calc_frame = LabelFrame(master, text="Calculation Parameters", padx=5, pady=5)
	lbl_calc_frame.pack(expand="yes")
	drop_list(lbl_calc_frame).calc_combobox()
	drop_list(lbl_calc_frame).calc_method()
	checkbox(lbl_calc_frame).opt_geo_cbox()

def gui_view():
        ##setting sheet dimension parameter variables from ASEgraphene.py to local variables
        horizontal_dimension = int(horizon_sheet_variable.get())
        vertical_dimension = int(vertical_sheet_variable.get())
        symmetry_int = int(symmetry_var.get())
        atoms = ASEgraphene.build_sheet(horizontal_dimension, vertical_dimension, symmetry=symmetry_int)
        unsat_int = int(unsat_var.get())
        if unsat_int==0:
        	ASEgraphene.daves_super_saturate(atoms)
        elif unsat_int==1:
        	pass
        ASEgraphene.view(atoms, viewer="avogadro")

def gui_calculate():
	if selected_calculator=="ORCA":

		horizontal_dimension = int(horizon_sheet_variable.get())
		vertical_dimension = int(vertical_sheet_variable.get())
		symmetry_int = int(symmetry_var.get())
		atoms = ASEgraphene.build_sheet(horizontal_dimension, vertical_dimension, symmetry=symmetry_int)
		unsat_int = int(unsat_var.get())

		if unsat_int==0:
			ASEgraphene.daves_super_saturate(atoms)
		elif unsat_int==1:
			pass

		if selected_geometry_opt==1:
			opt_geom_truefalse = True
		elif selected_geometry_opt==0:
			opt_geom_truefalse = False

		ASEgraphene.calc_edge_nitrogens(horizontal_dimension, vertical_dimension, method=selected_calc_method, optimize_geometry=opt_geom_truefalse, make_symmetric=selected_geometry_opt)
	else:
		pass

	



build_param_frame(root)
calc_param_frame(root)
button(root).bottom_buttons()

button(root).test_button()
root.mainloop()