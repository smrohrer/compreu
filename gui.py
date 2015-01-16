from Tkinter import *
import ASE-Graphene-Script

#Requires 'hovering cursor' tooltip for all major components

root = Tk()
root.title("Graphene Nitrogenator")
Label(root, text="GUI for Nitrogenated Graphene Energy Calculations").pack()


class sheet_param:
	def __init__(self, master):
		horizon_sheet_variable = StringVar()
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
		self.symmetry__var = IntVar()
		Checkbutton(self.cbox_frame, text="Make molecule symmetric", variable=self.symmetry__var).pack(anchor=W)

	def unsaturate_cbox(self):
		self.unsat_var = IntVar()
		Checkbutton(self.cbox_frame, text="Leave molecule unsaturated", variable=self.unsat_var, onvalue=1, offvalue=0).pack(anchor=W)

	def opt_geo_cbox(self):
		Checkbutton(self.cbox_frame, text="Optimize Geometry").pack(anchor=W)
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
		selected_calculator = self.calculator_var.get()

	def calc_method(self):
		self.method_var.set("am1")
		OptionMenu(self.drop_list_frame, self.method_var, "am1", "DFT").grid(row=1, column=1)
		Label(self.drop_list_frame, text="method").grid(row=1, sticky=E)



class button:
	def __init__(self, master):
		self.button_frame = Frame(master)
		self.button_frame.pack(expand="yes")

	def bottom_buttons(self):
		Button(self.button_frame, text="Calculate").grid(row=0, column=3)
		Button(self.button_frame, text="View in Avogadro").grid(row=0, column=0)


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
	


build_param_frame(root)
calc_param_frame(root)
button(root).bottom_buttons()
root.mainloop()