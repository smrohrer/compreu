from Tkinter import *


#button for instant open molecule in avogadro, 'Go' button shoud be called Calculate
root = Tk()
root.title("Graphene Nitrogenator")
Label(root, text="GUI for Nitrogenated Graphene Energy Calculations").pack()


class sheet_param:
	def __init__(self, master):
	
		param_frame = Frame(master)
		Label(param_frame, text="horizontal").grid(row=1, column=1)
		Label(param_frame, text="vertical").grid(row=1, column=3)
		Label(param_frame, text="Sheet Dimensions").grid(row=2, sticky=E)
		Label(param_frame, text="x").grid(row=2, column=2)
		Entry(param_frame, width=2).grid(row=2, column=1)
		Entry(param_frame, width=2).grid(row=2, column=3)
		param_frame.pack()
		

class checkbox:
	def __init__(self,master):
		self.cbox_frame = Frame(master)
		self.cbox_frame.pack()

	def symtrc_cbox(self):
		Checkbutton(self.cbox_frame, text="Make molecule symmetric").pack(anchor=W)

	def unsaturate_cbox(self):
		Checkbutton(self.cbox_frame, text="Leave molecule unsaturated").pack(anchor=W)

	def opt_geo_cbox(self):
		Checkbutton(self.cbox_frame, text="Optimize Geometry").pack(anchor=W)
#can't seem to get checkboxes to align


class drop_list:
	def __init__(self, master):
		self.drop_list_frame = Frame(master)
		self.drop_list_frame.pack()

	def calc_combobox(self):
		calculator_var = StringVar()
		calculator_var.set("ORCA")
		OptionMenu(self.drop_list_frame, calculator_var, "ORCA", "NWChem", "Gaussian").grid(row=0, column=1)
		Label(self.drop_list_frame, text="Calculator").grid(row=0, sticky=E)

	def calc_method(self):
		method_var = StringVar()
		method_var.set("am1")
		OptionMenu(self.drop_list_frame, method_var, "am1", "DFT").grid(row=1, column=1)
		Label(self.drop_list_frame, text="method").grid(row=1, sticky=E)


class button:
	def __init__(self, master):
		self.button_frame = Frame(master)
		self.button_frame.pack(expand="yes")

	def bottom_buttons(self):
		Button(self.button_frame, text="Calculate").grid(row=0, column=3)
		Button(self.button_frame, text="View in Avogadro").grid(row=0, column=0, padx="50")


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


<<<<<<< HEAD
build_param_frame(root)
calc_param_frame(root)
button(root).bottom_buttons()
root.mainloop()
=======
sh_pm = sheet_param(root)
symtrc_check = checkbox(root)
symtrc_check.symtrc_cbox()
root.mainloop()
>>>>>>> da48a0e0c761313a778069496e76d8e4e37421e9
