from Tkinter import *

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
		Checkbutton(self.cbox_frame, text="Make molecule symmetric").pack()
		

sh_pm = sheet_param(root)
symtrc_check = checkbox(root)
symtrc_check.symtrc_cbox()
root.mainloop()
