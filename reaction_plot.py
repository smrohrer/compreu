import pickle
import matplotlib.pyplot as plt
import numpy as np

nx = 5
nzs = [3,5]
colors = ['b','r','g','y']

data = pickle.load(open("rxndata"))

plt.figure()
for inz,nz in enumerate(nzs):
    delta_E = []
    homo = []
    for rxn in data:
        try:
            delta_E.append(rxn["nz"+str(nz)]["delta_E"])
            homo.append(rxn["nz"+str(nz)]["final_HOMO"])
        except KeyError:
            print "KeyError "+rxn["reactant"]+" "+rxn["product"]
    plt.plot(delta_E, homo, colors[inz]+".", label="nz="+str(nz))
plt.legend(loc="best")
plt.xlabel("Delta E rxn")
plt.ylabel("HOMO energy of product")
plt.savefig("delta_E_vs_final_HOMO.jpg")

plt.figure()
for inz,nz in enumerate(nzs):
    delta_E = []
    lumo = []
    for rxn in data:
        try:
            delta_E.append(rxn["nz"+str(nz)]["delta_E"])
            lumo.append(rxn["nz"+str(nz)]["initial_LUMO"])
        except KeyError:
            print "KeyError "+rxn["reactant"]+" "+rxn["product"]
    plt.plot(delta_E, lumo, colors[inz]+".", label="nz="+str(nz))
plt.legend(loc="best")
plt.xlabel("Delta E rxn")
plt.ylabel("LUMO energy of starting material")
plt.savefig("delta_E_vs_initial_LUMO.jpg")

plt.figure()
for inz,nz in enumerate(nzs):
    homo = []
    lumo = []
    for rxn in data:
        try:
            homo.append(rxn["nz"+str(nz)]["final_HOMO"])
            lumo.append(rxn["nz"+str(nz)]["initial_LUMO"])
        except KeyError:
            print "KeyError "+rxn["reactant"]+" "+rxn["product"]
    plt.plot(homo, lumo, colors[inz]+".", label="nz="+str(nz))
plt.legend(loc="best")
plt.xlabel("HOMO energy of product")
plt.ylabel("LUMO energy of starting material")
plt.savefig("final_HOMO_vs_initial_LUMO.jpg")

plt.figure()
for inz,nz in enumerate(nzs):
    nH = []
    lumo = []
    for rxn in data:
        try:
            nH.append(rxn["nz"+str(nz)]["initial_nH"])
            lumo.append(rxn["nz"+str(nz)]["initial_LUMO"])
        except KeyError:
            print "KeyError "+rxn["reactant"]+" "+rxn["product"]
    plt.plot(nH, lumo, colors[inz]+".", label="nz="+str(nz))
plt.legend(loc="best")
plt.xlabel("nH")
plt.ylabel("LUMO energy")
plt.savefig("nH_vs_LUMO.jpg")

plt.figure()
for inz,nz in enumerate(nzs):
    nH = []
    homo = []
    for rxn in data:
        try:
            nH.append(rxn["nz"+str(nz)]["initial_nH"]+2)
            homo.append(rxn["nz"+str(nz)]["final_HOMO"])
        except KeyError:
            print "KeyError "+rxn["reactant"]+" "+rxn["product"]
    plt.plot(nH, homo, colors[inz]+".", label="nz="+str(nz))
plt.legend(loc="best")
plt.xlabel("nH")
plt.ylabel("HOMO energy")
plt.savefig("nH_vs_HOMO.jpg")

sorted_rxns = [[] for i in range(nx)]
for rxn in data:
    sorted_rxns[rxn["nz"+str(nzs[0])]["initial_nH"]/2].append(rxn)
for nz in nzs:
    print "nz="+str(nz)
    reactant = "00000_00000"
    for set in sorted_rxns:
        minimum = 0.0
        product = ""
        for rxn in set:
            if rxn["reactant"] == reactant:
                if rxn["nz"+str(nz)]["delta_E"] < minimum:
                    minimum = rxn["nz"+str(nz)]["delta_E"]
                    product = rxn["product"]
        print reactant+" -> "+product+"\tdelta E = "+str(minimum)
        reactant = product