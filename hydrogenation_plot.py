import matplotlib.pyplot as plt
import pickle

H2_TOTAL_ENERGY = -27.536144
INITIAL_SHEET_TOTAL_ENERGY = -6519.32774

pickles = ["hydrogenate_random_final", "NH_combinations"]

N_indicies = [4,11,12,19,20,27,28,35,36,43]

nH = []
scf = []
ref_to_0 = []

for p in pickles:
    dataset = pickle.load(open(p))
    for data in dataset:
        scf.append(data["scf"])
        n_NH = 0
        for index in N_indicies:
            if data["H_attachment"][index] == 1:
                n_NH += 1
        nH.append(n_NH)
        sheet_nH2 = scf[-1] + (5-data["H_attachment"].count(1)/2.0)*H2_TOTAL_ENERGY
        ref_to_0.append(sheet_nH2 - INITIAL_SHEET_TOTAL_ENERGY - 5.0 * H2_TOTAL_ENERGY)

plt.figure()
plt.plot(nH, scf, 'o')
plt.savefig("scf.jpg")
plt.figure()
plt.plot(nH, ref_to_0, 'o')
plt.savefig("scf_ref_to_0.jpg")
