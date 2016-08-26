from sklearn import linear_model
import numpy as np
import pickle

H2_TOTAL_ENERGY = -27.536144
INITIAL_SHEET_TOTAL_ENERGY = -6519.32774

sheet_map = pickle.load(open("5x5map"))
N_indices = [6,18,30,42,54,17,29,41,53,65]

nH = []
H_attachment = []
H_gaussian = []
scf = []
homo = []
lumo = []
ref_to_0 = []

dataset = pickle.load(open("NH_combinations5x3"))
for data in dataset:
    scf.append(data["scf"])
    lumo.append(data["lumo"])
    homo.append(data["homo"])
    H_attachment.append([float(val) for val in data["name"].replace("_","")])
    H_gaussian.append([])

    for i1,site1 in enumerate(H_attachment[len(H_attachment)-1]):
        ri = sheet_map[N_indices[i1]]
        sum = 0
        for i2,site2 in enumerate(H_attachment[len(H_attachment)-1]):
            if site2 == 1:
                rH = sheet_map[N_indices[i2]]
                sum += np.exp(-(np.linalg.norm(ri-rH)**2))
        H_gaussian[len(H_gaussian)-1].append(sum)
        del sum

    sheet_nH2 = scf[-1] + (5-data["name"].count('1')/2.0)*H2_TOTAL_ENERGY
    ref_to_0.append(sheet_nH2 - INITIAL_SHEET_TOTAL_ENERGY - 5.0 * H2_TOTAL_ENERGY)

regr = linear_model.LinearRegression()
regr.fit(H_gaussian[:-20], scf[:-20])
prediction = regr.predict(H_gaussian[-20:])
print "--------------------------------------"
print "TOTAL ENERGY"
print "Prediction\tActual"
for i in range(0,20):
    print str(prediction[i])+"\t"+str(scf[len(scf)-20+i])
print "MAE="+str((1.0/len(prediction))*sum(np.absolute(np.subtract(prediction,scf[-20:]))))
print "RMSE="+str(np.sqrt((1.0/len(prediction))*sum(np.subtract(prediction,scf[-20:])**2)))

regr = linear_model.LinearRegression()
regr.fit(H_gaussian[:-20], homo[:-20])
prediction = regr.predict(H_gaussian[-20:])
print "--------------------------------------"
print "HOMO"
print "Prediction\tActual"
for i in range(0,20):
    print str(prediction[i])+"\t"+str(homo[len(homo)-20+i])
print "MAE="+str((1.0/len(prediction))*sum(np.absolute(np.subtract(prediction,homo[-20:]))))
print "RMSE="+str(np.sqrt((1.0/len(prediction))*sum(np.subtract(prediction,homo[-20:])**2)))

regr = linear_model.LinearRegression()
regr.fit(H_gaussian[:-20], lumo[:-20])
prediction = regr.predict(H_gaussian[-20:])
print "--------------------------------------"
print "LUMO"
print "Prediction\tActual"
for i in range(0,20):
    print str(prediction[i])+"\t"+str(lumo[len(lumo)-20+i])
print "MAE="+str((1.0/len(prediction))*sum(np.absolute(np.subtract(prediction,lumo[-20:]))))
print "RMSE="+str(np.sqrt((1.0/len(prediction))*sum(np.subtract(prediction,lumo[-20:])**2)))

regr = linear_model.LinearRegression()
regr.fit(H_gaussian[:-20], ref_to_0[:-20])
prediction = regr.predict(H_gaussian[-20:])
print "--------------------------------------"
print "REF TO 0"
print "Prediction\tActual"
for i in range(0,20):
    print str(prediction[i])+"\t"+str(ref_to_0[len(ref_to_0)-20+i])
print "MAE="+str((1.0/len(prediction))*sum(np.absolute(np.subtract(prediction,ref_to_0[-20:]))))
print "RMSE="+str(np.sqrt((1.0/len(prediction))*sum(np.subtract(prediction,ref_to_0[-20:])**2)))
