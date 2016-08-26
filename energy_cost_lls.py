import numpy as np
import pickle
data = pickle.load(open("NH_combinations5x3"))

y = []
X = []

# the known value, stored in y, is the total energy for each sheet.
# the unknowns are the constant term, c, and the "energy cost" per:
# p0, p1, p2 = number of hydrogens at the 0,1,2 positions when each set of edge nitrogens is mapped as 2 1 0 1 2
# s1, s2, s3, s4 = number of pairs of hydrogens placed 1, 2, 3, and 4 spaces apart on the same side
# t = number of groups of 3 next to each other
# a = number of pairs across from each other

def get_coefficients(pattern):
    p0 = [pattern[0][2], pattern[1][2]].count('1')
    p1 = [pattern[0][1], pattern[0][3], pattern[1][1], pattern[1][3]].count('1')
    p2 = [pattern[0][0], pattern[0][4], pattern[1][0], pattern[1][4]].count('1')

    s1 = 0
    s2 = 0
    s3 = 0
    s4 = 0
    t = 0

    for side in pattern:

        for i in range(0, len(side) - 1):
            if side[i] == '1' and side[i + 1] == '1':
                s1 += 1

        for i in range(0, len(side) - 2):
            if side[i] == '1' and side[i + 2] == '1':
                s2 += 1
                if side[i + 1] == '1':
                    t += 1

        for i in range(0, len(side) - 3):
            if side[i] == '1' and side[i + 3] == '1':
                s3 += 1

        for i in range(0, len(side) - 4):
            if side[i] == '1' and side[i + 4] == '1':
                s4 += 1

    a = 0
    for i in range(0, len(pattern[0])):
        if pattern[0][i] == '1' and pattern[1][i] == '1':
            a += 1

    return [1, p0, p1, p2, s1, s2, s3, s4, t, a]

for sheet in data[:-20]:

    y.append(sheet["scf"])

    pattern = sheet["name"].split('_')
    X.append(get_coefficients(pattern))

sol = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(X),X)),np.transpose(X)),y)

print "c = "+str(sol[0])
print "p0 = "+str(sol[1])
print "p1 = "+str(sol[2])
print "p2 = "+str(sol[3])
print "s1 = "+str(sol[4])
print "s2 = "+str(sol[5])
print "s3 = "+str(sol[6])
print "s4 = "+str(sol[7])
print "t = "+str(sol[8])
print "a = "+str(sol[9])

# print "--------------------------------------"
# print "Total Energy"
# print "Prediction\tActual"
#
# predictions = []
# actual = []
# for sheet in data[-20:]:
#     actual.append(sheet["scf"])
#     pattern = sheet["name"].split('_')
#     x = get_coefficients(pattern)
#     predictions.append(np.dot(sol,x))
#     print sheet["name"]+"\t"+str(x)
#
# for i in range(0,20):
#     print str(predictions[i])+"\t"+str(actual[i])
#
# print "MAE="+str((1.0/len(predictions))*sum(np.absolute(np.subtract(predictions,actual))))
# print "RMSE="+str(np.sqrt((1.0/len(predictions))*sum(np.subtract(predictions,actual)**2)))
