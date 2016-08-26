import pickle
import itertools

nx = 5
nzs = [3,5]

binaries = []
reactions = []

# find possible combinations, expressed in binary
# 0 = no H, 1 = H added
for i in range(0, 4 ** nx):
    binary = bin(i)[2:]
    while len(binary) < 2*nx:
        binary = '0' + binary
    mirror = binary[:nx][::-1]+"_"+binary[-nx:][::-1]
    if (not mirror in binaries) and (binary.count('1')%2 == 0):
        binaries.append(binary[:nx]+"_"+binary[-nx:])

# for each combination, make a list of all possible combinations formed by adding 2 H's
for binary in binaries:
    # find all spaces without H (i.e. 0 in the binary string)
    spaces = []
    for i,n in enumerate(binary):
        if n == '0':
            spaces.append(i)
    # find all combinations of 2 spaces to add H
    spaces_to_add = list(itertools.combinations(spaces, 2))
    products = []
    for combo in spaces_to_add:
        product = binary[:]
        for space in combo:
            product = product[:space]+'1'+product[space+1:]
        if product in binaries:
            products.append(product)
    for product in products:
        reactions.append({
            "reactant": binary,
            "product": product
        })

for nz in nzs:
    dataset = pickle.load(open("NH_combinations"+str(nx)+"x"+str(nz)))
    data = {}
    for sheet in dataset:
        data[sheet["name"]] = {
            "scf": sheet["scf"],
            "homo": sheet["homo"],
            "lumo": sheet["lumo"],
            "nH": list(sheet["Zs"]).count(1) - 2*(nz+1)
        }

    for reaction in reactions:
        try:
            reaction["nz"+str(nz)] = {
                "initial_nH": data[reaction["reactant"]]["nH"],
                "delta_E": data[reaction["product"]]["scf"] - data[reaction["reactant"]]["scf"],
                "initial_LUMO": data[reaction["reactant"]]["lumo"],
                "final_HOMO": data[reaction["product"]]["homo"]
            }
        except KeyError:
            print "KeyError "+reaction["reactant"]+" "+reaction["product"]

pickle.dump(reactions, open("rxndata", "w"))