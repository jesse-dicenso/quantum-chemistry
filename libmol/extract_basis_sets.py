# Format basis sets from basis set exchange (Q-Chem format)

basis_set = "STO-3G"

rawfile = basis_set + "_raw.txt"
lines = []
cutoffs = [0]
with open(rawfile, "r") as f:
    lines = f.readlines()
    k=1
    for l in range(0, len(lines)):
        if("****" in lines[l]): cutoffs.append(k)
        k+=1
        if (lines[l][0]!='D'):
           lines[l] = lines[l].replace('D', 'e')

for i in range(0, len(cutoffs)-1):
    element = lines[cutoffs[i]].split()[0]
    lcount = 0;
    for j in range(cutoffs[i]+1,cutoffs[i+1]-1):
        temp = lines[j].split()[0]
        if((temp=="S") or (temp=="SP") or (temp=="P") or (temp=="D") or (temp=="F")):
            lcount+=1

    with open(f"{basis_set}/{element}", "w") as g:
        g.write(f"{str(lcount)}\n")
        for j in range(cutoffs[i]+1,cutoffs[i+1]-1):
            g.write(lines[j])
