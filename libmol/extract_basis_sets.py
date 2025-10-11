# Format basis sets from basis set exchange (Q-Chem format)
basis_set = "def2-SVP"

rawfile = basis_set + "_raw.txt"
lines = []
cutoffs = [0]
with open(rawfile, "r") as f:
    lines = f.readlines()
    k=1
    for l in lines:
        if("****" in l): cutoffs.append(k)
        k+=1

for i in range(0, len(cutoffs)-1):
    element = lines[cutoffs[i]].split()[0]
    lcount = 0;
    for j in range(cutoffs[i]+1,cutoffs[i+1]-1):
        temp = lines[j].split()[0]
        if((temp=="S") or (temp=="SP") or (temp=="P") or (temp=="D") or (temp=="F")):
            lcount+=1

    with open(f"def2-SVP/{element}", "w") as g:
        g.write(f"{str(lcount)}\n")
        for j in range(cutoffs[i]+1,cutoffs[i+1]-1):
            g.write(lines[j])
