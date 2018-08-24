import urllib.request
import shutil
import sys
import subprocess
import re
import csv

#Stage I

def fetcher(url, pdb): #only tested on up to 108 sequences, theoretically works on any number. 
    try:
        urllib.request.urlretrieve(url, url.split("/").pop())#urllib3?
    except IOError as e:
        print("IO Error")
    except Exception as e:
        print("Exception Error")

def commander(pdb, chain):
    #subprocess.run(["pdbsplitchains " + pdb + ".pdb"], shell = True) #commented out as I was kindly given copies of all desired chains
    subprocess.run(["cat Fourecks.txt >> pro.txt"], shell = True)
    subprocess.run(["pdbtorsions " + pdb + chain + ".pdb | grep -C 4 PRO | cat >> pro.txt"], shell = True)
    subprocess.run(["pdbtorsions " + pdb + chain + ".pdb | wc -l | cat >> count_total.txt"], shell = True) #count-related
    subprocess.run(["echo -2 >> count_total.txt"], shell = True) #count-related
    subprocess.run(["pdbtorsions " + pdb + chain + ".pdb | awk '$5 >= -10 && $5 <= 10'| wc -l | cat >> count_tcis.txt"], shell = True) #count-related
    subprocess.run(["rm " + pdb + "*"], shell = True) #the part that removes processed chains for speed reasons
    subprocess.run(["cat Fourecks.txt >> pro.txt"], shell = True)
    
def protidy():
    subprocess.run(["sed -i '/--/d' pro.txt"], shell = True)
    subprocess.run(["sed -i '/#------------------------------------------/d' pro.txt"], shell = True)
    subprocess.run(["sed -i '/#Resnum  Resnam     PHI      PSI     OMEGA/d' pro.txt"], shell = True)

def identifier(culled):
    if culled[0].isdigit() == True:
        #base_url = "https://files.rcsb.org/download/"
        pdb_str = str(culled)
        pdb = pdb_str[0:4]
        chain = pdb_str[4]
        #url = base_url + pdb + ".pdb"
        print(pdb)
        #fetcher(url, pdb) #commented out as I was kindly given copies of all desired chains
        try:
            commander(pdb, chain)
        except Exception as e:
            print("Exception error")
        except IOError as e:
            print("IOError")
    else:
        pass

def windbuffer(n):
    n = n//2
    f = open("Fourecks.txt", 'a')
    for i in range(n):
        f.write("X000     XXX       999.999 999.999  999.999\n")

def coma(letter):
    l = list(letter)
    line = []
    for res in l:
        line.append(res)
    line = ','.join(map(str, line))
    return(line)

def city(line, windlen):
    fir = (5-22+(45*windlen//2))
    las = (8-22+(45*windlen//2))
    ctra = (31 + (45*windlen//2))
    ctrb = (39 + (45*windlen//2))
    first = line[:fir]
    last = line[las:]
    mid = line[fir:las]
    ct = line[ctra:ctrb]
    if abs(float(ct)) <= 10:
        mid = 'CIS'
        line = (first + mid + last)
    else:
        mid = 'TNS'
        line = (first + mid + last)
    return(line)

def steppe(mylist, windlen):
    st = ""
    nitems = len(mylist)
    for i in range(nitems):
        st += mylist[i]
    st += "---      ---\n"
    return(city(st, windlen))

def window(mylist, item, windlen):
    nitems = len(mylist)
    pos = int(windlen//2)
    for i in range(nitems-1):
        mylist[i] = mylist[i+1]
    mylist[nitems-1] = item
    if "PRO" in mylist[pos]:
        return steppe(mylist, windlen)
    else:
        return("")

def csvinit(windlen):
    with open("csvbase.csv", 'w')as newcsv:        
        window=''
        for i in range(windlen//2):
            window += "res" + str(i+1) + ','
        window += 'bond,'
        for i in range((windlen//2)-1):
            window += "res" + (str(i+1+windlen//2)) + ','
        window += "res" + str(windlen) + ",\n"
        newcsv.write(window)
    newcsv.close()

def arff(windlen):            
    window=''
    for i in range(windlen//2):
        window += "res" + str(i+1) + ','
    for i in range((windlen//2)):
        window += "res" + (str(i+1+windlen//2)) + ','
    window += "res" + str(windlen) + " "
    window += 'bond'
    return(window)

def throne(three):
    aa = {
    'ALA': 'A', 
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLY': 'G',
    'GLN': 'Q',
    'GLU': 'E',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'W',
    'TYR': 'Y',
    'VAL': 'V',
    'XXX': 'X',
    'CIS': '3',
    'TNS': '5',
    '---': "\n"}
    if len(three) == 3:
        try:
            return(aa[three])
        except KeyError as e:
            return('O') #non standard >> other
    else:
        return("")

def proppy(three):
    prop = {
    'ALA': 'H', #GAS
    'ARG': 'C', #HRK
    'ASN': 'B', #MQPCTN
    'ASP': 'D', #DE
    'CYS': 'B',
    'GLY': 'H',
    'GLN': 'B',
    'GLU': 'D',
    'HIS': 'C',
    'ILE': 'A', #first group according to paper by AM, seperation based on venn
    'LEU': 'A',
    'LYS': 'C',
    'MET': 'B', 
    'PHE': 'R', #ring
    'PRO': 'B',
    'SER': 'H',
    'THR': 'B',
    'TRP': 'R',
    'TYR': 'R',
    'VAL': 'A',
    'XXX': 'X',
    'CIS': '3',
    'TNS': '5',
    '---': "\n"} #go find willie taylor paper 1987 paper amino acid properties
    if len(three) == 3:
        try:
            return(prop[three])
        except KeyError as e:
            return('O')#non-standard gets marked as 'O', 'other'
    else:
        return("")

def splitter(redden, writea, writeb):
    for line in redden:
        writea.write(coma(str(throne(line))))
        writeb.write(coma(str(proppy(line))))


def ine(f,g):
    tick = 0
    for line in f:
        tick += 1
        if tick > 1:
            line = line.replace('5', 'False')
            line = line.replace('3', 'True')
        else:
            pass
        g.write(line)

'''
setting variables
_________________________________________________________________________________
'''   
windlen = 9 #has to be odd and integer for anticipated results
'''
code part
_________________________________________________________________________________
'''

#create fourecks.txt, to give a buffer in case of near-end PRO
windbuffer(windlen)
with open("culled.txt", 'r') as fi:
    f = fi.read().splitlines()
tick = 0
#read in culled.txt, identify protein sequences and chain
    #download said protein
    #split download into chains
    #get torsion angles and specifically those about PRO by windlen//2
    #delet downloaded files
    #add another buffer to the end
    #remove grep's "--" output
for line in f:
    identifier(line)

protidy()

#move a window through the above output, when the half+1th is a PRO then replaces it with bond type
f = open("pro.txt", 'r')
thing = open("thing.txt", 'w+')
l = []
for i in range(9):
    l += "i"
for line in f:
    thing.write(window(l, line, windlen))
f.close()

f = open("thing.txt",'r')
tmp = open("tmp.txt", 'w')
for line in f:
    cols = re.split(r'\s+', line)
    yy = str(cols[1] + '\n')
    tmp.write(yy)
f.close()
tmp.close()

csvinit(windlen)
naive = open("res_cod.txt", 'w+')
proper = open("res_trait.txt", 'w+')
with open("tmp.txt", 'r') as f:
    f = f.read().split()

splitter(f, naive, proper)

naive.close()
proper.close()

#csv
csvinit(windlen)
nv = open("naive.csv", 'w+')
tr = open("trait.csv", 'w+')
csv = open("csvbase.csv", 'r')
f = open("res_cod.txt", 'r')
g = open("res_trait.txt", 'r')


for line in csv:
    nv.write(line)
    tr.write(line)

for line in f:
    ph = coma(line)
    nv.write(ph)

for line in g:
    hp = coma(line)
    tr.write(hp)

f.close()
csv.close()
nv.close()
g.close()

tr.close()
sys_str_9 =("sed 's/.$//' naive.csv")
sys_str_10 =("sed 's/.$//' trait.csv")
subprocess.run([sys_str_9], shell = True)
subprocess.run([sys_str_10], shell = True)

g = open("na.csv", 'w+')
f =  open("naive.csv", 'r')
ine(f,g)
g.close()
f.close()

g = open("pr.csv", 'w+')
f = open("trait.csv", 'r')
ine(f,g)
g.close()
f.close()

arffstr = arff(windlen)
subprocess.run(["csv2arff -inputs=" + arffstr + " na.csv > naive.arff"], shell = True)
subprocess.run(["csv2arff -inputs=" + arffstr + " pr.csv > trait.arff"], shell = True)

#############
#quick count#
#############
'''
f =  open("count_tcis.txt", 'r')
counter = 0
for line in f:
    counter += int(line)
#print("total cis: ", counter)
f.close()
g = open("count_total.txt", 'r')
countD = 0
for line in g:
    countD += int(line)
#print("total count: ", countD)
g.close()

h = open("pr.csv", 'r')
countPC = 0
countPT = 0
for line in h:
    if "True" in line:
        countPC += 1
    if "False" in line:
        countPT += 1
    else:
        pass
#print("proCis count: ", countPC)
#print("proTrans count: ", countPT)
h.close()

#print("nproCis count: ", (counter-countPC))
#print("nproTrans count: ", (countD-countPT-counter))
#print("relative incidences")
print("Pro Cis incidence:      ", (100*countPC/countD), "%")
print("Pro Trans incidence:    ", (100*countPT/countD), "%")
print("nonPro Cis incidence:   ", (100*(counter-countPC)/countD), "%")
print("nonpro Trans incidence: ", (100*(countD-countPT-counter)/countD), "%")
print("Of cis:                 ", (100*(countPC/(counter-countPC))), "%")
'''