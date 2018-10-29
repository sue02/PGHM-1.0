#!/usr/bin/env python
import sys
import numpy as np
import os
import operator
import matplotlib.pyplot as plt
from matplotlib import *
import csv
from matplotlib.patches import Wedge
from pylab import *
import pylab
from matplotlib.projections import PolarAxes# The slices will be ordered and plotted counter-clockwise.
from matplotlib.colors import ListedColormap

import sys, getopt
 
ifile=''
samp1=''
samp2=''
ofile=''
 
###############################
# o == option
# a == argument passed to the o
###############################
# Cache an error with try..except 
# Note: options is the string of option letters that the script wants to recognize, with 
# options that require an argument followed by a colon (':') i.e. -i fileName
#
try:
    myopts, args = getopt.getopt(sys.argv[1:],"i:1:2:o:")
except getopt.GetoptError as e:
    print (str(e))
    print("Usage: %s -i input -1 sample1 -2 sample2 -o output" % sys.argv[0])
    sys.exit(2)
 
for o, a in myopts:
    #print o,a
    if o == '-i':
        ifile=a
    elif o == '-o':
        ofile=a
    elif o == '-1':
        samp1=a
    elif o == '-2':
        samp2=a
 
# Display input and output file name passed as the args
#print ("Input file : %s , Sample1 : %s, Sample2 : %s and output file: %s" % (ifile,samp1,samp2,ofile) )


inputfile = sys.argv[2]
y = sys.argv[4]
z = sys.argv[6]

outfile = 'corr.csv'
outfile1 = 'tmp'

nuctide= ['A','T','G','C']
difreq = []
trifreq = []
tetrafreq = []


# remove header from sequence file
myfile = open(inputfile+'_noHeader', 'w')
with open(inputfile) as reader:
    for line in reader:
        if line[0]!='>':
            myfile.write(str(line))
myfile.close()

# input file for calculating word frequecies
inputfile1 = inputfile+'_noHeader'

# caluclating di-tri-tetra word frequencies

def difrequency():
    for n in range(len(nuctide)):
        for i in range(4):
            difreq.append(nuctide[n]+nuctide[i])
    disort=sorted(difreq)
    return disort
def trifrequency():
    for n in range(len(nuctide)):
        for i in range(4):
            for j in range(4):
                trifreq.append(nuctide[n]+nuctide[i]+nuctide[j])
    trisort=sorted(trifreq)
    return trisort
def tetrafrequency():
    for n in range(len(nuctide)):
        for i in range(4):
            for j in range(4):
                for k in range(4):
                    tetrafreq.append(nuctide[n]+nuctide[i]+nuctide[j]+nuctide[k])
    tetrasort=sorted(tetrafreq)
    return tetrasort

di = difrequency()
tri = trifrequency()
tetra= tetrafrequency()
allFreq= di+tri+tetra

freqlst =[]
bigfreqlst = []
count1 = 0
diseqlst =[]
diseqlst1 =[]
triseqlst =[]
triseqlst1 =[]
tetraseqlst =[]
tetraseqlst1 =[]
di1 ={}
di2 ={}
tri1={}
tri2 ={}
tetra1={}
tetra2={}
ntidelst = []
bigntidelst = []
with open(inputfile1) as reader:
    for a in reader:
        a = a.strip()
        total_length = len(a)
        for x in di:
            for i in range(len(a)-1):
                if x ==a[i:i+2]:
                    diseqlst.append(a[i:i+2])
            for m in diseqlst:
                di1[m] = di1.get(m, 0) + 1
            diseqlst = []
            if x not in a:
                diseqlst1.append(x)
            for m in diseqlst1:
                di2[m] = di2.get(m, 0) +0
            diseqlst1 = []
        dAll1 = dict(di1, **di2)
        for k,v in sorted(dAll1.iteritems()):
            avg = round((float(v)/total_length),2)
            freqlst.append(avg)
            ntidelst.append(k)
        for x in tri:
            for i in range(len(a)-2):
                if a[i:i+3]==x:
                    triseqlst.append(a[i:i+3])
            for m in triseqlst:
                tri1[m] = tri1.get(m, 0) + 1
            triseqlst = []
            if x not in a:
                triseqlst1.append(x)
            for m in triseqlst1:
                tri2[m] = tri2.get(m, 0) +0
            triseqlst1 = []
        dAll2 = dict(tri1, **tri2)
        for k,v in sorted(dAll2.iteritems()):
            avg = (float(v)/total_length)
            ntidelst.append(k)
            freqlst.append(avg)
        for x in tetra:
            for i in range(len(a)-3):
                if a[i:i+4]==x:
                    tetraseqlst.append(a[i:i+4])
            for m in tetraseqlst:
                tetra1[m] = tetra1.get(m, 0) + 1
            tetraseqlst = []
            if x not in a:
                tetraseqlst1.append(x)
            for m in tetraseqlst1:
                tetra2[m] = tetra2.get(m, 0) +0
            tetraseqlst1 = []
        dAll3 = dict(tetra1, **tetra2)
        for k,v in sorted(dAll3.iteritems()):
            avg = (float(v)/total_length)
            ntidelst.append(k)
            freqlst.append(avg)
        bigfreqlst.append(freqlst)
        bigntidelst.append(ntidelst)
        freqlst = []
        ntidelst = []
        di1 ={}
        di2 ={}
        tri1={}
        tri2 ={}
        tetra1={}
        tetra2={}
        dAll1 ={}
        dAll2 ={}
        dAll3 ={}
        dAll4 ={}

data = np.array(bigfreqlst)


# giving initials to the sequences 
def nameCoeff():
    namelst1 = []
    namelst2 = []
    with open(inputfile) as f:
        for line in f:
            if line[0]=='>' and y in line:
                line = line.strip().split('\t')
                nm1 = line[1][0]+'1'
                namelst1.append(nm1)
            if line[0]=='>' and z in line:
                line = line.strip().split('\t')
                nm2 = line[1][0]+'2'
                namelst2.append(nm2)

    nm = namelst1+namelst2
    return nm
finalName = sorted(set(nameCoeff()))

# correlation file 

def corrCoeff():
    data = np.array(bigfreqlst)
    coeff = np.corrcoef(data, rowvar=1, bias=0, ddof=None)
    coeffarray = np.array(coeff)
    half  = np.triu(coeffarray)
   
    header = nameCoeff()
    header1 = [""]+nameCoeff()
    with open(outfile,'wb') as f:
        writer = csv.writer(f)
        writer.writerow(header1)
        for row, data_row in zip(header,half):
            writer.writerow([row]+data_row.tolist())
    return coeff


corrCoeff()

# tmp file with to get the correlation pairs

def first():
    with open(outfile) as f:
        writefile = file(outfile1, 'w')
        data = f.readlines()
        tb = [x for x in data[0].strip().split(',')]
        writefile.writelines(data[0])
        for i in data[1:]:
            line = i.strip().split(',')
            for n in range(1,len(line)):
                writefile.writelines(line[n]+'_'+line[0]+tb[n]+',')
            writefile.writelines('\n')
first()

# header for each correlation pair

def header():
    newlst = []
    lst = sorted(set(nameCoeff()))
    for i in range(len(lst)):
        for j in range(len(lst)):
            newlst.append(lst[i]+lst[j])

    return newlst

# correlation values with names

def correlationLst():
    biglst = []
    newlst = []
    name =  header()
    with open(outfile1) as f:
        data = f.readline()
        data = f.readlines()
        for nm in name:
            for i in data:
                line = i.strip().split(',')[:-1]
                for item in range(len(line)):
                    if nm in line[item]:
                        lst = line[item]
                        newlst.append(lst)
        
            biglst.append(newlst)
            newlst = []
    return biglst    

# fisher transformation

def fisherZ():
    R = correlationLst()
    rbar = []
    Zlst = []
    ZbarLst = []
    i = 0
    while i< len(R):
        for val in R[i]:
            j = val.split('_')[0]
            
            if float(j)!=1.0:
                r1 = 1 + float(j)
                r2 = 1 - float(j)
                Z = 0.5*np.log(r1/r2)
                Zlst.append(Z)
            elif float(j)==1.0:
                r1 = 1 + 0.999999
                r2 = 1 - 0.999999
                Z = 0.5*np.log(r1/r2)
                Zlst.append(Z)
        ZbarLst.append(np.mean(Zlst))
        i+=1
        Zlst = []
    for Zbar in ZbarLst:
        z1 = np.exp(Zbar)
        z2 = np.exp(-Zbar)
        rMean = (z1-z2)/(z1+z2)
        rbar.append(rMean)
    return rbar


# renaming with upper and lower cases instead of 1 and 2

def namelst():
    newlst = []
    for i in header():
        if i[1]=='1' and i[3]=='1':
            newlst.append(i[0]+i[2])
        elif i[1]=='1' and i[3]=='2':
            newlst.append(i[0]+i[2].lower())
        elif i[1]=='2' and i[3]=='2':
            newlst.append(i[0].lower()+i[2].lower())
        elif i[1]=='2' and i[3]=='1':
            newlst.append(i[0].lower()+i[2])

    return newlst



# get the names for the inner circle
def innerCircle():
    lst = []
    lst2 = []
    for i in finalName:
        lst.append(i[0]+i[0].lower())
        lst2.append(i[0])

    inner = list(set(lst))
    newinner = sorted(list(set(lst2)))
    return inner, newinner

d ={}
d2 ={}
innerDict ={}

name = []
val = []
innerNames =[]
innerVal = []

#make dictionary of Z values and names

values = fisherZ()
newlst = namelst()
for i, j in zip(newlst, values):
    d[i]=j

#make inner circle dictionary
innerLst = innerCircle()[0]
for i in innerLst:
    for k,v in d.iteritems():
        if k == i:
            nm = k
            val = v
            innerDict[nm]=v
            innerNames.append(k)
            innerVal.append(v)

# arranging inner circle
labelInner  = [innerNames[1],innerNames[0],innerNames[3],innerNames[2]]
valuesInner = [innerVal[1],innerVal[0],innerVal[3],innerVal[2]]


#sort full dictionary and make new dictioanry with values>0
for key, value in sorted(d.items(), key=lambda (k,v): (v,k)):
    if value!=0.0:
        name = key
        val = value
        d2[name]=val

# remove innercircle vaues from new dictionary
for i in innerLst:
    del d2[i]

    
upperName = []
upperVal = []
lowerName =[]
lowerVal =[]

#sort dictionary and get top 8 values form each platform
for k,v in sorted(d2.items(), key=lambda (k,v): (v,k)):
    if k.isupper():
        upperName.append(k)
        upperVal.append(v)
    if k.islower():
        lowerName.append(k)
        lowerVal.append(v)


list1, list2 = zip(*sorted(zip(upperVal[-10:],upperName[-10:])))
list3, list4 = zip(*sorted(zip(lowerVal[-10:],lowerName[-10:])))
                   
 
valOuter = list(list1 + list3)
LabelOuter = list(list2 + list4)

# remove same letter combination
newLabelOuter = []
for i in LabelOuter:
    if i[0]!=i[1]:
        newLabelOuter.append(i)

#sort the labels for outer circle and remove same groups values i.e. 'AA' or 'BB' for outer circle

def outerCircle():        
    val = valOuter
    lst1 = [None]*8
    lst2 = [None]*8
    list1 = sorted(newLabelOuter)
    list2 = sorted(newLabelOuter)
    nm = innerCircle()[1]
    for i in list1:
        if i[0] == nm[0] or i[1] == nm[0]:
            lst1[0] = i
            lst1[1] = i[::-1]
        if i[0] == nm[1] or i[1] == nm[1]: 
            lst1[2] = i
            lst1[3] = i[::-1]

        if i[0] == nm[2] or i[1] == nm[2]:
            lst1[4] = i
            lst1[5] = i[::-1]
        if i[0] == nm[3] or i[1] == nm[3]:
            lst1[6] = i
            lst1[7] = i[::-1]
    
    for i in list2:
        if i[0] == nm[0].lower() or i[1] == nm[0].lower():
            lst2[0] = i
            lst2[1] = i[::-1]
        if i[0] == nm[1].lower() or i[1] == nm[1].lower(): 
            lst2[2] = i
            lst2[3] = i[::-1]

        if i[0] == nm[2].lower() or i[1] == nm[2].lower():
            lst2[4] = i
            lst2[5] = i[::-1]
        if i[0] == nm[3].lower() or i[1] == nm[3].lower():
            lst2[6] = i
            lst2[7] = i[::-1]
    l1 = sorted(list(set(lst1)))
    l2 = sorted(list(set(lst2)))
    
    newlst = l1+l2
    finalLst = [None] * 8


    for i in newlst:
        if i[0]== nm[0]:
            finalLst[0] = i
        if i[0]== nm[1]:
            finalLst[1] = i
        if i[0]== nm[2]:
            finalLst[2] = i
        if i[0]== nm[3]:
            finalLst[3] = i 
        if i[0]== nm[0].lower():
            finalLst[4] = i
        if i[0]== nm[1].lower():
            finalLst[5] = i
        if i[0]== nm[2].lower():
            finalLst[6] = i
        if i[0]== nm[3].lower():
            finalLst[7] = i
        

    return finalLst
finalnameLst = outerCircle()
# reaaraging values for outer circle

def arrangingValLabelOuterCircle():
    d = {}
    finalD = {}
    finalVal = []
    finalLabel = []
    LabelOuter = list(list2 + list4)
    for i, j in zip(LabelOuter, valOuter):
        d[i] = j
    for a in finalnameLst:
        for k, v in d.iteritems():
            if k==a or k==a[::-1]:
                finalVal.append(v)
                finalLabel.append(a)

        
    return finalVal, finalLabel

finalOuterVal = arrangingValLabelOuterCircle()[0]
finalOuterLabel = arrangingValLabelOuterCircle()[1]



newfinalOuterVal = [finalOuterVal[0],finalOuterVal[4],finalOuterVal[1],finalOuterVal[5],finalOuterVal[2],finalOuterVal[6],finalOuterVal[3],finalOuterVal[7]]
newfinalOuterLabel = [finalOuterLabel[0],finalOuterLabel[4],finalOuterLabel[1],finalOuterLabel[5],finalOuterLabel[2],finalOuterLabel[6],finalOuterLabel[3],finalOuterLabel[7]]

II = []
OO = []

for x in valuesInner:
    val = float(x)
    r = 1-(val**2)
    newR = -1*(np.log(r))
    II.append(newR)
for x in newfinalOuterVal:
    val = float(x)
    r = 1-(val**2)
    newR = -1*(np.log(r))
    OO.append(newR)


fullLst = II + OO
fullLst1 = valuesInner + newfinalOuterVal

valuesI = fullLst[:4]
valuesO = fullLst[4:]

newN = []
for i in fullLst:
    j = float(i)
    newN.append(j*1000)

maxN = round(max(newN))


fig = figure()
ax = fig.add_subplot(111, projection = 'polar')
subplot(111,projection='polar')
theta = arange(maxN)*2*pi/maxN


#plot circle

plot(theta,0.25*ones(maxN),'k') #A to B Circle
plot(theta, 0.65*ones(maxN), 'k') # B to C Circle
plot([0, 0],[0,0.65],'k') # A to B 0 degrees line
plot([pi/4.,pi/4.],[0.25, 0.65], 'k') # A to B 45 degrees line
plot([pi/2, pi/2.],[0, 0.65], 'k') # A to B 90 degrees line
plot([3*pi/4.,3*pi/4],[0.25, 0.65], 'k') # A to B 135 degrees line
plot([pi,pi],[0,0.65],'k') # A to B 180 degrees line
plot([5*pi/4,5*pi/4],[0.25, 0.65], 'k') # A to B 225 degrees line
plot([3*pi/2,3*pi/2],[0, 0.65], 'k') # A to B 270 degrees line
plot([7*pi/4,7*pi/4],[0.25, 0.65], 'k') # A to B 315 degrees line


# make color lst

minima = min(fullLst)
maxima = max(fullLst)

#print minima,maxima

norm = matplotlib.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.Greys_r)


I = [mapper.to_rgba(v)[0] for v in fullLst[:4]]

O = [mapper.to_rgba(v)[0] for v in fullLst[4:]]
newO = [O[6],O[7],O[0],O[1],O[2],O[3],O[4],O[5]]


halfwedge = np.pi/8

start = 0
end = 45

for i in range(8):
    if O[i]>=0.5:
        ax.add_artist(Wedge((0.5,0.5), 0.325, start,end, width=0.2, transform=ax.transAxes, color=str(newO[i])))
        ax.text(np.pi/0.25 - halfwedge*(2*i+1), 0.425,newfinalOuterLabel [i], fontsize =15,verticalalignment='center', horizontalalignment='center')
        #ax.text(np.pi/0.25 - halfwedge*(2*i+1), 0.555, round(O[i],3),fontsize =10,verticalalignment='bottom', horizontalalignment='center')
    if O[i]<=0.5:
        ax.add_artist(Wedge((0.5,0.5), 0.325, start,end, width=0.2, transform=ax.transAxes, color=str(newO[i])))
        ax.text(np.pi/0.25 - halfwedge*(2*i+1), 0.425,newfinalOuterLabel [i], fontsize =15,verticalalignment='center', horizontalalignment='center',color='w')
        #ax.text(np.pi/0.25 - halfwedge*(2*i+1), 0.555, round(O[i],3),fontsize =10,verticalalignment='bottom', horizontalalignment='center',color='w')
    start = start + 45
    end = end+45

start1 = 0
end1 = 90
for i in range(4):
    if I[i]>=0.5:
        ax.add_artist(Wedge((0.5,0.5), 0.125, start1,end1, width=0.125, transform=ax.transAxes, color = str(I[i])))
        ax.text(np.pi/2 - halfwedge*2*(2*i+1), 0.145, labelInner[i],fontsize =15,verticalalignment='center', horizontalalignment='center')
        #ax.text(np.pi/2 - halfwedge*2*(2*i+1), 0.210, round(I[i],3),fontsize =10,verticalalignment='bottom', horizontalalignment='center')
    if I[i]<0.5:
        ax.add_artist(Wedge((0.5,0.5), 0.125, start1,end1, width=0.125, transform=ax.transAxes, color = str(I[i])))
        ax.text(np.pi/2 - halfwedge*2*(2*i+1), 0.145, labelInner[i],fontsize =15,verticalalignment='center', horizontalalignment='center', color='w')
        #ax.text(np.pi/2 - halfwedge*2*(2*i+1), 0.210, round(I[i],3),fontsize =10,verticalalignment='bottom', horizontalalignment='center', color = 'w')

       
    start1 = start1 + 90
    end1 = end1 + 90


# remove grids and coordinated
ax.grid(False)
rgrids((1,1),('',''))
ax.spines['polar'].set_visible(False)
ax.xaxis.set_ticks([])
ax.set_theta_zero_location("N")
ax.set_theta_direction('clockwise')         

#show()

plt.savefig(sys.argv[8], dpi=300)
os.remove(outfile1)
os.remove(outfile)
os.remove(inputfile+'_noHeader')
