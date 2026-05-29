'''
Graph functions on complex domain by colour coding;
and search for roots by analysing the graphical results (2D)
'''
import sys
import time
print ("   WC =", time.perf_counter(), "  PT = ", time.process_time())
stimeWC = time.perf_counter()   # wall clock
stimePT = time.process_time()   # processor time

import numpy as np
import matplotlib.pyplot as plt
from mpmath import *
import cmath
import cpxmods
                                    #@@@@@ specify Low value, High value, and Resolution here...
lovalR = -6;  hivalR = 6 
lovalI = -6;  hivalI = 6
resol  =  401
thldG  =  0.5                       #@@@@@ graph threshold - fade if fz exceeds (zero = no fading)
thldS  =  0.25                      #@@@@@ solution threshold - disregard if fz exceeds
deltaR = (hivalR - lovalR) / (resol-1)
refnct = 4                          #@@@@@ number of Refining iterations (0 - 5)

solnDsep = abs((hivalR - lovalR) / 20)          # min separation distance for roots
#solnDsep = 0.5

#funtitle = "f(z) = z^pi + 1"                          #@@@@@ set title here...
#funtitle = "f(z) = z^3.5 + 1" 
#funtitle = "f(z) = z^2.9 + 1"
#funtitle = "f(z) = e^z - 1"   # 7^z + 2^z - 130"
funtitle = "f(z) = 2^z + 2^-z - 1"
#A = 2.7444
#funtitle = "f(z) = z/log(z) - " + str(A)
#funtitle = "f(z) = 4^z + 5^z - 6^z"
#funtitle = "f(z)=zeta(z)" 

def myFunc(z):                                          #@@@@@  put function here:  fz = f(z)
    fz = 2**z + 2**-z - 1
    #fz = zetaFun(z)
    #fz = np.power(z, np.pi) + 1
    return fz

print("Title :", funtitle) 
print("Range =", lovalR, hivalR, ";  resol =", resol, "; delta =", deltaR,
                                   "; thldG =", thldG, "; thldS =", thldS)
print()

domR = np.linspace(lovalR, hivalR, resol)    #real domain
domI = np.linspace(hivalI, lovalI, resol)    #imag domain

def RGBfun1(xR, xI):
    z1 = complex(xR, xI)
    z2 = myFunc(z1)                                     
    
    z2colr, Quad = RGBmap(z2)
    if thldG > 0 and abs(z2) > thldG:
        pass
        #z2colr = (1,1,1,1.0)           # blackout larger values of f(z)
    return z2, z2colr

def RGBmap(z):
    Quadrant = 0; QuadColr = "k"
    
    magZ = abs(z); magLo = 0.5; magHi = 2.0
    alpha = 0.1/(magZ/5.0+0.1)
    if thldG > 0 and magZ > thldG:
        alpha = alpha * 0.1             # fade larger values of f(z)
    if z.real == 0.0 and z.imag == 0.0:
        RGBval = (1,0,1,alpha)
        QuadColr = "w"
        pass
        #print ("@@@ Zero value @@@ ", z)
    else:
      if z.real < 0:
        if z.imag < 0:
            RGBval = (1,0,0,alpha)  # 3rd quad = red
            #Rval = shade
            Quadrant = 3
            QuadColr = "red"
            if magZ > magHi: QuadColr = "darkred"
            if magZ < magLo: QuadColr = "lightsalmon"
            QuadColr = "r"
        else:
            RGBval = (1,1,0,alpha)  # 2nd quad = yellow
            #Rval = shade
            #Gval = shade
            Quadrant = 2
            QuadColr = "gold"
            if magZ > magHi: QuadColr = "orange"
            if magZ < magLo: QuadColr = "yellow"
            QuadColr = "y"
      else:
        if z.imag < 0:
            RGBval = (0,1,0,alpha)  # 4th quad = green
            #Gval = shade
            Quadrant = 4
            QuadColr = "limegreen"
            if magZ > magHi: QuadColr = "darkgreen"
            if magZ < magLo: QuadColr = "lime"
            QuadColr = "g"
        else:
            RGBval = (0,0,1,alpha)  # 1st quad = blue
            #Bval = shade
            Quadrant = 1
            QuadColr = "blue"
            if magZ > magHi: QuadColr = "darkblue"
            if magZ < magLo: QuadColr = "aqua"
            QuadColr = "b"
        
    return RGBval, Quadrant

def zetaFun(z):
    #print ("z =", z)
    if z.real == 1:
        #print ("   z.real = 1:", z) 
        z2 = complex(1000000000, 0)
        #zetaZ = altzeta(z) 
    else:  
        zetaZ = zeta(z)        # error for real z = 1
        #zetaZ = altzeta(z)    # OK for all z
        z2 = complex(zetaZ) 
    return (z2);
    
def natLog(z):
    if z.real == 0 and z.imag == 0:
        z2 = -1000 - 1000j
    else:
        z2 = cmath.log(z)
    return (z2)

zedR = np.zeros ((resol*resol), dtype=np.float32)   # array of z.real values
zedI = np.zeros ((resol*resol), dtype=np.float32)   # array of z.imag values
#fClr = np.zeros ((resol*resol), dtype=np.float32)   # array of f(z) colours
fzColr = []     # list of f(z) colours

plt.close('all') 
fig = plt.figure()
ax = fig.add_subplot(111)
plt.subplots_adjust(right=0.75)
#ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('z.Real')
ax.set_ylabel('z.Imag')
#ax.set_zlabel('f(z)')

solnList = []       # list of possible solutions (z1,z2,|z2|)

print("Solution proximity >", f"{solnDsep:.4f}")

I = J = 0
idx = 0
for xI in domI:
  J = 0
  for xR in domR:
    z1 = complex(xR, xI)
    z2, z2colr = RGBfun1(xR, xI)
    zedR[idx] = xR
    zedI[idx] = xI
    fzColr.append(z2colr) 
    idx += 1
    
    z2abs = abs(z2)
    if I == 0 and J == 0:
        minZ1 = z1
        minZ2 = z2
        minZ2a = z2abs
        minI = I; minJ = J
        #print ("..min at", I, J, ":", z1, "     ~", z2, z2abs)
    else:
        if z2abs < minZ2a:
            minZ1 = z1
            minZ2 = z2
            minZ2a = z2abs
            minI = I; minJ = J
            #print ("  min at", I, J, ":", z1, "     ~", z2, z2abs)

    if z2abs < thldS:                               # within solution threshold
        newsoln = True
        #print("   possible new soln at", z1, z2)
        for solidx, soln in enumerate(solnList):
                #print ("      soln#", solidx, soln, len(solnList) )
            soln_z1 = soln[0]; soln_z2 = soln[1]; soln_z2abs = soln[2] 
            if abs(z1 - soln_z1) < solnDsep:        # close to an existing soln?
                newsoln = False
                if z2abs < soln_z2abs:              # a better soln?
                    solupd = [z1, z2, z2abs]
                    solnList[solidx] = solupd       # then update it
                    #print("     update soln", solidx, f"{z1:.4f}", f"{z2:.4f}", f"{z2abs:.4f}")
                break 
        if newsoln and len(solnList) < 10:
            solnew = [z1, z2, z2abs]
            solnList.append(solnew)
            print("       new soln at", f"{z1:.4f}", f" {z2:.4f}", f" {z2abs:.4f}")        
    
    J += 1
  I += 1

print()
print (len(solnList), "solutions before culling...")
for idx, soln in enumerate(solnList):
    print ("      soln#", idx, f"  {soln[0]:.4f}", f"   {soln[1]:.4f}", f"   {soln[2]:.4f}") 
for idxA, solnA in enumerate (solnList):
    idxB = len(solnList)
    for solnB in reversed(solnList):
        idxB -= 1
        if idxB > idxA:
            solnA_z1 = solnA[0]; solnA_z2 = solnA[1]; solnA_z2abs = solnA[2]
            solnB_z1 = solnB[0]; solnB_z2 = solnB[1]; solnB_z2abs = solnB[2]
            if abs(solnA_z1 - solnB_z1) < solnDsep:         # solutions close by?
                print("  duplicates...", idxA, idxB) 
                del solnList[idxB]
                print("     delete soln", idxB)           

print()
print(len(solnList), "after 1st cull...    z               fz             |fz|")
for solidx, soln in enumerate(solnList):
    print ("      soln#", solidx, f"  {soln[0]:.4f}", f"   {soln[1]:.4f}", f"   {soln[2]:.4f}") 
print()

refnReach = deltaR                          # refine solutions (to minimise |fz|)
print("Refine solutions...   Reach =", f"{refnReach:f}")
for rct in range(refnct):
    #print("Refine solutions...   Reach =", f"{refnReach:f}")
    print("  refloop", rct+1)
    for solidx, soln in enumerate(solnList):
        zVal = soln[0]
        refSoln = cpxmods.cpxExtractSoln(myFunc, zVal, refnReach)
        if refSoln != zVal:
            soln[0] = refSoln
            soln[1] = myFunc(refSoln)
            soln[2] = abs(soln[1])
            #print ("    repl soln#", solidx, f"  {soln[0]:f}", f"   {soln[1]:f}", f"   {soln[2]:f}") 
    for solidx, soln in enumerate(solnList):
        print ("      soln#", solidx, f"  {soln[0]:f}", f"   {soln[1]:f}", f"   {soln[2]:f}") 
    #print()
    refnReach /= 10.0

print()
solidx = len(solnList)
for soln in reversed(solnList):
    solidx -= 1 
#for solidx, soln in enumerate(solnList):
    #print ("      soln#", solidx, soln, len(solnList) )
    soln_z1 = soln[0]; soln_z2 = soln[1]; soln_z2abs = soln[2]
    if (soln_z2abs > thldS / 20.0 
            or soln_z1.real < lovalR or soln_z1.real > hivalR
            or soln_z1.imag < lovalI or soln_z1.imag > hivalI): 
        print("  remove soln#", solidx)
        del solnList[solidx]
print(len(solnList), "after 2nd cull...    z               fz             |fz|")
for solidx, soln in enumerate(solnList):
    print ("      soln#", solidx, f"  {soln[0]:.4f}", f"   {soln[1]:.4f}", f"   {soln[2]:.4f}") 
print()

print ("   zedR:", zedR.dtype, zedR.shape, zedR.size)
#print (zedR)
print ("   zedI:", zedI.dtype, zedI.shape, zedI.size)
#print (zedI) 
print (" fzColr:", len(fzColr), sys.getsizeof(fzColr), "; max =", max(fzColr))
#print (fzColr)

ax.scatter(zedR, zedI, c=fzColr, s=300/resol)
ax.set_aspect('equal')

print ("minimum z1:", minZ1, abs(minZ1), ";\n        z2:", minZ2, abs(minZ2))
print ("    at I,J:", minI, minJ)
print ("         Z:", domR[minJ], domI[minI])

print ("   WC =", time.perf_counter(), "  PT = ", time.process_time())
etimeWC = time.perf_counter()   # wall clock
etimePT = time.process_time()   # processor time
print ("Durtn WC =", (etimeWC - stimeWC))
print ("Durtn PT =", (etimePT - stimePT))

ax.set_title(funtitle)                      # set title above
ax.set_facecolor('whitesmoke')
plt.grid(True)

if thldG > 0:
    ax.text(lovalR, lovalI, "Threshold=" + str(round(thldG,3)), size='small', style='italic')
plt.gcf().text(0.05, 0.01, "Resolution = " + str(resol), fontsize=6)
plt.gcf().text(0.80, 0.80, "Minimum at...", fontsize=8)
plt.gcf().text(0.84, 0.77, str(' zR = {0:.4f}'.format(domR[minJ])), fontsize=8)
plt.gcf().text(0.85, 0.74, str(' zI = {0:.4f}'.format(domI[minI])), fontsize=8)
plt.gcf().text(0.82, 0.71, str(' |f(z)| = {:.2e}'.format(abs(minZ2))), fontsize=8)

plt.gcf().text(0.76, 0.56, str('{} solutions detected...'.format(len(solnList))), fontsize=8)
for solidx, soln in enumerate(solnList):
    solx = soln[0]; fx= soln[1]; fxAbs=soln[2]
    solxR = solx.real; solxI = solx.imag
    solxFmt = ('{:+.4f}'.format(solx))
    if fxAbs < 0.00001:
        solnQual = ">5"
    else:
        solxQ = -np.log10(fxAbs)           # Quality of the solution - inv log scale
        solnQual = ('{:.1f}'.format(solxQ))
    plt.gcf().text(0.78, 0.52-0.04*solidx,
                   str('{})  {}   |{}|'.format(solidx+1, solxFmt, solnQual)), fontsize=7)
            #str('{0:}) {1:.3f}+{2:.3f}j {3:.1}'.format(solidx+1, solxR, solxI, solxQ)), fontsize=7)

#np.set_printoptions(precision=4)
#plt.gcf().text(0.85, 0.60, " zI = "+str(domI[minI]), fontsize=8)
#plt.gcf().text(0.82, 0.45, " |f(z)| = "+str(abs(minZ2)), fontsize=8)
#plt.gcf().text(0.82, 0.45, " |f(z)| = "+str('{:.2e}'.format(abs(minZ2))), fontsize=8)

#plt.subplots_adjust(right=0.9)
#plt.imshow(fun1array, extent=(loval,hival,loval,hival))     # matplotlib.pyplot.imshow displays RGB array
plt.show()
