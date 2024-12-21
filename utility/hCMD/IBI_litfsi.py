import numpy as np
from scipy.interpolate import CubicSpline
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.signal import savgol_filter
from scipy.signal.signaltools import wiener
from scipy.optimize import curve_fit
from scipy.optimize import fmin
Temp = 298.0
beta = 1.0/(np.float64(1.987204259e-3)*Temp)
exact_file_name = "../ref/run/Average_RDF.d"
smoothing_ratio = 1e-5
epsilon = 5.0
delta_r = 0.01
RDF = np.loadtxt("run/Average_RDF.d")
RDFexact = np.loadtxt(exact_file_name)

npair=15
cutoff=10.0
nrdf=480
pairtype=["C OW",
        "C HW",
        "N OW",
        "N HW",
        "O OW",
        "O HW",
        "F OW",
        "F HW",
        "S OW",
        "S HW",
        "Li OW",
        "Li HW",
        "OW OW",
        "OW HW",
        "HW HW"
          ]

nmols=2
nbonds=[4,1]
totbonds=sum(nbonds)
bonds=[
        ["harm",4],
        ["harm",6],
        ["harm",2],
        ["harm",2],
        ["quar",2]
        ]
nangs=[7,1]
totangs=sum(nangs)
angs=[
        ["harm",4],
        ["harm",2],
        ["harm",4],
        ["harm",6],
        ["harm",6],
        ["harm",2],
        ["harm",1],
        ["harm",1],
        ]

Correction = RDFexact*1.0
difference = RDFexact*1.0
for pair_to_correct in range(1,npair+1):
  gmax = np.amax([RDFexact[:,pair_to_correct],RDF[:,pair_to_correct]])
  print("gmax",pair_to_correct,gmax)
  Correction[:,pair_to_correct] = -(1.0/beta)*np.log((epsilon*gmax+RDFexact[:,pair_to_correct])/(epsilon*gmax+RDF[:,pair_to_correct]))
  difference[:,pair_to_correct] = RDFexact[:,pair_to_correct] - RDF[:,pair_to_correct]

np.savetxt("New_Correction_epsilon.dat",Correction)  
np.savetxt("difference.dat",difference)  

Previous_Potential = np.loadtxt("Previous_Potential.dat")


r_vals = np.linspace(0.0,Correction[-1,0],int(Correction[-1,0]/delta_r)+1)
New_Potential = np.zeros((int(Correction[-1,0]/delta_r)+1,(RDF[0,:].size)))
New_Potential[:,0] = r_vals*1.0

for i in np.arange(1,(RDF[0,:].size)):
  spl_old = interp1d(Previous_Potential[:,0],Previous_Potential[:,i],kind='linear',fill_value='extrapolate')
  spl_correc = interp1d(Correction[:,0],Correction[:,i],kind='linear',fill_value='extrapolate')  
  New_Potential[:,i] = spl_old(r_vals) + spl_correc(r_vals)
  spl = UnivariateSpline(New_Potential[:,0],New_Potential[:,i],k=3,s=smoothing_ratio)
  New_Potential[:,i] = spl(New_Potential[:,0])
  deriv = spl.derivative(n=1)
  f = open("Pair_"+str(i)+"_pot.table","w")
  f.write(pairtype[i-1])
  f.write("\n")
  for j in np.arange(1,np.size(New_Potential[:,0])):
    f.write(str(j)+" "+str(New_Potential[j,0])+" "+str(New_Potential[j,i])+" "+str(-deriv(New_Potential[j,0]))+"\n")
  f.close()

np.savetxt("New_Potential.dat",New_Potential)

New_Potential_Raw = New_Potential*1.0
def quartic(x,c0,x0,c2,c3,c4):
  return c0 + 0.5*c2*(x-x0)**2 + (1.0/3.0)*c3*(x-x0)**3 + 0.25*c4*(x-x0)**4

def quartic_mod(x,c0,x0,c2,c3,cmod):
  return c0 + c2*(x-x0)**2 + c3*(x-x0)**3 + c3**2*cmod/(2.0*c2)*(x-x0)**4

def quadratic(x,c0,x0,c2):
  return c0 + 0.5*c2*(x-x0)**2

print("Now computing the new bond parameters")
pair = 2
exact_bond_file = "../ref/run/Average_Intra.d" 

smoothing_ratio = 1e-8
Ex_Intrafull = np.loadtxt(exact_bond_file)
Ap_Intrafull = np.loadtxt("run/Average_Intra.d")

nrdf=480
Ex_Intra = np.zeros((nrdf,totbonds+1))
Ap_Intra = np.zeros((nrdf,totbonds+1))
inp = open("old_FIELD","r")
lines=inp.readlines()
fieldsplit=[]
for line in lines:
    fieldsplit.append(line.split())
f = open("new_FIELD","w")
fieldnew = fieldsplit

bondparamsfull=[]
bondparams=[]
bondparamsnew=[]
angparamsfull=[]
angparams=[]
angparamsnew=[]
linenum=0
bondline=[]
angline=[]
molbnd=[]
molang=[]
for line in fieldsplit:
    if line[0].lower() == "bonds":
        molbnd.append(int(line[1]))
        bondline.append(linenum+1)
    if line[0].lower() == "angles":
        molang.append(int(line[1]))
        angline.append(linenum+1)
    linenum = linenum+1

for i in np.arange(0,nmols):
    bndl = bondline[i]
    for j in np.arange(0,molbnd[i]):
        bondparamsfull.append(fieldsplit[bndl])
        bndl = bndl + 1
    bondparamsfull.append("")
    angl = angline[i]
    for j in np.arange(0,molang[i]):
        angparamsfull.append(fieldsplit[angl])
        angl = angl + 1
    angparamsfull.append("")
nbnd=0
ntyp=0

bondparamsfullnew=bondparamsfull
angparamsfullnew=angparamsfull
#for line in fieldsplit:
#    f.write(' '.join(line) + '\n')
for i in np.arange(1,nmols+1):
    for j in np.arange(1,nbonds[i-1]+1):
        bondparams.append(bondparamsfull[nbnd][3:7])
        nbnd=nbnd+bonds[j-1][1]
    nbnd=nbnd+1
for i in np.arange(0,totbonds):
    print("Bond number: ",i+1)
    Ex_Intra[:,0] = Ex_Intrafull[:,0]
    Ap_Intra[:,0] = Ap_Intrafull[:,0]
    Ex_Intra[:,1] = Ex_Intrafull[:,i+1]
    Ap_Intra[:,1] = Ap_Intrafull[:,i+1]
    
    Ap_Intra_pot = -(1.0/beta)*np.log(Ap_Intra)
    Ex_Intra_pot = -(1.0/beta)*np.log(Ex_Intra)
    Ap_Intra_pot[:,0] = Ap_Intra[:,0]*1.0
    Ex_Intra_pot[:,0] = Ex_Intra[:,0]*1.0
     
    Ap_Intra_pot_trim = Ap_Intra_pot[~np.isnan(Ap_Intra_pot[:,1])]
    Ap_Intra_pot_trim = Ap_Intra_pot_trim[~np.isinf(Ap_Intra_pot_trim[:,1])]
    Ex_Intra_pot_trim = Ex_Intra_pot[~np.isnan(Ex_Intra_pot[:,1])]
    Ex_Intra_pot_trim = Ex_Intra_pot_trim[~np.isinf(Ex_Intra_pot_trim[:,1])]
    
    ex_spl = UnivariateSpline(Ex_Intra_pot_trim[:,0],Ex_Intra_pot_trim[:,1],s=smoothing_ratio,ext=1)
    ap_spl = UnivariateSpline(Ap_Intra_pot_trim[:,0],Ap_Intra_pot_trim[:,1],s=smoothing_ratio,ext=1)
    #ex_spl = UnivariateSpline(Ex_Intra_pot_trim[:,0],Ex_Intra_pot_trim[:,1],ext=1)
    #ap_spl = UnivariateSpline(Ap_Intra_pot_trim[:,0],Ap_Intra_pot_trim[:,1],ext=1)
    
    Ex_Intra_pot[:,1] = ex_spl(Ex_Intra_pot[:,0])
    Ap_Intra_pot[:,1] = ap_spl(Ap_Intra_pot[:,0])
    overlap = Ex_Intra[:,1] * Ap_Intra[:,1] / np.amax(np.concatenate((Ex_Intra[:,1],Ap_Intra[:,1])))
    fit_region_intra = overlap>np.exp(-10.0)
    fit_region_intra = Ex_Intra_pot[fit_region_intra,0]
    Intra_Correction = Ex_Intra_pot*1.0
    Intra_Correction[:,1] = Intra_Correction[:,1]-Ap_Intra_pot[:,1]
    New_Intra_Pot = Intra_Correction*1.0
    k2=float(bondparams[i][0])
    r0=float(bondparams[i][1])
    k3=float(bondparams[i][2])
    k4=float(bondparams[i][3])
    typ=bonds[i][0]
    print(i+1,typ,k2,r0,k3,k4)
    if(typ=="quar"):
        New_Intra_Pot[:,1] = 0.5*k2*(Intra_Correction[:,0]-r0)**2+(1.0/3.0)*k3*(Intra_Correction[:,0]-r0)**3 + 0.25*k4*(Intra_Correction[:,0]-r0)**4  
    elif(typ=="harm"):
        New_Intra_Pot[:,1] = 0.5*k2*(Intra_Correction[:,0]-r0)**2 
    spl_correc = CubicSpline(Intra_Correction[:,0],Intra_Correction[:,1])
    New_Intra_Pot[:,1] = New_Intra_Pot[:,1] + spl_correc(New_Intra_Pot[:,0])
    spl = CubicSpline(New_Intra_Pot[:,0], New_Intra_Pot[:,1])
    deriv = spl.derivative
    minimum = fmin(spl,r0)
    rmax = np.sqrt(3.0/(beta*k2))
    rs = np.linspace(fit_region_intra[0],fit_region_intra[-1],50)
    Fitting_Pot = spl(rs)
    if(typ=="quar"):
        coefs, covar = curve_fit(quartic, rs, Fitting_Pot,p0=(0,r0,k2,k3,k4),maxfev=5000)
        if(i==totbonds-1):
            k2=coefs[2]
            r0=coefs[1]
            k3=coefs[3]
            k4=coefs[4]
        print(fit_region_intra[0],fit_region_intra[-1])
        print("New k2,r0,k3,k4",coefs[2],coefs[1],coefs[3],coefs[4])
    elif(typ=="harm"):
        coefs, covar = curve_fit(quadratic, rs, Fitting_Pot,p0=(0,r0,k2),maxfev=5000)
        if(i==totbonds):
            k2=coefs[2]
            r0=coefs[1]
        k3=0.0
        k4=0.0
        print(fit_region_intra[0],fit_region_intra[-1])
        print("New k2,r0",coefs[2],coefs[1])
    #New_Intra_Raw = 1.0*New_Intra_Pot
    #New_Intra_Pot[:,1] = quartic(New_Intra_Pot[:,0],0.0,coefs[1],coefs[2],coefs[3],coefs[4])  
    bondparamsnew.append([k2,r0,k3,k4])
bnd=0
pbnd=0
for i in np.arange(0,nmols):
    for j in np.arange(0,nbonds[i]):
        if(i>0):
            bnds=bonds[j+nbonds[i-1]][1]
        else:
            bnds=bonds[j][1]
        for k in np.arange(0,bnds):
            bondparamsfullnew[bnd][3:7]=bondparamsnew[pbnd][:]
#            for l in np.arange(0,7):
#                f.write(str(paramsfullnew[bnd][l])+" ")
#            f.write("\n")
            bnd = bnd+1
        pbnd=pbnd+1
#    f.write("\n")
    bnd=bnd+1
#
#f.close()
#
print("Now computing the new angle parameters")
exact_angle_file = "../ref/run/Average_Angle.d" 

fp = open("angle_correct.dat","w")
smoothing_ratio = 1e-8
Ex_Angfull = np.loadtxt(exact_angle_file)
Ap_Angfull = np.loadtxt("run/Average_Angle.d")

Ex_Ang = np.zeros((nrdf,2))
Ap_Ang = np.zeros((nrdf,2))

nang=0
ntyp=0
for i in np.arange(1,nmols+1):
    for j in np.arange(1,nangs[i-1]+1):
        angparams.append(angparamsfull[nang][4:7])
        nang=nang+angs[j-1][1]
    nang=nang+1
for i in np.arange(0,totangs):
    print("Angle number: ",i+1)
    Ex_Ang[:,0] = Ex_Angfull[:,0]
    Ap_Ang[:,0] = Ap_Angfull[:,0]
    Ex_Ang[:,1] = Ex_Angfull[:,i+1]
    Ap_Ang[:,1] = Ap_Angfull[:,i+1]
   
    #print(Ex_Ang[280:290,:])
    Ap_Ang_pot = -(1.0/beta)*np.log(Ap_Ang)
    Ex_Ang_pot = -(1.0/beta)*np.log(Ex_Ang)
    Ap_Ang_pot[:,0] = Ap_Ang[:,0]*1.0
    Ex_Ang_pot[:,0] = Ex_Ang[:,0]*1.0
    #print(Ex_Ang_pot[280:290,:])
     
    Ap_Ang_pot_trim = Ap_Ang_pot[~np.isnan(Ap_Ang_pot[:,1])]
    Ap_Ang_pot_trim = Ap_Ang_pot_trim[~np.isinf(Ap_Ang_pot_trim[:,1])]
    Ex_Ang_pot_trim = Ex_Ang_pot[~np.isnan(Ex_Ang_pot[:,1])]
    Ex_Ang_pot_trim = Ex_Ang_pot_trim[~np.isinf(Ex_Ang_pot_trim[:,1])]
     
    #ex_spl = UnivariateSpline(Ex_Ang_pot_trim[:,0],Ex_Ang_pot_trim[:,1],s=smoothing_ratio,ext=3)
    #ap_spl = UnivariateSpline(Ap_Ang_pot_trim[:,0],Ap_Ang_pot_trim[:,1],s=smoothing_ratio,ext=3)
    ex_spl = UnivariateSpline(Ex_Ang_pot_trim[:,0],Ex_Ang_pot_trim[:,1],ext=3)
    ap_spl = UnivariateSpline(Ap_Ang_pot_trim[:,0],Ap_Ang_pot_trim[:,1],ext=3)
    
    #print(Ex_Ang_pot_trim[:,:])
    #print(ex_spl.get_coeffs())
    Ex_Ang_pot[:,1] = ex_spl(Ex_Ang_pot[:,0])
    Ap_Ang_pot[:,1] = ap_spl(Ap_Ang_pot[:,0])
    #print(Ex_Ang_pot[280:310,:])
    overlap = Ex_Ang[:,1] * Ap_Ang[:,1] / np.amax(np.concatenate((Ex_Ang[:,1],Ap_Ang[:,1])))
    fit_region_Ang = overlap>np.exp(-10.0)
    fit_region_Ang = Ex_Ang_pot[fit_region_Ang,0]
    Ang_Correction = Ex_Ang_pot*1.0
    Ang_Correction[:,1] = Ang_Correction[:,1]-Ap_Ang_pot[:,1]
    New_Angle_Pot = Ang_Correction*1.0
    k2=float(angparams[i][0])
    theta0=float(angparams[i][1])
    theta0=theta0*np.pi/180.0
    typ=angs[i][0]
    print(i+1,typ,k2,theta0*180.0/np.pi)
    New_Angle_Pot[:,1] = 0.5*k2*(Ang_Correction[:,0]-theta0)**2 
    spl_correc = CubicSpline(Ang_Correction[:,0],Ang_Correction[:,1])
    New_Angle_Pot[:,1] = New_Angle_Pot[:,1] + spl_correc(New_Angle_Pot[:,0])
    spl = CubicSpline(New_Angle_Pot[:,0], New_Angle_Pot[:,1])
    deriv = spl.derivative
    minimum = fmin(spl,r0)
    thetamax = np.sqrt(3.0/(beta*k2))
    thetas = np.linspace(fit_region_Ang[0],fit_region_Ang[-1],50)
    Fitting_Pot = spl(thetas)
    deltav = spl_correc(thetas)
    ex_theta = ex_spl(thetas)
    ap_theta = ap_spl(thetas)
    for j in np.arange(0,50):
        fp.write(str(thetas[j])+" "+str(Fitting_Pot[j])+" "+str(deltav[j])+" "+str(ex_theta[j])+" "+str(ap_theta[j])+"\n")
    fp.write("\n\n")
    coefs, covar = curve_fit(quadratic, thetas, Fitting_Pot,p0=(0,theta0,k2),maxfev=5000)
    k2old=k2
    thetaold=theta0*180.0/np.pi
    if(i==totangs-1):
        k2=coefs[2]
        theta0=coefs[1]
    theta0=theta0*180.0/np.pi
    k3=0.0
    k4=0.0
    print(fit_region_Ang[0],fit_region_Ang[-1])
    print("New k2,theta0",coefs[2],180.0*coefs[1]/np.pi)
    print("Percent Difference k2: ",abs((k2-k2old)/((k2+k2old)/2.0))*100.0," theta0: ",abs((theta0-thetaold)/((theta0+thetaold)/2.0))*100.0)
    print(theta0,thetaold)
    #New_Intra_Raw = 1.0*New_Intra_Pot
    #New_Intra_Pot[:,1] = quartic(New_Intra_Pot[:,0],0.0,coefs[1],coefs[2],coefs[3],coefs[4])  
    angparamsnew.append([k2,theta0,k3,k4])
fp.close()

ang=0
pang=0
for i in np.arange(0,nmols):
    for j in np.arange(0,nangs[i]):
        if(i>0):
            angls=angs[j+nangs[i-1]][1]
        else:
            angls=angs[j][1]
        for k in np.arange(0,angls):
            angparamsfullnew[ang][4:7]=angparamsnew[pang][:]
#            for l in np.arange(0,8):
#                f.write(str(paramsfullnew[ang][l])+" ")
#            f.write("\n")
            ang = ang+1
        pang=pang+1
#    f.write("\n")
    ang=ang+1
#print(angparamsfullnew[:][:])
for i in np.arange(0,nmols):
    bndl = bondline[i]
    if i==0:
        jmin=0
        jmax=molbnd[i]
    else:
        jmin=jmax+1
        jmax=jmin+molbnd[i]
    for j in np.arange(jmin,jmax):
        fieldsplit[bndl][3]=str(bondparamsfullnew[j][3])
        fieldsplit[bndl][4]=str(bondparamsfullnew[j][4])
        fieldsplit[bndl][5]=str(bondparamsfullnew[j][5])
        fieldsplit[bndl][6]=str(bondparamsfullnew[j][6])
        #bondparamsfull.append(fieldsplit[bndl])
        bndl = bndl + 1
    angl = angline[i]
    if i==0:
        kmin=0
        kmax=molang[i]
    else:
        kmin=kmax+1
        kmax=kmin+molang[i]
    for k in np.arange(kmin,kmax):
        fieldsplit[angl][4]=str(angparamsfullnew[k][4])
        fieldsplit[angl][5]=str(angparamsfullnew[k][5])
        fieldsplit[angl][6]=str(angparamsfullnew[k][6])
        fieldsplit[angl][7]=str(angparamsfullnew[k][7])
        #bondparamsfull.append(fieldsplit[bndl])
        angl = angl + 1
    #bondparamsfull.append("")
#    angl = angline[i]
#    for j in np.arange(0,molang[i]):
#        angparamsfull.append(fieldsplit[angl])
#        angl = angl + 1
#    angparamsfull.append("")

for line in fieldsplit:
    f.write(' '.join(line) + '\n')
f.close()

exit()

