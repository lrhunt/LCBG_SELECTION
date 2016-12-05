import numpy as np
import kcorrect
import astropy as ap
import matplotlib.pyplot as plt
import kcorrect
import kcorrect.utils as ut
import numpy as np
from astropy.cosmology import WMAP9 as cosmo

filein='CHILESSOURCELIST.txt'


IDNUM,RA,DECL,b06,b06err,v06,v06err,r06,r06err,i06,iflag,iSerr,iCerr,z06,z06err,zbest,rh,stellarmass,zuse=np.loadtxt(filein,skiprows=1,unpack=True)

wheresub=np.where(iflag==1)[0]
wherecfht=np.where(iflag==2)[0]

b06=b06+0.189
v06=v06+0.04
r06=r06-0.04
i06[wheresub]=i06[wheresub]-0.02
z06=z06+0.054

wheresub=np.where(iflag==1)[0]
wherecfht=np.where(iflag==2)[0]

bmaggies=ut.mag2maggies(b06)
vmaggies=ut.mag2maggies(r06)
rmaggies=ut.mag2maggies(v06)
imaggies=ut.mag2maggies(i06)
zmaggies=ut.mag2maggies(z06)

binvervar=ut.invariance(bmaggies,b06err)
vinvervar=ut.invariance(vmaggies,v06err)
rinvervar=ut.invariance(rmaggies,r06err)
iSinvervar=ut.invariance(imaggies,iSerr)
iCinvervar=ut.invariance(imaggies,iCerr)
zinvervar=ut.invariance(zmaggies,z06err)

acfht=np.stack((zbest,bmaggies,vmaggies,rmaggies,imaggies,zmaggies,binvervar,vinvervar,rinvervar,iCinvervar,zinvervar),axis=-1)
asub=np.stack((zbest,bmaggies,vmaggies,rmaggies,imaggies,zmaggies,binvervar,vinvervar,rinvervar,iSinvervar,zinvervar),axis=-1)
carr=np.ndarray((len(bmaggies),6))
carrsub=np.ndarray((len(bmaggies),6))
carrcfht=np.ndarray((len(bmaggies),6))
rmarr=np.ndarray((len(bmaggies),5))
rmarrsub=np.ndarray((len(bmaggies),5))
rmarrcfht=np.ndarray((len(bmaggies),5))
rmarr0=np.ndarray((len(bmaggies),5))

kcorrect.load_templates()
kcorrect.load_filters('/home/lhunt/programs/kcorrect/data/templates/LumFuncFiltersS.dat')


for i in range(0,len(asub)):
     c=kcorrect.fit_coeffs(asub[i])
     for j in range(0,len(carrsub[i])):
          carrsub[i][j]=c[j]

for i in range(0,len(carrsub)):
     rm=kcorrect.reconstruct_maggies(carrsub[i])
     for j in range(0,len(rmarrsub[i])):
          rmarrsub[i][j]=rm[j]


kcorrect.load_templates()
kcorrect.load_filters('/home/lhunt/programs/kcorrect/data/templates/LumFuncFiltersC.dat')


for i in range(0,len(wherecfht)):
     c=kcorrect.fit_coeffs(acfht[i])
     for j in range(0,len(carrcfht[i])):
          carrcfht[i][j]=c[j]

for i in range(0,len(wherecfht)):
     rm=kcorrect.reconstruct_maggies(carrcfht[i])
     for j in range(0,len(rmarr[i])):
          rmarrcfht[i][j]=rm[j]

carr[wherecfht]=carrcfht[wherecfht]
carr[wheresub]=carrsub[wheresub]
rmarr[wherecfht]=rmarrcfht[wherecfht]
rmarr[wheresub]=rmarrsub[wheresub]

kcorrect.load_templates()
kcorrect.load_filters('/home/lhunt/programs/kcorrect/data/templates/LumFuncFiltersB.dat')

for i in range(0,len(rmarr)):
     rm0=kcorrect.reconstruct_maggies(carr[i],redshift=0)
     for j in range(0,len(rmarr[i])):
          rmarr0[i][j]=rm0[j]

kcorr=-2.5*np.log10(rmarr/rmarr0)

M=np.zeros_like(b06)

for i in range(0,len(zbest)):
     if zbest[i]<=0.1:
          M[i]=b06[i]-5*(np.log10(100000*cosmo.luminosity_distance(zbest[i]).value))-kcorr[i][1]
     if zbest[i]<=0.32 and zbest[i]>0.1:
          M[i]=v06[i]-5*(np.log10(100000*cosmo.luminosity_distance(zbest[i]).value))-kcorr[i][2]
     if zbest[i]<=0.5 and zbest[i]>0.32:
          M[i]=r06[i]-5*(np.log10(100000*cosmo.luminosity_distance(zbest[i]).value))-kcorr[i][3]

bv=b06-v06-kcorr[:,1]+kcorr[:,2]

SBe=np.zeros_like(b06)

for i in range(0,len(zbest)):
     if zbest[i]<=0.1:
          SBe[i]=b06[i]-kcorr[i][1]+0.726+2.5*np.log10(2*np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
     if zbest[i]<=0.32 and zbest[i]>0.1:
          SBe[i]=v06[i]-kcorr[i][2]+0.726+2.5*np.log10(2*np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
     if zbest[i]<=0.5 and zbest[i]>0.32:
          SBe[i]=r06[i]-kcorr[i][3]+0.726+2.5*np.log10(2*np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])


SBlist=np.where((SBe<21) & (SBe>10))[0]
LCBGS=np.array([])
for i in range(0,len(SBlist)):
     if M[SBlist[i]]<-18.5:
             if bv[SBlist[i]]<0.6:
                     LCBGS=np.append(LCBGS,SBlist[i])

LCBGS=np.ndarray.astype(LCBGS,dtype=int)


