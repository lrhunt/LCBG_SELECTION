#Remove K-Band, filter response may not be accurate, throws off other photometry. LH 2/21/2017
#Lines changed at 62/63,70/75,91,98,103:
#	allmaggies=np.stack((umaggies,kmaggies,bmaggies,vmaggies,rmaggies,imaggies,zmaggies),axis=-1)
#	allinvervar=np.stack((uinvervar,kinvervar,binvervar,vinvervar,rinvervar,iinvervar,zinvervar),axis=-1)
#	carr=np.ndarray((len(bmaggies),6))
#	rmarr=np.ndarray((len(bmaggies),8))
#	rmarr0=np.ndarray((len(bmaggies),8))
#	rmarr0B=np.ndarray((len(bmaggies),8))
#	rmarr0V=np.ndarray((len(bmaggies),8))
#	rmarr0U=np.ndarray((len(bmaggies),8))
#	kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/Lum_Func_Filters_UKS.dat')
#	kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/BESSEL_B.dat')
#	kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/BESSEL_V.dat')
#	kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/BESSEL_U.dat')

import numpy as np
import kcorrect
import astropy as ap
import matplotlib.pyplot as plt
import kcorrect
import kcorrect.utils as ut
import numpy as np
from astropy.cosmology import WMAP9 as cosmo



filein='/home/lrhunt/CATALOGS/PHOT/LAMBDAR_MAG_R.txt'


IDtot,RAtot,DECtot,utot,uerrtot,btot,berrtot,vtot,verrtot,rtot,rerrtot,itot,ierrtot,ztot,zerrtot,ktot,kerrtot,NUVtot,NUVerrtot,rh,zbesttot,zusetot,zphottot,SGCLASStot=np.loadtxt(filein,unpack=True)

ID=IDtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
RA=RAtot[np.where((zbesttot>0)&(zbesttot<1.3)&(zusetot<=2)&(SGCLASStot==0))[0]]
DECL=DECtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
SGCLASS=SGCLASStot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
umag=utot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
uerr=uerrtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
bmag=btot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
berr=berrtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
vmag=vtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
verr=verrtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
imag=itot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
ierr=ierrtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
zmag=ztot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
zerr=zerrtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
kmag=ktot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
kerr=kerrtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
zbest=zbesttot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
zuse=zusetot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
zphot=zphottot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
rmag=rtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
rerr=rerrtot[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]
rh=rh[np.where((zbesttot>0)&(zbesttot<1.2)&(zusetot<=2)&(SGCLASStot==0))[0]]

umaggies=ut.mag2maggies(umag)
kmaggies=ut.mag2maggies(kmag)
bmaggies=ut.mag2maggies(bmag)
vmaggies=ut.mag2maggies(vmag)
rmaggies=ut.mag2maggies(rmag)
imaggies=ut.mag2maggies(imag)
zmaggies=ut.mag2maggies(zmag)

uinvervar=ut.invariance(umaggies,uerr)
kinvervar=ut.invariance(kmaggies,kerr)
binvervar=ut.invariance(bmaggies,berr)
vinvervar=ut.invariance(vmaggies,verr)
rinvervar=ut.invariance(rmaggies,rerr)
iinvervar=ut.invariance(imaggies,ierr)
zinvervar=ut.invariance(zmaggies,zerr)

allmaggies=np.stack((umaggies,bmaggies,vmaggies,rmaggies,imaggies,zmaggies),axis=-1)
allinvervar=np.stack((uinvervar,binvervar,vinvervar,rinvervar,iinvervar,zinvervar),axis=-1)

carr=np.ndarray((len(bmaggies),6))
rmarr=np.ndarray((len(bmaggies),7))
rmarr0=np.ndarray((len(bmaggies),7))
rmarr0B=np.ndarray((len(bmaggies),7))
rmarr0V=np.ndarray((len(bmaggies),7))
rmarr0U=np.ndarray((len(bmaggies),7))
print('Computing k-corrections')
kcorrect.load_templates()
kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/Lum_Func_Filters_US.dat')

for i in range(0,len(carr)):
	carr[i]=kcorrect.fit_nonneg(zbest[i],allmaggies[i],allinvervar[i])
for i in range(0,len(carr)):
	rmarr[i]=kcorrect.reconstruct_maggies(carr[i])
	rmarr0[i]=kcorrect.reconstruct_maggies(carr[i],redshift=0)

kcorrect.load_templates()
kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/BESSEL_B2.dat')

for i in range(0,len(carr)):
	rmarr0B[i]=kcorrect.reconstruct_maggies(carr[i],redshift=0)

kcorrect.load_templates()
kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/BESSEL_V2.dat')

for i in range(0,len(carr)):
	rmarr0V[i]=kcorrect.reconstruct_maggies(carr[i],redshift=0)

kcorrect.load_templates()
kcorrect.load_filters('/home/lrhunt/programs/kcorrect/data/templates/BESSEL_U2.dat')

for i in range(0,len(carr)):
	rmarr0U[i]=kcorrect.reconstruct_maggies(carr[i],redshift=0)


kcorr=-2.5*np.log10(rmarr/rmarr0)
kcorrM=-2.5*np.log10(rmarr/rmarr0B)
corrB=-2.5*np.log10(rmarr0B)
corrV=-2.5*np.log10(rmarr0V)
corrU=-2.5*np.log10(rmarr0U)

M=np.zeros_like(zbest)
bv=corrB[:,3]-corrV[:,4]
for i in range(0,len(zbest)):
	if zbest[i]<=0.1:
		M[i]=bmag[i]-cosmo.distmod(zbest[i]).value-kcorrM[i][3]
	if zbest[i]<=0.35 and zbest[i]>0.1:
		M[i]=vmag[i]-cosmo.distmod(zbest[i]).value-kcorrM[i][4]
	if zbest[i]<=0.55 and zbest[i]>0.35:
		M[i]=rmag[i]-cosmo.distmod(zbest[i]).value-kcorrM[i][5]
	if zbest[i]<=0.75 and zbest[i]>0.55:
		M[i]=imag[i]-cosmo.distmod(zbest[i]).value-kcorrM[i][6]
	if zbest[i]>0.75:
		M[i]=zmag[i]-cosmo.distmod(zbest[i]).value-kcorrM[i][7]

SBe=np.zeros_like(bmag)

for i in range(0,len(zbest)):
     if zbest[i]<=0.1:
          SBe[i]=bmag[i]-kcorrM[i][3]+0.753+2.5*np.log10(np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
     if zbest[i]<=0.35 and zbest[i]>0.1:
          SBe[i]=vmag[i]-kcorrM[i][4]+0.753+2.5*np.log10(np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
     if zbest[i]<=0.55 and zbest[i]>0.35:
          SBe[i]=rmag[i]-kcorrM[i][5]+0.753+2.5*np.log10(np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
     if zbest[i]<=0.75 and zbest[i]>0.55:
          SBe[i]=imag[i]-kcorrM[i][6]+0.753+2.5*np.log10(np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
     if zbest[i]>0.75:
          SBe[i]=zmag[i]-kcorrM[i][7]+0.753+2.5*np.log10(np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])

print(bv)
LCBGS=np.where((M<=-18.5)&(SBe<=21)&(bv<0.6))[0]

LCBGS=np.ndarray.astype(LCBGS,dtype=int)
#MSTAR=np.power(10,stellarmass)/(2*np.pi*(cosmo.kpc_proper_per_arcmin(zbest).value*rh*0.03*2/60)*cosmo.kpc_proper_per_arcmin(zbest).value*rh*0.03*2/60)
MHIB=2.89-0.34*M
MHID=8.17+1.32*np.log10(cosmo.kpc_proper_per_arcmin(zbest).value*rh*0.03*2/60)
#MHINUVR=np.log10(np.power(10,stellarmass)*np.power(10,((-0.322*np.log10(MSTAR))-(0.234*(NUV-r06))+2.817)))
stellarmass=np.zeros_like(MHID)
MHINUVR=np.zeros_like(MHID)
FLUXB=(np.power(10,MHIB)*1+zbest)/(236000*np.power(cosmo.luminosity_distance(zbest).value,2)*(15.625/1420405.7*(1+zbest)*300000)*np.sqrt(200/(15.625/1420405.7*(1+zbest)*300000)))
FLUXD=(np.power(10,MHID)*1+zbest)/(236000*np.power(cosmo.luminosity_distance(zbest).value,2)*(15.625/1420405.7*(1+zbest)*300000)*np.sqrt(200/(15.625/1420405.7*(1+zbest)*300000)))
FLUXNUVR=(np.power(10,MHINUVR)*1+zbest)/(236000*np.power(cosmo.luminosity_distance(zbest).value,2)*(15.625/1420405.7*(1+zbest)*300000)*np.sqrt(200/(15.625/1420405.7*(1+zbest)*300000)))
RMSB=FLUXB/np.power(((7*rh*0.03)/6),2)
RMSD=FLUXD/np.power(((7*rh*0.03)/6),2)
RMSNUVR=FLUXNUVR/np.power(((7*rh*0.03)/6),2)
LCBGFILEARRAY=np.stack((ID[LCBGS],RA[LCBGS],DECL[LCBGS],zbest[LCBGS],zuse[LCBGS],bmag[LCBGS],vmag[LCBGS],rh[LCBGS],stellarmass[LCBGS],MHIB[LCBGS],RMSB[LCBGS],MHID[LCBGS],RMSD[LCBGS],MHINUVR[LCBGS],RMSNUVR[LCBGS]),axis=-1)
np.savetxt('COSMOS_LCBGS.txt',LCBGFILEARRAY,header='COSMOSID	RA	DEC	REDSHIFT	REDSHIFT_FLAG	m_{B}	m_{V}	RADIUS	STELLARMASS	MHIB	RMSB	MHID	RMSD	MHINUVR	RMSNUVR')

