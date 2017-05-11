#Used Lange 2015 to convert F814W re to B re (log(RB)=0.108*log(lambda_814/lambda_B)+log(R814)===>RB=1.0659*R814

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

ID=IDtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
RA=RAtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
DECL=DECtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
SGCLASS=SGCLASStot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
umag=utot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
uerr=uerrtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
bmag=btot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
berr=berrtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
vmag=vtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
verr=verrtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
imag=itot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
ierr=ierrtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
zmag=ztot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
zerr=zerrtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
kmag=ktot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
kerr=kerrtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
zbest=zbesttot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
zuse=zusetot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
zphot=zphottot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
rmag=rtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
rerr=rerrtot[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]
rh=rh[np.where((zbesttot>0)&(zbesttot<0.45)&((zusetot==3)|(zusetot==4))&(SGCLASStot==0)&(np.arccos(np.cos((np.pi/2)-(DECtot*np.pi/180))*np.cos((np.pi/2)-(2.35*np.pi/180))+np.sin((np.pi/2)-(DECtot*np.pi/180))*np.sin((np.pi/2)-(2.35*np.pi/180))*np.cos((RAtot-150.35)*np.pi/180))<0.007653))[0]]

print(rh)

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
kcorrV=-2.5*np.log10(rmarr/rmarr0V)
corrB=-2.5*np.log10(rmarr0B)+0.09
corrV=-2.5*np.log10(rmarr0V)-0.02
corrU=-2.5*np.log10(rmarr0U)-0.79

M=np.zeros_like(zbest)
bv=corrB[:,3]-corrV[:,4]
#M=corrB[:,3]-cosmo.distmod(zbest).value

for i in range(0,len(zbest)):
	if zbest[i]<=0.1:
		M[i]=bmag[i]-0.05122-cosmo.distmod(zbest[i]).value-kcorrM[i][2]
	if zbest[i]<=0.35 and zbest[i]>0.1:
		M[i]=vmag[i]+0.069802-cosmo.distmod(zbest[i]).value-kcorrM[i][3]
	if zbest[i]<=0.55 and zbest[i]>0.35:
		M[i]=rmag[i]-0.01267-cosmo.distmod(zbest[i]).value-kcorrM[i][4]
	if zbest[i]<=0.75 and zbest[i]>0.55:
		M[i]=imag[i]-0.004512-cosmo.distmod(zbest[i]).value-kcorrM[i][5]
	if zbest[i]>0.75:
		M[i]=zmag[i]-0.00177-cosmo.distmod(zbest[i]).value-kcorrM[i][6]

M=M+0.09

#bv=bmag-vmag-kcorrM[:,3]+kcorrV[:,4]

#SBe=np.zeros_like(bmag)
SBe=M+2.5*np.log10((2*np.pi*np.power(cosmo.angular_diameter_distance(zbest).value*np.tan(rh*0.03*4.84814e-6)*1e3*1.0659,2)))+2.5*np.log10((360*60*60/(2*np.pi*0.01))**2)

#for i in range(0,len(zbest)):
#     if zbest[i]<=0.1:
#          SBe[i]=bmag[i]-kcorrM[i][3]+0.726+2.5*np.log10(2*np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
#     if zbest[i]<=0.35 and zbest[i]>0.1:
#          SBe[i]=vmag[i]-kcorrM[i][4]+0.726+2.5*np.log10(2*np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
#     if zbest[i]<=0.55 and zbest[i]>0.35:
#          SBe[i]=rmag[i]-kcorrM[i][5]+0.726+2.5*np.log10(2*np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
#     if zbest[i]<=0.75 and zbest[i]>0.55:
#          SBe[i]=imag[i]-kcorrM[i][6]+0.726+2.5*np.log10(2*np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])
#     if zbest[i]>0.75:
#          SBe[i]=zmag[i]-kcorrM[i][7]+0.726+2.5*np.log10(2*np.pi*np.power(rh[i]*0.03,2))-10*np.log10(1+zbest[i])

print(bv)
LCBGS=np.where((M<=-18)&(SBe<=21.5)&(bv<0.7))[0]

LCBGS=np.ndarray.astype(LCBGS,dtype=int)
#MSTAR=np.power(10,stellarmass)/(2*np.pi*(cosmo.kpc_proper_per_arcmin(zbest).value*rh*0.03*2/60)*cosmo.kpc_proper_per_arcmin(zbest).value*rh*0.03*2/60)
MHIB=2.89-0.34*M
MHID=8.17+1.32*np.log10(cosmo.kpc_proper_per_arcmin(zbest).value*rh*0.03*2*2/60)
detectionlimit1000=2.36*np.power(10,5)*np.power(cosmo.luminosity_distance(zbest).value,2)*5e-5*(31.6/1420405.7*(1+zbest)*300000)*np.sqrt(150/(31.6/1420405.7*(1+zbest)*300000))*3/(1+zbest)
detectionlimit178=2.36*np.power(10,5)*np.power(cosmo.luminosity_distance(zbest).value,2)*7.5e-5*(62.5/1420405.7*(1+zbest)*300000)*np.sqrt(150/(62.5/1420405.7*(1+zbest)*300000))*3/(1+zbest)
detectedB178=np.zeros_like(zbest)
detectedD178=np.zeros_like(zbest)
detectedB1000=np.zeros_like(zbest)
detectedD1000=np.zeros_like(zbest)
p=0
q=0
r=0
s=0
for i in range(0,len(detectedB178)):
	if MHIB[i]>np.log10(detectionlimit178[i]):
		print('detected in 178 hours')
		detectedB178[i]=1
		p=p+1
	if MHIB[i]>np.log10(detectionlimit1000[i]):
		print('detected in 1000 hours')
		detectedB1000[i]=1
		q=q+1
	if MHID[i]>np.log10(detectionlimit178[i]):
		detectedD178[i]=1
		r=r+1
	if MHID[i]>np.log10(detectionlimit1000[i]):
		detectedD1000[i]=1
		s=s+1

print('{} galaxies collected in 178 hours, {} detected in 1000 hours (magnitude)'.format(p,q))
print('{} galaxies collected in 178 hours, {} detected in 1000 hours (diameter)'.format(r,s))
np.set_printoptions(suppress=True)
LCBGFILEARRAY=np.stack((ID[LCBGS],RA[LCBGS],DECL[LCBGS],zbest[LCBGS],zuse[LCBGS],bmag[LCBGS],vmag[LCBGS],rh[LCBGS],corrU[LCBGS,3],corrB[LCBGS,3],corrV[LCBGS,3],MHIB[LCBGS],MHID[LCBGS],detectedB178[LCBGS],detectedD178[LCBGS],detectedB1000[LCBGS],detectedD1000[LCBGS]),axis=-1)
np.savetxt('CHILES_LCBGS_VEGA_HIMASS.txt',LCBGFILEARRAY,header='COSMOSID	RA	DEC	REDSHIFT	REDSHIFT_FLAG	m_{B}	m_{V}	RADIUS	corrU	corrB	corrV	MHIB	MHID	DetectedB178	DetectedD178	DetectedB1000	DetectedD1000')
LCBGTEXARRAY=np.stack((ID[LCBGS],RA[LCBGS],DECL[LCBGS],zbest[LCBGS],rh[LCBGS],corrB[LCBGS,3],corrV[LCBGS,3],MHIB[LCBGS],MHID[LCBGS],detectedB178[LCBGS],detectedD178[LCBGS],detectedB1000[LCBGS],detectedD1000[LCBGS]),axis=-1)
np.savetxt('CHILES_LCBGS_VEGA_HIMASS_TEX.txt',LCBGTEXARRAY,header='COSMOSID	RA	DEC	REDSHIFT	RADIUS	B	V	MHIB	MHID	DetectedB178	DetectedD178	DetectedB1000	DetectedD1000',delimiter='$&$',newline='\\\ \n$',fmt='%5.6f')

