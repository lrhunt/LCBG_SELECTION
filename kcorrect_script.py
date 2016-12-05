# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import WMAP9 as cosmo
import math as m
zarrrm,b06rm,v06rm,r06rm,f81406rm,z06rm=np.loadtxt('outmaggies.txt',unpack=True)
zarr0,b06rm0,v06rm0,r06rm0,f81406rm0,z06rm0=np.loadtxt('outmaggies_0.txt',unpack=True)
b06,berr,v06,verr,r06,rerr,i06,ierr,z06,zerr,zbest,classif,zuse=np.loadtxt('zCOSMOS_REGION_ALL.txt',skiprows=1,unpack=True)
ranged=np.where(np.logical_and(zarrrm>0,zarrrm<1.3))[0]
kb=-2.5*np.log10(b06rm[ranged]/b06rm0[ranged])
kv=-2.5*np.log10(v06rm[ranged]/v06rm0[ranged])
kr=-2.5*np.log10(r06rm[ranged]/r06rm0[ranged])
kf814=-2.5*np.log10(f81406rm[ranged]/f81406rm0[ranged])
kz=-2.5*np.log10(z06rm[ranged]/z06rm0[ranged])
b06r=b06[ranged]
v06r=v06[ranged]
r06r=r06[ranged]
i06r=i06[ranged]
z06r=z06[ranged]
zbestr=zbest[ranged]
selection=np.random.random_integers(0,len(ranged)-1,100)
Mbcorr=b06r[selection]-5*(np.log10(100000*cosmo.luminosity_distance(zbestr[selection]).value))-(0.1055*np.power(zbestr[selection],6)-0.6393*np.power(zbestr[selection],5)+0.977*np.power(zbestr[selection],4)+1.1534*np.power(zbestr[selection],3)-4.6431*np.power(zbestr[selection],2)+3.8905*zbestr[selection]-0.0483)
Mb=b06r[selection]-5*(np.log10(100000*cosmo.luminosity_distance(zbestr[selection]).value))-kb[selection]
Mv=v06r[selection]-5*(np.log10(100000*cosmo.luminosity_distance(zbestr[selection]).value))-kv[selection]
Mr=r06r[selection]-5*(np.log10(100000*cosmo.luminosity_distance(zbestr[selection]).value))-kr[selection]
Mi=i06r[selection]-5*(np.log10(100000*cosmo.luminosity_distance(zbestr[selection]).value))-kf814[selection]
Mz=z06r[selection]-5*(np.log10(100000*cosmo.luminosity_distance(zbestr[selection]).value))-kz[selection]
if len(np.where(Mb<-30)[0])>0:
	Mb[np.where(Mb<-30)[0]]=np.nan

if len(np.where(Mv<-30)[0])>0:
	Mv[np.where(Mv<-30)[0]]=np.nan

if len(np.where(Mr<-30)[0])>0:
	Mr[np.where(Mr<-30)[0]]=np.nan

if len(np.where(Mi<-30)[0])>0:
	Mi[np.where(Mi<-30)[0]]=np.nan

if len(np.where(Mz<-30)[0])>0:
	Mz[np.where(Mz<-30)
plt.scatter(zbestr[selection],Mb-Mv)
plt.ylabel('Difference in  M$_{B}$ (calculated using m$_{b}$, k$_{b}$ and m$_{v}$, k$_{v}$')
plt.xlabel('Redshift (z)')
plt.savefig('Mb-Mv.png')
plt.scatter(zbestr[selection],Mb-Mr)
plt.ylabel('Difference in  M$_{B}$ (calculated using m$_{B}$, k$_{B}$ and m$_{r}$, k$_{r}$)')
plt.xlabel('Redshift (z)')
plt.savefig('Mb-Mr.png')
plt.scatter(zbestr[selection],Mb-Mi)
plt.ylabel('Difference in  M$_{B}$ (calculated using m$_{B}$, k$_{B}$ and m$_{F814w}$, k$_{F814w}$)')
plt.xlabel('Redshift (z)')
plt.savefig('Mb-Mi.png')
plt.scatter(zbestr[selection],Mb-Mz)
plt.ylabel('Difference in  M$_{B}$ (calculated using m$_{B}$, k$_{B}$ and m$_{z}$, k$_{z}$)')
plt.xlabel('Redshift (z)')
plt.savefig('Mb-Mz.png')
plt.scatter(zbestr[selection],Mv-Mr)
plt.ylabel('Difference in  M$_{B}$ (calculated using m$_{V}$, k$_{V}$ and m$_{r}$, k$_{r}$)')
plt.xlabel('Redshift (z)')
plt.savefig('Mv-Mr.png')
plt.scatter(zbestr[selection],Mv-Mi)
plt.ylabel('Difference in  M$_{B}$ (calculated using m$_{V}$, k$_{V}$ and m$_{F814W}$, k$_{F814W}$)')
plt.xlabel('Redshift (z)')
plt.savefig('Mv-Mi.png')
plt.scatter(zbestr[selection],Mv-Mz)
plt.ylabel('Difference in  M$_{B}$ (calculated using m$_{V}$, k$_{V}$ and m$_{z}$, k$_{z}$)')
plt.xlabel('Redshift (z)')
plt.savefig('Mv-Mz.png')
plt.scatter(zbestr[selection],Mb-Mbcorr)
plt.ylabel('Difference in  M$_{B}$ (using old k-correction and new k-correction')
plt.xlabel('Redshift (z)')
plt.savefig('difoldnewk.png')
Poggianti=plt.scatter(zbestr[selection],Mbcorr,color='blue',label='M$_{B}$ using k-correction from Poggianti, 1997')
Blanton=plt.scatter(zbestr[selection],Mb,color='red',label='M$_{B}$ using k-correction from Blanton')
plt.legend(handles=[Poggianti,Blanton])
plt.ylabel('M$_{B}$')
plt.xlabel('Redshift (z)')
plt.savefig('oldkvsnewk.png')


