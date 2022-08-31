#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The main functions are find2c, hselect, normalize, rvbina, rvextract, setrvs, spbina, splot, vexplore and vgrid. 
For a more detailed description of how to use this package, please read the ReadMe.
COMPLETAR--https://github.com/NickMilsonPhysics/BinaryStarSolver/blob/master/README.md--
"""
import os
import glob
import os.path 
import numpy as np
import math
from astropy.io import fits
from astropy import units as u
from astropy.io.fits.verify import VerifyWarning
from PyAstronomy import pyasl
from specutils import Spectrum1D
from specutils.manipulation import SplineInterpolatedResampler
from specutils.fitting import fit_continuum
from astropy.modeling.fitting import LinearLSQFitter
from astropy.modeling.polynomial import Chebyshev1D
from scipy import signal
from astropy.constants import c
from scipy.optimize import curve_fit
import matplotlib.pylab as plt
from matplotlib import gridspec
from matplotlib.widgets import Button
from numba import jit
from progress.bar import ChargingBar
from multiprocessing import Pool
import time
from mpl_toolkits.mplot3d import Axes3D
from os import scandir, getcwd
from os.path import abspath
import operator
spline3 = SplineInterpolatedResampler()

print("====================================\nBinary Spectral Disentangling (v0.1)\n====================================\n")
print("Available functions list:\n")
print("\tfind2c\n\thselect\n\tnormalize\n\trvbina\n\trvextract\n\tsetrvs\n\tspbina\n\tsplot\n\tvgrid\n\tvexplore\n\n")
def help(function):
    if function.lower() == 'find2c':
        print('\nFIND2C: detect mass ratio and effective temperature of secondary companion for a spectra dataset.\n')
        print('Mandatory arguments:')
        print('       lis - list of observed spectra to process (use @ for file list)')
        print('       lit - path to folder containing spectra templates')
        print('    vgamma - best known value for systemic radial velocity (km/s)\n')
        print('Optional arguments:')
        print('       spa - string for primary component name')
        print('       spb - string for secondary component name')
        print('      qmin - minimum mass ratio for cross correlation analysis')
        print('      qmax - maximum mass ratio for cross correlation analysis')    
        print('    deltaq - mass increments for cross correlation analysis')
        print('      wreg - spectral regions for cross correlation analysis separated by "-"')
        print('     nproc - number of CPU cores to use in computing process\n')
    elif function.lower() == 'hselect':
        print('\nHSELECT: extract keyword values from file or file list.\n')
        print('Mandatory arguments:')
        print('       img - file or file list to process')
        print('   keyword - keyword to be extracted from each selected file\n')
    elif function.lower() == 'normalize':
        print('\nNORMALIZE: scale each spectrum from dataset to a mean continuum factor.\n')
        print('Mandatory arguments:')
        print('       lis - list of observed spectra to process (use @ for file list)\n')
        print('Optional arguments:')
        print('   interac - show corrected spectra (interactive mode)\n')
    elif function.lower() == 'rvbina':
        print('\nRVBINA: compute radial velocities for spectra dataset\n')
        print('Mandatory arguments:')
        print('       lis - list of observed spectra to process (use @ for file list)\n')
        print('Optional arguments:')
        print('       spa - string for primary component name')
        print('       spb - string for secondary component name')
        print('        ta - spectrum template to use as primary component')
        print('        tb - spectrum template to use as secondary component')
        print('      wreg - spectral regions for cross correlation analysis separated by "-"')
        print('     aconv - absorver convergence factor')
        print('     keyjd - header keyword for Julian Date')
        print('   fitcont - continuum subtract the spectra prior to correlation\n')
        print('   interac - show cross correlation function (interactive mode)\n')
    elif function.lower() == 'rvextract':
        print('\nRVEXTRACT: extract radial velocities from spectra dataset and graphicate RV convergence as function of interations.\n')
        print('Mandatory arguments:')
        print('       lis - list of observed spectra to process (use @ for file list)\n')
        print('Optional arguments:')
        print('    output - output file name')
        print('     graph - show convergence graphics\n')
    elif function.lower() == 'setrvs':
        print('\nSETRVS: set radial velocities for each spectrum from dataset.\n')
        print('Mandatory arguments:')
        print('       lis - list of observed spectra to process (use @ for file list)\n')
        print('Optional arguments:')
        print('        ta - spectrum template to use as primary component')
        print('        tb - spectrum template to use as secondary component')
        print('      wreg - spectral regions for cross correlation analysis separated by "-"')
        print('     keyjd - header keyword for Julian Date')
        print('   fitcont - continuum subtract the spectra prior to correlation\n')
        print('   interac - show cross correlation function (interactive mode)\n')
    elif function.lower() == 'spbina':
        print('\nSPBINA: computing spectra for primary and secondary components.\n')
        print('Mandatory arguments:')
        print('       lis - list of observed spectra to process (use @ for file list)\n')
        print('Optional arguments:')
        print('       spa - string for primary component name')
        print('       spb - string for secondary component name')
        print('       nit - number of iterations')
        print('      frat - flux ratio among components')
        print('    reject - reject pixels using a sigma clipping algorithm')     
        print('         q - mass ratio among components')
        print('    vgamma - best known value for systemic radial velocity (km/s)')      
        print('    obspha - calculate spectra for all phases')
        print('   showtit - enable user interface')
    elif function.lower() == 'splot':
        print('\nSPLOT: plot spectrum from a file in FITS format.\n')
        print('Mandatory arguments:')
        print('      file - file name\n')
        print('Optional arguments:')
        print('      xmin - lower wavelength limit')
        print('      xmax - upper wavelength limit')
        print('      ymin - lower flux limit')
        print('      ymax - upper flux limit')
        print('     scale - flux scale factor')
        print('   markpix - mark pixel values')
        print('    newfig - show spectrum in a new window') 
        print('     color - graph color') 
    elif function.lower() == 'vexplore':
        print('\nVEXPLORE: show results for systemic radial velocities grid analysis\n')
        print('Mandatory arguments:')
        print('       obj - object name to show (extracted from header keyword in VGRID')
    elif function.lower() == 'vgrid':
        print('\nVGRID: detect mass ratio and effective temperature of secondary companion using a systemic radial velocities grid.\n')
        print('Mandatory arguments:')
        print('       lis - list of observed spectra to process (use @ for file list)')
        print('       lit - path to spectra templates folder')
        print('Optional arguments:')
        print('     svmin - lower systemic radial velocity (km/s)')
        print('     svmax - upper systemic radial velocity (km/s)')
        print('      step - radial velocity step for meshing grid (km/s)')
        print('      qmin - minimum mass ratio for cross correlation analysis')
        print('      qmax - maximum mass ratio for cross correlation analysis')    
        print('    deltaq - mass incrememnts forcross correlation analysis')
        print('      wreg - spectral regions for cross correlation analysis separated by "-"')
        print('     nproc - number of CPU cores to use in computing process\n')
    else:
        print('\nUnknown function. Please check the name from functions available list:\n')
        print("\tfft\tfind2c\n\tnormalize\n\trvbina\n\trvextract\n\tsetrvs\n\tspbina\n\tsplot\n\tvgrid\n\tvexplore\n")

################################################################
################################################################
################################################################
def find2c(lis, lit, spa='A', spb='B', vgamma=0, qmin=0.02, qmax=0.5, deltaq=0.01, wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900',nproc=6):
    plt.ion()
    print('\n\tRunning FIND2C (v0.1)\n')
    VerifyWarning('ignore')
    larch=makelist(lis)
    if lit[len(lit)-1] == '/':
        lit=lit[0:len(lit)-1]
    ltemp=sorted(listmp(lit))
    q_array=np.arange(round(qmin,7),round(qmax+deltaq,7),round(deltaq,7))
    path1=os.getcwd()
    #compute B spectrum for each element from q_array
    print('Calculating spectra...')
    pool=Pool(processes=nproc)
    qa2=np.array_split(q_array,nproc)
    pres=[pool.apply_async(qparallel, args= [chk,lis,larch,spa,spb,deltaq,vgamma]) for chk in qa2]
    time.sleep(1)
    pool.close()
    pool.join()
    time.sleep(1)
    pres = [chk.get() for chk in pres]
    pool.terminate()
    print('\t\t\t\t\tDone!')
    print('')
    #create fading mask
    dlist=[]
    hobs = fits.open(larch[0], 'update')
    d1 = hobs[0].header['CDELT1']
    dlist.append(d1)
    try:
        obj1 = hobs[0].header['OBJECT']
        obj1=obj1.replace(' ','')
    except KeyError:
        obj1='NoObject'
    hobs.close(output_verify='ignore')
    htmp = fits.open(ltemp[0], 'update')
    d2 = htmp[0].header['CDELT1']
    dlist.append(d2)
    htmp.close(output_verify='ignore')
    waux1,faux1 = pyasl.read1dFitsSpec(larch[0])
    waux2,faux2 = pyasl.read1dFitsSpec(ltemp[0])
    #gap expresed wavelenght magins in pixels
    gap=50
    winf=max(waux1[0],waux2[0])+gap
    wsup=min(waux1[-1],waux2[-1])-gap
    new_disp_grid,fmask = setregion(wreg,np.max(dlist),winf,wsup)
    #Apply mask over B//q spectra (filtering and apodizing)
    matrix_sq=np.zeros(shape=(len(q_array),len(new_disp_grid)))
    bar1 = ChargingBar('Loading B_q spectra:', max=len(q_array))
    for j,xq in enumerate(q_array):
        aux1=str(round(xq,len(str(deltaq+1))-2)).replace('.','')
        wimg,fimg = pyasl.read1dFitsSpec(spb+aux1+'.fits')
        spec_cont=continuum(wimg, fimg, type='diff',lo=2.5,hi=3.5, graph=False)
        aux_img = Spectrum1D(flux=spec_cont*u.Jy, spectral_axis=wimg*0.1*u.nm)
        aux2_img = spline3(aux_img, new_disp_grid*0.1*u.nm)
        matrix_sq[j] = splineclean(aux2_img.flux.value)*fmask
        bar1.next()
    bar1.finish()
    print('')
    #Cross correlation between each spectrum and B//q
    vector_t=np.zeros(len(ltemp))
    matrix_cc=np.zeros(shape=(len(ltemp),len(q_array)))
    bar2 = ChargingBar('Comparing templates:', max=len(ltemp))
    for k,tmp in enumerate(ltemp):
        htmp = fits.open(tmp, 'update')
        teff= htmp[0].header['TEFF']
        vector_t[k]=teff
        htmp.close(output_verify='ignore')
        wt1,ft1 = pyasl.read1dFitsSpec(tmp)
        temp_cont = continuum(wt1, ft1, type='diff', lo=2.5,hi=6, graph=False)
        aux_tmp1 = Spectrum1D(flux=temp_cont*u.Jy, spectral_axis=wt1*0.1*u.nm)
        aux2_tmp1 = spline3(aux_tmp1, new_disp_grid*0.1*u.nm)
        template1=splineclean(aux2_tmp1.flux.value*fmask)
        tt=np.mean(template1**2)
        for l,xq in enumerate(q_array):
            bb=np.mean(matrix_sq[l]**2)
            tb=np.mean((template1*matrix_sq[l]))
            cc=tb/(np.sqrt(bb)*np.sqrt(tt))
            matrix_cc[k,l]=cc
        bar2.next()
    #save matrix_cc in FITS format
    vgindex=str(vgamma)
    auxout = os.path.isfile(obj1+'_vg_'+vgindex+'.fits')
    if auxout:
        os.remove(obj1+'_vg_'+vgindex+'.fits')
    fits.writeto(obj1+'_vg_'+vgindex+'.fits',matrix_cc,overwrite=True)
    hcc = fits.open(obj1+'_vg_'+vgindex+'.fits', 'update')
    hcc[0].header['Q0'] = qmin
    hcc[0].header['Q1'] = qmax
    hcc[0].header['DELTA_Q'] = deltaq
    hcc[0].header['T0'] = vector_t[0]
    hcc[0].header['T1'] = vector_t[-1]
    hcc[0].header['DELTA_T'] = vector_t[1]-vector_t[0]
    hcc[0].header['VGAMMA'] = vgamma
    hcc.close(output_verify='ignore')
    bar2.finish()
    print('')
    #Estimate mass ratio q and effective temperature Teff for according to a parabole function
    qmed=np.max(matrix_cc,axis=0)
    iq2=int(np.where(qmed == np.max(qmed))[0])
    if qmed[iq2] !=qmed[0] and qmed[iq2] !=qmed[-1]:
        best_q = q_array[iq2]-(qmed[iq2+1]-qmed[iq2-1])*deltaq/2/(qmed[iq2-1]+qmed[iq2+1]-2*qmed[iq2])
    else:
        best_q=q_array[iq2]
    tmed=np.max(matrix_cc,axis=1)
    jt2=np.argmax(tmed)
    #best_B=matrix_sq[iq2]
    rvblist=[]
    for img in larch:
        hdul = fits.open(img, 'update')
        vra = hdul[0].header['VRA']
        vrb = vgamma - (vra-vgamma)/q_array[iq2]
        rvblist.append(round(vrb,3))
        hdul.flush()
        hdul.close()
    print('\t· · · · · · · · · · · · · ·')
    print('\t  Teff='+str(int(vector_t[jt2]))+' K\tq = '+str(round(best_q,2))+'  ')
    print('\t· · · · · · · · · · · · · ·')
    #Graph results for q vs Teff
    fig=plt.figure(figsize=[8,5])
    ax = Axes3D(fig)
    ax.set_xlabel("Mass ratio ($q$)", fontsize=10)
    ax.set_ylabel("Temp. [1000 K]", fontsize=10)
    X, Y = np.meshgrid(q_array, vector_t/1000)
    Z = matrix_cc
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='jet')
    aux3 = os.path.isfile(obj1+'_CC.jpeg')
    if aux3:
        os.remove(obj1+'_CC.jpeg')
    plt.savefig(path1+'/'+obj1+'_CC.jpeg')
    #Graph results for best q and Teff
    fig=plt.figure(figsize=[8,2])
    ax1=plt.subplot(1, 2, 1) 
    ax1.plot(q_array,qmed,marker='',ls='-', color='blue')
    ax1.set_ylabel("Correlation", fontsize=10)
    ax1.set_xlabel("Mass ratio ($q$)", fontsize=10)
    ax1.set(xlim=(min(q_array), max(q_array)))
    plt.axvline(x=best_q, color='black', linestyle='--',linewidth=1)
    ax2=plt.subplot(1, 2, 2) 
    plt.setp(ax2.get_yticklabels(), visible=False)
    ax2.plot(vector_t,tmed,marker='',ls='-', color='red')
    ax2.set_xlabel("Temp. [1000 K]", fontsize=10)
    ax2.set(xlim=(2300, 8400))
    plt.axvline(x=vector_t[jt2], color='black', linestyle='--',linewidth=1)
    plt.subplots_adjust(top=0.995,bottom=0.15,left=0.08,right=0.995,hspace=0.25,wspace=0.025)
    aux4 = os.path.isfile(obj1+'_qT.jpeg')
    if aux4:
        os.remove(obj1+'_qT.jpeg')
    plt.savefig(path1+'/'+obj1+'_qT.jpeg')

################################################################
################################################################
################################################################
def hselect(img,keyword):
    VerifyWarning('ignore')
    larch=makelist(img)
    for i,img in enumerate(larch):
        hdul = fits.open(img)
        try:
            print(img+': '+str(hdul[0].header[keyword]))
        except KeyError:
            print(img)
################################################################
################################################################
################################################################   
def normalize(lis,interac=True):
    #Apply doppler correction for RV=0
    global yours
    plt.close()
    plt.ioff()
    VerifyWarning('ignore')
    larch=makelist(lis)
    #Calculate mean spectrum
    hobs = fits.open(larch[0], 'update')
    dx= hobs[0].header['CDELT1']
    wx1 = hobs[0].header['CRVAL1']
    wx2 = hobs[0].header['CRVAL1'] + dx *(hobs[0].header['NAXIS1']-1)    
    new_disp_grid=np.arange(wx1-50,wx2+50,dx)
    nwvl = len(new_disp_grid)
    obs_matrix = np.zeros(shape=(len(larch),nwvl))
    ldel=[]
    print('\n\tProcessing spectra list...\n')
    larch2=larch.copy()
    for k,img in enumerate(larch):
        wimg,fimg = pyasl.read1dFitsSpec(img)
        rp=np.where(fimg ==0)
        fimg[rp]=np.mean(fimg)
        hdul = fits.open(img, 'update')
        vra = hdul[0].header['VRA']
        hdul.close(output_verify='ignore')
        w2 = wimg * np.sqrt((1.-vra/299792.458)/(1.+vra/299792.458))
        aux_img = Spectrum1D(flux=fimg *u.Jy, spectral_axis=w2 *0.1*u.nm)
        aux_img2 = spline3(aux_img,new_disp_grid*0.1*u.nm)
        if interac:
            fig = plt.figure(figsize=[15,8])
            plt.title('Spectrum '+img)
            plt.xlabel('Wavelenght [A]')
            plt.ylabel('Flux')
            plt.plot(wimg,fimg, color='red',ls='-')  
            plt.plot(new_disp_grid,aux_img2.flux.value,ls='--',marker='',color='blue')
            plt.legend(('Original','RV corrected'))
            plt.tight_layout()
            fig.canvas.mpl_connect('key_press_event', on_key)
            print('\t............................................')
            print('\t:      Press Y to save and continue        :')
            print('\t:  Press any button to discard spectrum    :')
            print('\t············································')
            plt.show()
            if yours == 'y' or yours == 'Y':
                obs_matrix[k]=splineclean(aux_img2.flux.value)
            else:
                ldel.append(k)
                larch2.remove(img)
        else:
            obs_matrix[k]=splineclean(aux_img2.flux.value)
    #delete spectra unselected
    obs_matrix = np.delete(obs_matrix,ldel,axis=0)
    wei=np.mean(obs_matrix,axis=1)/np.max(np.mean(obs_matrix,axis=1))
    smean1 = sfcomb(obs_matrix,wei)
    if interac:      
        fig = plt.figure(figsize=[15,8])
        plt.plot(new_disp_grid,smean1,c='black',ls='--')
        plt.show()
    wtr1 = np.abs(new_disp_grid - (new_disp_grid[0]+50)).argmin(0)
    wtr2 = np.abs(new_disp_grid - (new_disp_grid[-1]-50)).argmin(0)
    grid2=new_disp_grid[wtr1:wtr2]
    for i,img in enumerate(larch2):
        wv,fl = pyasl.read1dFitsSpec(img)
        rp=np.where(fl ==0)
        fl[rp]=np.mean(fl)
        aux_spec = Spectrum1D(flux=fl *u.Jy, spectral_axis=wv *0.1*u.nm)
        spec = spline3(aux_spec, grid2*0.1*u.nm)
        fs=splineclean(spec.flux.value)
        hdul = fits.open(img, 'update')
        vra = hdul[0].header['VRA']
        hdul.close(output_verify='ignore')
        grid_aux = new_disp_grid * np.sqrt((1.+vra/299792.458)/(1.-vra/299792.458))
        aux_smean = Spectrum1D(flux=smean1 *u.Jy, spectral_axis=grid_aux *0.1*u.nm)
        aux_smean2 = spline3(aux_smean,grid2*0.1*u.nm)
        smean2=splineclean(aux_smean2.flux.value)
        s1=fs/smean2
        f_cont = continuum(grid2, s1 ,type='fit',graph=False)
        f_norm = fs / f_cont
        naux = img.replace('.fits','')
        print('')
        print('\tSaving '+img+'...')
        pyasl.write1dFitsSpec(naux+'_NORM.fits', f_norm, grid2, clobber=True)
        copyheader(img,naux+'_NORM.fits')        
################################################################
################################################################
################################################################
def rvbina(lis, spa='A', spb='B', ta='templateA', tb='templateB', 
           wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900', 
           aconv=0.5, keyjd='MJD-OBS', fitcont=True,interac=True):
    plt.ioff()
    spline3 = SplineInterpolatedResampler()
    VerifyWarning('ignore')
    larch=makelist(lis)
    #read spectra
    ta = ta.replace('.fits','')
    tb = tb.replace('.fits','')
    aaux1 = os.path.isfile(ta+'.fits')
    baux1 = os.path.isfile(tb+'.fits')
    if aaux1 == False:
        print('Can not access to primary template spectrum')
        print('END')
    if baux1 == False:
        print('Can not access to secondary template spectrum')
        print('END')
    aaux2 = os.path.isfile(spa+'.fits')
    baux2 = os.path.isfile(spb+'.fits')
    if aaux2 == False:
        print('Can not access to '+spa+'.fits')
        print('END')
    if baux2 == False:
        print('Can not access to '+spb+'.fits')
        print('END')
    #check each spectrum access
    if aaux1 and baux1 and aaux2 and baux2:
        for k,img in enumerate(larch):
            print('\t·············································')
            print('\tProcessing '+img+'...\n')
            wimg,fimg = pyasl.read1dFitsSpec(img)
            hdul = fits.open(img, 'update')
            if k==0:
                #read components mean spectra
                wa,fa = pyasl.read1dFitsSpec(spa+'.fits')
                wb,fb = pyasl.read1dFitsSpec(spb+'.fits')
                #read templates
                wta, fta = pyasl.read1dFitsSpec(ta+'.fits')
                tadul = fits.open(ta+'.fits', mode='readonly')
                del_a = tadul[0].header['CDELT1']
                tadul.close(output_verify='ignore')
                wtb, ftb = pyasl.read1dFitsSpec(tb+'.fits')
                tbdul = fits.open(tb+'.fits', mode='readonly')
                del_b =  tbdul[0].header['CDELT1']
                tbdul.close(output_verify='ignore')
                del_img = hdul[0].header['CDELT1']
                delta=np.min([del_a,del_b,del_img])
                winf=max(wimg[0],wta[0],wtb[0])
                wsup=min(wimg[-1],wta[-1],wtb[-1])
                new_disp_grid,fmask = setregion(wreg,delta,winf,wsup)
                try:
                    obj1 = hdul[0].header['OBJECT']
                except KeyError:
                    obj1='NoObject'
            #read RV for each component
            try:
                vra = hdul[0].header['VRA']
            except KeyError:
                vra = None
            try:
                vrb = hdul[0].header['VRB']
            except KeyError:
                vrb = None
            try:
                xjd = hdul[0].header[keyjd]
            except KeyError:
                xjd = np.nan
#Set RVA if them do not exist
            if vra==None:
                print('\tRV for primary component not found.')
                setrvs(img,ta=ta,wreg=wreg,fintcont=fitcont,keyjd=keyjd) 
            if vrb==None:
                print('\tRV for secondary component not found.')
                setrvs(img,tb=tb,wreg=wreg,fintcont=fitcont,keyjd=keyjd) 
            wlprime_A = wa * np.sqrt((1.+vra/299792.458)/(1.-vra/299792.458))
            wlprime_B = wb * np.sqrt((1.+vrb/299792.458)/(1.-vrb/299792.458))
            aux_img = Spectrum1D(flux=fimg*u.Jy, spectral_axis=wimg*0.1*u.nm)
            aux_sa = Spectrum1D(flux=fa*u.Jy, spectral_axis=wlprime_A*0.1*u.nm)
            aux_sb = Spectrum1D(flux=fb*u.Jy, spectral_axis=wlprime_B*0.1*u.nm)
            #Resampling of each spectrum (lineal interpolation) with the template grid dispersion
            aux2_img = spline3(aux_img, new_disp_grid*0.1*u.nm)
            aux2_sa = spline3(aux_sa, new_disp_grid*0.1*u.nm)
            aux2_sb = spline3(aux_sb, new_disp_grid*0.1*u.nm)
            dsA = aux2_img.flux.value - aux2_sa.flux.value
            dsA = splineclean(dsA)
            dsB = aux2_img.flux.value - aux2_sb.flux.value
            dsB = splineclean(dsB)
            print('Cross-correlation for primary component\n')
            best_vra,aerr,asave=fxcor(new_disp_grid,dsB,wta,fta,fmask,fitcont=fitcont)
            print('Cross-correlation for secondary component\n')
            best_vrb,berr,bsave=fxcor(new_disp_grid,dsA,wtb,ftb,fmask,fitcont=fitcont)
            #Save new RVA values to the FITS files
            if asave:
                new_vra = best_vra * aconv + vra * (1 - aconv)
                err_a = np.abs(new_vra - vra)
                hdul[0].header['VRA'] = round(new_vra,6)
                print(img+', VRA: '+str(round(vra,3))+' km/s --> '+str(round(new_vra,3))+' km/s')
            else:
                new_vra=vra
                err_a=0
            if bsave:
                new_vrb = best_vrb * aconv + vrb * (1 - aconv)
                err_b = np.abs(new_vrb - vrb)
                hdul[0].header['VRB'] = round(new_vrb,6)
                print(img+', VRB: '+str(round(vrb,3))+' km/s --> '+str(round(new_vrb,3))+' km/s')
            else:
                new_vrb=vrb
                err_b=0
            hdul.flush(output_verify='ignore')
            hdul.close(output_verify='ignore')
            name=img.replace('.fits','')
            aux2 = os.path.isfile(name+'.log')
            if aux2 == False:
                flog = open(name+'.log', 'w')
                flog.write('#RV_A\te_it_A\tRV_B\te_it_B\n')
            else:
                flog = open(name+'.log', 'a')
            flog.write(str(round(new_vra,6))+'\t'+str(round(err_a,6))+'\t'+str(round(new_vrb,6))+'\t'+str(round(err_b,6))+'\n')
            flog.close()
            #create vr.txt output file
            aux3 = os.path.isfile(obj1+'_RV.txt')
            if aux3 == False:
                frv = open(obj1+'_RV.txt', 'w')
                frv.write('#JD\tRV_A\te_A\tRV_B\te_B\n')
            else:
                frv = open(obj1+'_RV.txt', 'a')
            frv.write(str(xjd)+'\t'+str(round(new_vra,6))+'\t'+str(round(err_a,6))+'\t'+str(round(new_vrb,6))+'\t'+str(round(err_b,6))+'\n')
            frv.close()
################################################################
################################################################
################################################################
def rvextract(lis, output='rv.txt', graph=True):
    plt.ion()
    plt.close()
    larch=makelist(lis)
    print('\t······························')
    print('\t    Press ENTER to continue') 
    print('\t      or press Q to exit')
    print('\t······························')
    f2 = open(output, 'w')
    f2.write('#RV_A\tRV_B\n')
    for img in larch:
        print('RV convergence for: '+img)
        name=img.replace('.fits','')
        a=np.loadtxt(name+'.log')
        nit=np.arange(1,len(a)+1,1)
        if graph:
            fig = plt.figure(figsize=[8,6])
            ax1 = fig.add_subplot(2, 1, 1)  
            ax1.set_title('Primary RV convergence')
            ax1.set_xlabel('Iteration')
            ax1.set_ylabel('Radial Velocity [km/s]')
            ax1.plot(nit, a[:,0], 'b-')
            ax2 = fig.add_subplot(2, 1, 2, sharex=ax1,sharey=ax1)  
            ax2.set_title('Secondary RV convergence')
            ax2.set_xlabel('Iteration')
            ax2.set_ylabel('Radial Velocity [km/s]')
            ax2.plot(nit, a[:,2], 'r-')
            fig.tight_layout()
            yours = input()
            plt.close()
            if yours == 'q' or yours == 'Q':
                break
        f2.write(str(a[len(a)-1,0])+'\t'+str(a[len(a)-1,2])+'\n')
    f2.close()
################################################################
################################################################
################################################################
def setrvs(lis, ta='templateA', tb=None, wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900', 
           keyjd='MJD-OBS', fitcont=True, interac=True):
    plt.ioff()
    spline3 = SplineInterpolatedResampler()
    VerifyWarning('ignore')
    larch=makelist(lis)
    wmins=[]
    wmaxs=[]
    deltamin=[]
    #if template A exists read 
    if ta !=None:
        ta = ta.replace('.fits','')
        aaux1 = os.path.isfile(ta+'.fits')
        if aaux1:
            wta, fta = pyasl.read1dFitsSpec(ta+'.fits')
            wmins.append(wta[0])
            wmaxs.append(wta[-1])
            tadul = fits.open(ta+'.fits', mode='readonly')
            de1 = tadul[0].header['CDELT1']
            tadul.close(output_verify='ignore')
            deltamin.append(de1)
    #if template B exists read 
    if tb !=None:
        tb = tb.replace('.fits','')
        baux1 = os.path.isfile(tb+'.fits')
        if baux1:
            wtb, ftb = pyasl.read1dFitsSpec(tb+'.fits')
            wmins.append(wtb[0])
            wmaxs.append(wtb[-1])  
            tbdul = fits.open(tb+'.fits', mode='readonly')
            de2 = tbdul[0].header['CDELT1']
            tbdul.close(output_verify='ignore')
            deltamin.append(de2)
    #star analysis
    try:
        for k,img in enumerate(larch):
            print('\t·············································')
            print('\tSet RVs for '+img+'\n')
            wimg,fimg = pyasl.read1dFitsSpec(img)
            hdul = fits.open(img, 'update')
            if k==0:
                de3 = hdul[0].header['CDELT1']
                deltamin.append(de3)
                wmins.append(wimg[0])
                wmaxs.append(wimg[-1])      
                winf=max(wmins)
                wsup=min(wmaxs)
                delta=np.min([deltamin])
                new_disp_grid,fmask = setregion(wreg,delta,winf,wsup)
            aux_img = Spectrum1D(flux=fimg*u.Jy, spectral_axis=wimg*0.1*u.nm)
            aux2_img = spline3(aux_img, new_disp_grid*0.1*u.nm)
            fnew=aux2_img.flux.value
            if ta !=None:
                print('Cross-correlation for primary component\n')
                best_vra,aerr,asave=fxcor(new_disp_grid,fnew,wta,fta,fmask,fitcont=fitcont)
                if asave:
                    hdul[0].header['VRA'] = round(best_vra,6)
                    print(img+', VRA: '+str(round(best_vra,3))+' km/s')
            if tb !=None:
                print('Cross-correlation for secondary component\n')
                best_vrb,berr,bsave=fxcor(new_disp_grid,fnew,wtb,ftb,fmask,fitcont=fitcont)
                if bsave:
                    hdul[0].header['VRB'] = round(best_vrb,6)
                    print(img+', VRB: '+str(round(best_vrb,3))+' km/s')
            hdul.flush(output_verify='ignore')
            hdul.close(output_verify='ignore')
    except:
        print('Can not access to template spectrum')
        print('END')
################################################################
################################################################
################################################################
def spbina(lis, spa='A', spb='B', nit=5, frat=1.0, reject=True,q=None,vgamma=None,obspha=False,showtit=True):
    if showtit:
        print('')
        print('\t  Running SPBINA')
        print('')
    spline3 = SplineInterpolatedResampler()
    VerifyWarning('ignore')
    larch=makelist(lis)
    nimg=len(larch)
    haux = fits.open(larch[0], 'update')
    delta = haux[0].header['CDELT1']
    xwmin=[]
    xwmax=[]
    for i in range(nimg):
        hx = fits.open(larch[i], 'update')
        xdel = hx[0].header['CDELT1']
        xw0 = hx[0].header['CRVAL1']
        xn = hx[0].header['NAXIS1']
        xw1=xw0+xdel*(xn-1)
        xwmin.append(xw0)
        xwmax.append(xw1)
    new_disp_grid=np.arange(np.max(xwmin),np.min(xwmax),delta)
    nwvl = len(new_disp_grid)
    baux1 = os.path.isfile(spb+'.fits')
#STEP 1: create B.fits
    if baux1==False:
        fa = 1.0/(1.0 + frat)
        fb = frat*fa
        tmp_matrix = np.zeros(shape=(nimg,nwvl))
        cont=0
        for img in larch:
            wimg,fimg = pyasl.read1dFitsSpec(img)
            aux_img = Spectrum1D(flux=fimg*u.Jy, spectral_axis=wimg*0.1*u.nm)
            aux2_img = spline3(aux_img, new_disp_grid*0.1*u.nm)
            tmp = aux2_img.flux.value
            tmp_matrix[cont] = splineclean(tmp)
            cont+=1
        f_mean = np.zeros(nwvl)
        for i in range(nwvl):
            f_mean[i] = np.mean(tmp_matrix[:,i])
        f_cont = continuum(new_disp_grid, f_mean,type='fit',graph=False)
        f_cont=splineclean(f_cont)
        B = f_cont* fb
    else:
        waux, faux = pyasl.read1dFitsSpec(spb+'.fits')
        aux_B = Spectrum1D(flux=faux*u.Jy, spectral_axis=waux*0.1*u.nm)
        aux2_B = spline3(aux_B, new_disp_grid*0.1*u.nm)
        tmp2 = aux2_B.flux.value
        B = splineclean(tmp2)
#STEP 2: obs - B.fits
    obs_matrix = np.zeros(shape=(nimg,nwvl))
    dsA_matrix = np.zeros(shape=(nimg,nwvl))
    dsB_matrix = np.zeros(shape=(nimg,nwvl))
    za_matrix = np.zeros(shape=(nimg,nwvl))
    zb_matrix = np.zeros(shape=(nimg,nwvl))
    cont=0
    vra_array = np.zeros(nimg)
    vrb_array = np.zeros(nimg)
    for img in larch:
        hdul = fits.open(img, 'update')
        vra = hdul[0].header['VRA']
        hdul.close(output_verify='ignore')
        if q==None:
            vrb = hdul[0].header['VRB']
        elif q>0:
            if vgamma==None:
                vgamma=0
            vrb = vgamma - (vra-vgamma)/q
        vra_array[cont]=vra
        vrb_array[cont]=vrb
        wimg,fimg = pyasl.read1dFitsSpec(img)
        #doppler correction for B.fits
        wlprime_B = new_disp_grid * np.sqrt((1.+vrb/299792.458)/(1.-vrb/299792.458))
        aux_sb = Spectrum1D(flux=B*u.Jy, spectral_axis=wlprime_B *0.1*u.nm)
        aux2_sb = spline3(aux_sb, new_disp_grid*0.1*u.nm)
        fb_dop = aux2_sb.flux.value
        #Replace np.nan values for the nearest element
        aux_img = Spectrum1D(flux=fimg*u.Jy, spectral_axis=wimg*0.1*u.nm)
        aux2_img = spline3(aux_img, new_disp_grid*0.1*u.nm)
        dsB = aux2_img.flux.value - fb_dop
        dsB_matrix[cont] = splineclean(dsB)
        obs_matrix[cont] = splineclean(aux2_img.flux.value)
        cont+=1
        wei=np.mean(obs_matrix,axis=1)/np.max(np.mean(obs_matrix,axis=1))
#STEP 3: calculate A.fits
    for i in range(nit):
        if showtit==True and i==0:
            bar3 = ChargingBar('Calculating spectra:', max=nit)
        for j in range(nimg):
            wlprime_A = new_disp_grid * np.sqrt((1.-vra_array[j]/299792.458)/(1.+vra_array[j]/299792.458))
            aux_sa = Spectrum1D(flux=dsB_matrix[j] *u.Jy, spectral_axis=wlprime_A *0.1*u.nm)
            aux2_sa = spline3(aux_sa, new_disp_grid*0.1*u.nm)
            fa_dop = aux2_sa.flux.value
            #Replace np.nan values for the nearest element
            za_matrix[j] = splineclean(fa_dop)
        if reject:
            A = sfcomb(za_matrix,wei)
        else:
            A = np.average(za_matrix,axis=0,weights=wei)
#STEP 4: obs - A.fits
#STEP 5: calculate B.fits
        for j in range(nimg):
            wlprime_A = new_disp_grid * np.sqrt((1.+vra_array[j]/299792.458)/(1.-vra_array[j]/299792.458))
            aux_sa = Spectrum1D(flux=A*u.Jy, spectral_axis=wlprime_A *0.1*u.nm)
            aux2_sa = spline3(aux_sa, new_disp_grid*0.1*u.nm)
            fa_dop = aux2_sa.flux.value
            fa_dop =  splineclean(fa_dop)
            dsA_matrix[j] = obs_matrix[j] - fa_dop
            wlprime_B = new_disp_grid * np.sqrt((1.-vrb_array[j]/299792.458)/(1.+vrb_array[j]/299792.458))
            aux_sb = Spectrum1D(flux=dsA_matrix[j]*u.Jy, spectral_axis=wlprime_B*0.1*u.nm)
            aux2_sb = spline3(aux_sb, new_disp_grid*0.1*u.nm)
            fb_dop = aux2_sb.flux.value
            zb_matrix[j] = splineclean(fb_dop)
        if reject:
            B = sfcomb(zb_matrix,wei)
        else:
            B = np.average(zb_matrix,axis=0,weights=wei)
#STEP 6: obs - B.fits
        for j in range(nimg):
            wlprime_B = new_disp_grid * np.sqrt((1.+vrb_array[j]/299792.458)/(1.-vrb_array[j]/299792.458))
            aux_sb = Spectrum1D(flux=B*u.Jy, spectral_axis=wlprime_B *0.1*u.nm)
            aux2_sb = spline3(aux_sb, new_disp_grid*0.1*u.nm)
            fb_dop = aux2_sb.flux.value
            dsB_matrix[j] = obs_matrix[j] - splineclean(fb_dop)
        if showtit:
            bar3.next()
    if showtit:
        bar3.finish()
    pyasl.write1dFitsSpec(spa+'.fits', A, wvl=new_disp_grid, clobber=True)
    pyasl.write1dFitsSpec(spb+'.fits', B, wvl=new_disp_grid, clobber=True)
    if obspha:
        for i,img in enumerate(larch):
            pyasl.write1dFitsSpec('ds-B_'+img, dsB_matrix[i], wvl=new_disp_grid, clobber=True)
            copyheader(img,'ds-B_'+img)
            pyasl.write1dFitsSpec('ds-A_'+img, dsA_matrix[i], wvl=new_disp_grid, clobber=True)
            copyheader(img,'ds-A_'+img)
################################################################
################################################################
################################################################
def splot(file,xmin=None,xmax=None,ymin=None,ymax=None, scale= 1., markpix=False, newfig=True, color='r'):
    plt.ion()
    w,f = pyasl.read1dFitsSpec(file)
    if newfig:
        plt.figure(figsize=[20,10])
    if xmin==None:
        x1=np.min(w)
    else:
        x1=xmin
    if xmax==None:
        x2=np.max(w)
    else:
        x2=xmax
    if ymin==None:
        y1=np.min(f)*scale
    else:
        y1=ymin
    if ymax==None:
        y2=np.max(f)*scale
    else:
        y2=ymax*scale
    plt.axis([x1,x2,y1,y2])  
    plt.ylabel('Flux')
    plt.xlabel('Wavelength')
    plt.title(file)
    plt.plot(w,f*scale,marker='',color=color,linewidth=1)
    if markpix:
        plt.plot(w,f*scale,marker='.',markersize=2,color='black',linestyle='')
    plt.tight_layout()
################################################################
################################################################
################################################################
def vgrid(lis, lit, svmin=-1, svmax=1, step=0.1, qmin=0.02, qmax=0.5, deltaq=0.01, wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900',nproc=6):
    istr='00'
    if os.path.isdir('output_00')==False:
        os.mkdir('output_00')
    else:
        it=1
        while os.path.isdir('output_'+istr):
            if it <10:
                istr='0'+str(it)
            else:
                istr=str(it)
            it+=1
        os.mkdir('output_'+istr)
    outfolder='output_'+istr
    print('Output folder: '+outfolder)
    svrange=np.arange(svmin,svmax+step,step) 
    if math.modf(step)[0] == 0:
        nrd=0
    else:
        nrd=len(str(step))-str(step).find('.')-1
    plt.ion()
    print('\n\t  Running VGRID (v0.1)\n')
    VerifyWarning('ignore')
    larch=makelist(lis)
    if lit[len(lit)-1] == '/':
        lit=lit[0:len(lit)-1]
    ltemp=sorted(listmp(lit))
    q_array=np.arange(round(qmin,7),round(qmax+deltaq,7),round(deltaq,7))
    #Create fading mask
    dlist=[]
    hobs = fits.open(larch[0], 'update')
    d1 = hobs[0].header['CDELT1']
    dlist.append(d1)
    try:
        obj1 = hobs[0].header['OBJECT']
        obj1=obj1.replace(' ','')
    except KeyError:
        obj1='NoObject'
    hobs.close(output_verify='ignore')
    htmp = fits.open(ltemp[0], 'update')
    d2 = htmp[0].header['CDELT1']
    dlist.append(d2)
    htmp.close(output_verify='ignore')
    waux1,faux1 = pyasl.read1dFitsSpec(larch[0])
    waux2,faux2 = pyasl.read1dFitsSpec(ltemp[0])
    gap=50
    winf=max(waux1[0],waux2[0])+gap
    wsup=min(waux1[-1],waux2[-1])-gap
    new_disp_grid,fmask = setregion(wreg,np.max(dlist),winf,wsup)
    #Load templates and create array for temperatures
    vector_t=np.zeros(len(ltemp))
    tt_array=np.zeros(len(ltemp))
    matrix_tmp=np.zeros(shape=(len(ltemp),len(new_disp_grid)))
    bar2 = ChargingBar('Loading templates:', max=len(ltemp))
    for k,tmp in enumerate(ltemp):
        htmp = fits.open(tmp, 'update')
        teff= htmp[0].header['TEFF']
        vector_t[k]=teff
        htmp.close(output_verify='ignore')
        wt1,ft1 = pyasl.read1dFitsSpec(tmp)
        temp_cont = continuum(wt1, ft1, type='diff', lo=2.5,hi=6, graph=False)
        aux_tmp1 = Spectrum1D(flux=temp_cont*u.Jy, spectral_axis=wt1*0.1*u.nm)
        aux2_tmp1 = spline3(aux_tmp1, new_disp_grid*0.1*u.nm)
        template1=splineclean(aux2_tmp1.flux.value*fmask)
        matrix_tmp[k]=template1
        tt_array[k]=np.mean(template1**2)
        bar2.next()
    bar2.finish()
    print('')
    #Execute FIND2C analysis for vgamma grid
    bar0 = ChargingBar('Calculating syst. rv:', max=len(svrange))
    for vgamma in svrange:
        matrix_cc=np.zeros(shape=(len(ltemp),len(q_array)))
        vgindex=str(round(vgamma,nrd))
        for xq in q_array:
            aux1=str(round(xq,len(str(deltaq+1))-2)).replace('.','')
            aux2 = os.path.isfile('B'+aux1+'.fits')
            if aux2:
                os.remove('B'+aux1+'.fits')
            aux3 = os.path.isfile('A'+aux1+'.fits')
            if aux3:
                os.remove('A'+aux1+'.fits') 
        pool=Pool(processes=nproc)
        qa2=np.array_split(q_array,nproc)
        pres=[pool.apply_async(qparallel, args= [chk,lis,larch,'A','B',deltaq,vgamma]) for chk in qa2]
        time.sleep(1)
        pool.close()
        pool.join()
        time.sleep(1)
        pres = [chk.get() for chk in pres]
        pool.terminate()
        matrix_sq=np.zeros(shape=(len(q_array),len(new_disp_grid)))
        #Load calculated B_q spectra
        for j,xq in enumerate(q_array):
            aux1=str(round(xq,len(str(deltaq+1))-2)).replace('.','')
            wimg,fimg = pyasl.read1dFitsSpec('B'+aux1+'.fits')
            spec_cont=continuum(wimg, fimg, type='diff',lo=2.5,hi=3.5, graph=False)
            aux_img = Spectrum1D(flux=spec_cont*u.Jy, spectral_axis=wimg*0.1*u.nm)
            aux2_img = spline3(aux_img, new_disp_grid*0.1*u.nm)
            matrix_sq[j] = splineclean(aux2_img.flux.value)*fmask
        #Load calculated tt values from templates
        for k,tt in enumerate(tt_array):
            for l,xq in enumerate(q_array):
                bb=np.mean(matrix_sq[l]**2)
                tb=np.mean(matrix_tmp[k]*matrix_sq[l])
                cc=tb/(np.sqrt(bb)*np.sqrt(tt))
                matrix_cc[k,l]=cc
        bar0.next()
        #save matrix_cc in FITS file
        fits.writeto(outfolder+'/'+obj1+'_vg_'+vgindex+'.fits',matrix_cc,overwrite=True)
        hcc = fits.open(outfolder+'/'+obj1+'_vg_'+vgindex+'.fits', 'update')
        hcc[0].header['Q0'] = qmin
        hcc[0].header['Q1'] = qmax
        hcc[0].header['DELTA_Q'] = deltaq
        hcc[0].header['T0'] = vector_t[0]
        hcc[0].header['T1'] = vector_t[-1]
        hcc[0].header['DELTA_T'] = vector_t[1]-vector_t[0]
        hcc[0].header['VGAMMA'] = vgamma
        hcc.close(output_verify='ignore') 
        #Clean A.fits and B.fits files
        for xq in q_array:
            aux1=str(round(xq,len(str(deltaq+1))-2)).replace('.','')
            aux2 = os.path.isfile('B'+aux1+'.fits')
            if aux2:
                os.remove('B'+aux1+'.fits')
            aux3 = os.path.isfile('A'+aux1+'.fits')
            if aux3:
                os.remove('A'+aux1+'.fits')
    vexplore(outfolder)
    bar0.finish()
################################################################
################################################################
################################################################
def vexplore(folder):
    plt.ioff()
    folder=folder.replace('/','')
    #list vgamma files into output
    nfiles = os.listdir(folder)
    vgdicc = {}
    for xf in nfiles:
        if os.path.isfile(os.path.join(folder, xf)) and xf.endswith('.fits'):
            hdul = fits.open(folder+'/'+xf, 'update')
            xvg=hdul[0].header['VGAMMA']
            hdul.close(output_verify='ignore')
            vgdicc[xf]=xvg
    #sort list
    vgd_sort=sorted(vgdicc.items(), key=operator.itemgetter(1), reverse=False)
    n=len(vgd_sort)
    #read pixels values and q_array and Teff array
    qall=[]
    tall=[]
    ccall=[]
    for i,img in enumerate(vgd_sort):
        hdul = fits.open(folder+'/'+img[0], 'update')
        #make array for q
        qmin=hdul[0].header['Q0']
        qmax=hdul[0].header['Q1']
        deltaq=hdul[0].header['DELTA_Q']
        q_array=np.arange(round(qmin,7),round(qmax+deltaq,7),round(deltaq,7))
        #make array for Teff
        t0=hdul[0].header['T0']  
        t1=hdul[0].header['T1']
        deltat=hdul[0].header['DELTA_T']
        vector_t=np.arange(t0,t1+deltat,deltat)
        #read pixels values
        ccx = hdul[0].data
        hdul.close(output_verify='ignore')
        qall.append(q_array)
        tall.append(vector_t)
        ccall.append(ccx)
    fig=plt.figure(figsize=[8,6])
    ax = Axes3D(fig)
    ax.set_autoscalez_on(False)
    ax.set_zlim3d(bottom=np.min(ccall), top=np.max(ccall))
    ax.set_xlabel("Mass ratio ($q$)", fontsize=10)
    ax.set_ylabel("Temp. [1000 K]", fontsize=10)
    ax.set_title('vgamma = '+str(vgd_sort[0][1])+' km/s',y=1)
    X, Y = np.meshgrid(qall[0], tall[0]/1000)
    ax.plot_surface(X, Y, ccall[0], rstride=1, cstride=1, cmap='jet',vmin=np.min(ccall), vmax=np.max(ccall))
    class Index:
        ind = 0
        def next(self, event):
            self.ind += 1
            i=self.ind  % n
            ax.cla()
            ax.set_zlim3d(bottom=np.min(ccall), top=np.max(ccall))
            ax.plot_surface(X, Y, ccall[i], rstride=1, cstride=1, cmap='jet',vmin=np.min(ccall), vmax=np.max(ccall))
            ax.set_title('vgamma = '+str(vgd_sort[i][1])+' km/s',y=1)
            ax.set_xlabel("Mass ratio ($q$)", fontsize=10)
            ax.set_ylabel("Temp. [1000 K]", fontsize=10)
            plt.draw()
    
        def prev(self, event):
            self.ind -= 1 % n
            i=self.ind % n
            ax.cla()
            ax.set_zlim3d(bottom=np.min(ccall), top=np.max(ccall))
            ax.plot_surface(X, Y, ccall[i], rstride=1, cstride=1, cmap='jet',vmin=np.min(ccall), vmax=np.max(ccall))
            ax.set_title('vgamma = '+str(vgd_sort[i][1])+' km/s',y=1)
            ax.set_xlabel("Mass ratio ($q$)", fontsize=10)
            ax.set_ylabel("Temperature [x1000 K]", fontsize=10)
            plt.draw()
    callback = Index()
    axprev = plt.axes([0.05, 0.88, 0.1, 0.075])
    axnext = plt.axes([0.85, 0.88, 0.1, 0.075])
    bnext = Button(axnext, 'Next')
    bnext.on_clicked(callback.next)
    bprev = Button(axprev, 'Previous')
    bprev.on_clicked(callback.prev)
    plt.show()
################################################################
################################################################
################################################################
# INTERNAL FUNCTIONS
#Parallel process
def qparallel(q_array,lis,larch,spa,spb,deltaq,vgamma):
    for xq in q_array:
        aux1=str(round(xq,len(str(deltaq+1))-2)).replace('.','')
        aux2 = os.path.isfile(spb+aux1+'.fits')
        if aux2==False:
            spbina(lis, spa=spa+aux1, spb=spb+aux1, nit=1,q=xq,vgamma=vgamma,showtit=False)

def makelist(lis):
    aux1=lis.find('@')
    if aux1 == -1:
        path=os.getcwd()
        aux2=glob.glob(path+'/'+lis)
        newlist=[]
        for i in range(len(aux2)):
            newlist.append(aux2[i].replace(path+'/',''))
    else:
        aux2=lis.replace('@','')
        f=open(aux2,'r')
        aux3=f.readlines()
        f.close()
        newlist=[]
        for i in range(len(aux3)):
            newlist.append(aux3[i].rstrip('\n'))
    return(newlist)

@jit(nopython=True)
def sfcomb(aflux,wei,nit=5,sigma=3):
    nwvl=aflux.shape[1]
    s_mean = []
    for w in range(nwvl):
        fdata = aflux[:,w]/wei
        weights=wei
        for i in range(nit):
            fmean = np.sum(fdata*weights)/np.sum(weights)
            fsig=np.sqrt(np.sum((fdata-fmean)**2*weights)/np.sum(weights))
            faux=fdata
            cont=0
            for j,fj in enumerate(faux):
                if fj < (fmean-sigma*fsig) or fj > (fmean+sigma*fsig):
                    fdata=np.delete(fdata, j-cont)
                    weights=np.delete(weights,j-cont)
                    cont+=1
        fmean = np.sum(fdata*weights)/np.sum(weights)
        fsig=np.sqrt(np.sum((fdata-fmean)**2)/np.sum(weights))
        s_mean.append(fmean)
    return(np.array(s_mean))
            
def setregion(wreg,delta,winf,wsup,amort=0.1):
    reg1=wreg.split(',')
    reg2=[]
    for i,str1 in enumerate(reg1):
        reg2.append([int(str1.split('-')[0]),int(str1.split('-')[1])])
    reg3=[]
    stat1=True
    for j,wvx in enumerate(reg2):
        x1=wvx[0]
        x2=wvx[1]
        if stat1:
            if x1>=winf:
                reg3.append(wvx)
                stat1=False
            elif x1<winf and x2>winf:
                wvx[0]=winf
                reg3.append(wvx)
                stat1=False
            elif x1<winf and x2<=winf:
                stat1=True
        else:
            if x1>wsup and x2>wsup:
                break
            elif x1<wsup and x2>=wsup:
                wvx[1]=wsup
                reg3.append(wvx)
            elif x1<wsup and x2<wsup:
                reg3.append(wvx)
    wvl=np.arange(reg3[0][0],reg3[-1][1]+delta,delta)
    f=np.zeros(len(wvl))
    for k,interv in enumerate(reg3):
        x1=interv[0]
        x2=interv[1]
        i1 = np.abs(wvl - x1).argmin(0)
        i2 = np.abs(wvl - x2).argmin(0)
        am2=amort*(x2-x1)
        xarr=wvl[i1:i2]
        mask=np.zeros(len(xarr))
        for k,w in enumerate(xarr):
            if w<=(x1+am2):
                mask[k]=np.sin(np.pi*(w-x1)/(2*am2))
            elif w>(x1+am2) and w<(x2-am2):
                mask[k]=1
            else:
                mask[k]=np.cos(np.pi*(w-x2+am2)/(2*am2))
        f[i1:i2]=mask
    return(wvl,f)

def continuum(w,f, order=12, type='fit', lo=2, hi=3, nit=20, graph=True):
    w_cont=w.copy()
    f_cont=f.copy()
    sigma0=np.std(f_cont)
    wei=~np.isnan(f_cont)*1
    i=1
    nrej1=0
    while i < nit:
        c0=np.polynomial.chebyshev.Chebyshev.fit(w_cont,f_cont,order,w=wei)(w_cont)
        resid=f_cont-c0
        sigma0=np.sqrt(np.average((resid)**2, weights=wei))
        wei = 1*np.logical_and(resid>-lo*sigma0,resid<sigma0*hi)
        nrej=len(wei)-np.sum(wei)
        if nrej==nrej1:
            break
        nrej1=nrej
        i=i+1
    s1=Spectrum1D(flux=c0*u.Jy, spectral_axis=w_cont*0.1*u.nm) 
    c1= fit_continuum(s1, model=Chebyshev1D(order),fitter=LinearLSQFitter())
    if type=='fit':
        fout=c1(w*0.1*u.nm).value
    elif type=='ratio':
        fout=f_cont/c1(w*0.1*u.nm).value
    elif type=='diff':
        fout=f_cont-c1(w*0.1*u.nm).value
    if graph:
        fig = plt.figure(figsize=[20,10])
        ngrid = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[5, 1])
        ax1 = fig.add_subplot(ngrid[0])
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.ylabel('Flux')
        ax1.plot(w,f,color='gray')
        ax1.plot(w_cont,f_cont,color='blue',linestyle='',marker='.',markersize=2)
        ax1.plot(w_cont,c1(w_cont*0.1*u.nm).value,c='r',linestyle='--')
        ax2 = fig.add_subplot(ngrid[1], sharex=ax1,sharey=ax1)
        plt.xlabel('Wavelength [nm]')
        ax2.plot(w,(f-c1(w*0.1*u.nm).value),color='gray',linestyle='',marker='.',markersize=1)
        ax2.plot(w_cont,(f_cont-c1(w_cont*0.1*u.nm).value),color='blue',linestyle='',marker='.',markersize=1)
        ax2.axhline(y=0, color='red', linestyle='--',linewidth=1)
        plt.tight_layout()
        plt.show()
    return(fout)
    
def on_key(event):
    global yours,nmin,rv1,rv2,gbase,xcent
    #print('You pressed: ', event.key)
    yours=event.key
    #print('x = ',event.xdata)
    #print('y = ',event.ydata)
    if event.key=='y' or event.key=='Y':
        plt.close()
    elif event.key=='m' or event.key=='M':
        if nmin:
            rv1=event.xdata
            print('RV min: '+str(round(rv1,3))+' km/s')
            nmin=False
        else:
            rv2=event.xdata
            print('RV max: '+str(round(rv2,3))+' km/s')
            nmin=True
            yours='r'
            plt.close()
    elif event.key=='b' or event.key=='B':
        gbase=event.ydata
        plt.close()
    elif event.key=='c' or event.key=='C':
        xcent=event.xdata
        rv1=xcent-10
        rv2=xcent+10
        plt.close()
    event.canvas.draw()

def Gauss(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))
    
def fxcor(w, f, wt, ft, mask, fitcont=True, rvcent=None):
    global yours,nmin,rv1,rv2,gbase,xcent
    nmin=True
    plt.ioff()
    #adjust continuum
    if fitcont:
        fci=continuum(w=w, f=f, type='diff', graph=False)
        fct=continuum(w=wt, f=ft, type='diff', graph=False)
    else:
        fci=f
        fct=f
    #convert grid to ln w
    aux_grid=np.log(w)
    dlog = aux_grid[-1] - aux_grid[-2]
    new_log_grid=np.arange(aux_grid[0],aux_grid[-1],dlog)
    #rescale mask
    aux_mask = Spectrum1D(flux=mask*u.Jy, spectral_axis=np.log(w)*0.1*u.nm)
    aux2_mask= spline3(aux_mask, new_log_grid*0.1*u.nm)
    log_mask = splineclean(aux2_mask.flux.value)
    #rescale image spectrum
    aux_img = Spectrum1D(flux=fci*u.Jy, spectral_axis=np.log(w)*0.1*u.nm)
    aux2_img = spline3(aux_img, new_log_grid*0.1*u.nm)
    log_img=splineclean(aux2_img.flux.value)
    fi2=log_img*log_mask
    #rescale template spectrum
    aux_sa = Spectrum1D(flux=fct*u.Jy, spectral_axis=np.log(wt)*0.1*u.nm)
    aux2_sa = spline3(aux_sa, new_log_grid*0.1*u.nm)
    log_tmp = splineclean(aux2_sa.flux.value)
    ft2=log_tmp*log_mask
    #execute correlate
    cc1=signal.correlate(fi2,ft2,method='fft')
    #normalize cc1
    cc1=cc1/(np.sqrt(np.sum(np.power(fi2,2)))*(np.sqrt(np.sum(np.power(ft2,2)))))
    #find peak value
    i1=int(np.where(cc1==max(cc1))[0])
    #x-lags for cc1
    lamlog1=new_log_grid-new_log_grid[0]
    lamlog2=-new_log_grid+new_log_grid[0]
    lamlog2.sort()
    llog=np.concatenate((lamlog2[0:-1],lamlog1),axis=0)
    axisrv=llog*c.value/1000
    #adjust 
    if rvcent==None:
        ccmax=np.argmax(cc1)
        xcent=axisrv[ccmax]
    else:
        xcent=rvcent
    stat=True
    rv1=xcent-10
    rv2=xcent+10
    gbase=np.mean(cc1)
    while stat:
        near0=axisrv.flat[np.abs(axisrv - xcent).argmin()]
        i0=int(np.where(axisrv == near0)[0])
        near1=axisrv.flat[np.abs(axisrv - min(rv1,rv2)).argmin()]
        i1=int(np.where(axisrv == near1)[0])
        near2=axisrv.flat[np.abs(axisrv - max(rv1,rv2)).argmin()]
        i2=int(np.where(axisrv == near2)[0])
        try:
            xb=axisrv[i1:i2]
            yb=cc1[i1:i2]
            mb = np.sum(xb * (yb-gbase)) / np.sum(yb-gbase)
            sigb = np.sqrt(np.abs(np.sum((yb-gbase) * (xb - mb)**2) / np.sum(yb-gbase)))
            pb1,pb2 = curve_fit(Gauss, xb, yb-gbase, p0=[np.max(yb-gbase), mb, sigb])
            ybgauss=Gauss(xb, *pb1)+gbase
            xrv=pb1[1]
            rverr = np.sqrt(np.diag(pb2))[1]
            sfit=True
            serr=str(round(rverr,3))
        except Exception:
            xrv=axisrv[np.argmax(cc1)]
            rverr=np.nan
            sfit=False
            serr=str(rverr)
        print('Radial Velocity: '+str(round(xrv,3))+' +/- '+serr+' km/s')
        #graficar
        fig, (a0, a1) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[1, 5]})
        plt.setp(a0.get_yticklabels(), visible=False)
        plt.setp(a0.get_xticklabels(), visible=False)
        a0.yaxis.offsetText.set_visible(False)
        a0.plot(axisrv,cc1,color='blue')
        a0.axvline(axisrv[i0], color='red', linestyle='--',linewidth=1)
        a1.set_xlabel("Radial Velocity [km/s]", fontsize=10)
        a1.set_ylabel("Correlation", fontsize=10)
        a1.plot(axisrv,cc1,color='blue')
        a1.set(xlim=(axisrv[i1-10], axisrv[i2+10]))
        ymin2=min(cc1[i1:i2])
        ymax2=max(cc1[i1:i2])
        ex2=(ymax2-ymin2)*0.1
        a1.set(ylim=(ymin2-ex2, ymax2+ex2))
        if sfit:
            plt.plot(xb, yb, color='black',marker='.',linestyle='')
            plt.plot(xb, ybgauss, color='green', label='fit',linestyle='--')
        else:
            plt.legend(('No fit'))
        plt.tight_layout()
        fig.canvas.mpl_connect('key_press_event', on_key)
        print('\t..........................................')
        print('\t:   Press M to mark gaussian regions     :')
        print('\t:   Press B to set background level      :')
        print('\t:     Press Y to save and continue       :')
        print('\t:     Press Q to discard spectrum        :') 
        print('\t··········································')
        plt.show()
        if yours == 'y' or yours == 'Y':
            fsave = True
            print('\t\t\t\t\tDone!')
            plt.close()
            break
        elif yours == 'q' or yours == 'Q':
            fsave = False
            plt.close()
            break
    return xrv, serr, fsave

@jit(nopython=True)
def splineclean(fspl):
    if np.isnan(fspl[0]):
        cinit=1
        while np.isnan(fspl[cinit]):
            cinit+=1
        finit=fspl[cinit]
        fspl[0:cinit+1]=finit
    if np.isnan(fspl[-1]):
        cend=-2
        while np.isnan(fspl[cend]):
            cend-=1
        fend=fspl[cend]
        fspl[cend:len(fspl)]=fend
    return(fspl)

def copyheader(img1,imgout):
    hdul = fits.open(img1, 'update')
    hnorm =  fits.open(imgout, 'update')
    listk=('CDELT1','CTYPE1','BUNIT','ORIGIN','DATE','TELESCOP','INSTRUME',
           'OBJECT','RA','DEC','EQUINOX','RADECSYS','EXPTIME','MJD-OBS','DATE-OBS','UTC','LST',
           'PI-COI','CTYPE1','CTYPE2','ORIGFILE','UT','ST','AIRMASS','VRA','VRB')
    for i,k in enumerate(listk):
        try:
            hnorm[0].header[k] = hdul[0].header[k]
        except KeyError:
            print('Keyword '+k+' not found')
    hnorm.flush(output_verify='ignore')
    hnorm.close(output_verify='ignore')
    hdul.close(output_verify='ignore')

def listmp(ruta = getcwd()):
    return [abspath(arch.path) for arch in scandir(ruta) if arch.is_file()]
