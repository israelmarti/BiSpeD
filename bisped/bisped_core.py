#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Módulo interno de bisped - Contiene todas las funciones auxiliares y de procesamiento (version 1.7)
No debe ser importado directamente por el usuario
"""

import os
import glob
import sys
import time
import math
import operator
import gc
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.io.fits.verify import VerifyWarning
from PyAstronomy import pyasl
from specutils import Spectrum
from specutils.manipulation import SplineInterpolatedResampler
from specutils.fitting import fit_continuum
from astropy.modeling.fitting import LinearLSQFitter
from astropy.modeling.polynomial import Chebyshev1D
from scipy import signal
from astropy.constants import c
from scipy.optimize import curve_fit
import matplotlib
import matplotlib.pylab as plt
from matplotlib import gridspec
from matplotlib.widgets import Button
from numba import jit
from progress.bar import ChargingBar
from multiprocessing import Pool
from os import scandir, getcwd
from os.path import abspath
from scipy.fft import fft
import importlib
from contextlib import contextmanager
import csv
import psutil
import multiprocessing as mp

# ============================================================
# Usar Agg solo para find2c y vgrid
# ============================================================
@contextmanager
def temp_backend(backend_name='Agg'):
    """
    Cambia temporalmente el backend de matplotlib.
    Útil para procesamiento en paralelo que no necesita GUI.
    """
    original_backend = matplotlib.get_backend()
    try:
        matplotlib.use(backend_name)
        # Recargar pyplot para que use el nuevo backend
        import matplotlib.pyplot as plt
        importlib.reload(plt)
        yield
    finally:
        # Restaurar el backend original
        matplotlib.use(original_backend)
        import matplotlib.pyplot as plt
        importlib.reload(plt)

# ============================================================
# Variables globales para eventos de teclado
# ============================================================
yours = None
nmin = True
rv1 = None
rv2 = None
gbase = None
xcent = None

# ============================================================
# Funciones internas auxiliares
# ============================================================
def suggest_nproc(larch, ram_fraction=0.9):
    """
    Sugiere el número óptimo de procesos para Pool basado en RAM disponible,
    número de espectros y tamaño de cada espectro.
    """
    # 1. Obtener RAM total (en bytes)
    mem = psutil.virtual_memory()
    ram_total = mem.total
    ram_limit = mem.available*ram_fraction
    N = len(larch)
    # 3. Leer NAXIS1 del primer espectro
    with fits.open(larch[0], mode='readonly') as hdul:
        naxis1 = hdul[0].header.get('NAXIS1', 0)
    if naxis1 == 0:
        # Si no tiene NAXIS1, estimar con pyasl
        w, f = pyasl.read1dFitsSpec(larch[0])
        naxis1 = len(w)
    # 4. Estimar memoria por proceso hijo
    # spbina requiere 5 matrices en simultaneo
    num_matrices = 5
    mem_per_process = N * naxis1 * 8 * num_matrices
    # 6. Calcular número máximo de procesos por RAM
    if ram_limit <= 0:
        # Si no hay suficiente RAM para ningún proceso, devolver 1
        return 1
    nproc_by_ram = int(ram_limit / mem_per_process)
    if nproc_by_ram < 1:
        nproc_by_ram = 1
    # 7. Limitar por número de núcleos (físicos o lógicos)
    cpu_count = mp.cpu_count()
    nproc = min(nproc_by_ram, int(cpu_count*ram_fraction))
    print(f"Recommended number of processes: {nproc}")
    # 9. Devolver al menos 1
    return max(1, nproc)
    

def plot_find2c_results(obj1, q_array, vector_t, matrix_cc, path1, best_q, jt2, tmed, qmed):
    """
    Genera, guarda y muestra el gráfico de resultados de find2c en modo interactivo.
    """
    import matplotlib.pyplot as plt
    from matplotlib import gridspec
    import numpy as np
    import os

    # plt.ion() ya está activo en find2c, pero por si acaso:
    plt.ion()
    plt.close('all')

    fig = plt.figure(figsize=[6, 7.4])
    gs = gridspec.GridSpec(nrows=2, ncols=2, height_ratios=[7, 1])
    ax0 = fig.add_subplot(gs[0, :], projection='3d')
    ax0.set_title(obj1)
    plt.rc('xtick', labelsize=9)
    plt.rc('ytick', labelsize=9)
    ax0.set_xlabel("Mass ratio ($q$)", fontsize=9)
    ax0.set_ylabel("Temp. [1000 K]", fontsize=9)
    X, Y = np.meshgrid(q_array, vector_t / 1000)
    Z = matrix_cc
    ax0.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap='jet')

    ax1 = fig.add_subplot(gs[1, 0])
    ax1.plot(q_array, qmed, marker='', ls='-', color='blue')
    ax1.set_ylabel("Correlation", fontsize=9)
    ax1.set_xlabel("Mass ratio ($q$)", fontsize=9)
    ax1.set(xlim=(min(q_array), max(q_array)))
    plt.axvline(x=best_q, color='black', linestyle='--', linewidth=1)

    ax2 = fig.add_subplot(gs[1, 1])
    ax2.plot(vector_t, tmed, marker='', ls='-', color='red')
    ax2.set_xlabel("Temp. [1000 K]", fontsize=9)
    ax2.set(xlim=(vector_t[0], vector_t[-1]))
    plt.axvline(x=vector_t[jt2], color='black', linestyle='--', linewidth=1)

    # Guardar imagen
    outfile = os.path.join(path1, obj1 + '_CC.jpeg')
    if os.path.isfile(outfile):
        os.remove(outfile)
    plt.savefig(outfile, dpi=100, bbox_inches='tight')

    # Mostrar la figura (no bloquea porque plt.ion() está activo)
    plt.show()

    # Si quieres que la figura se cierre automáticamente después de un tiempo,
    # podrías usar plt.pause(0.1) o simplemente dejarla abierta.
    # Dejarla abierta permite interactuar (zoom, pan, etc.).

def makelist(lis):
    aux1 = lis.find('@')
    if aux1 == -1:
        path = os.getcwd()
        aux2 = glob.glob(path + '/' + lis)
        newlist = []
        for i in range(len(aux2)):
            newlist.append(aux2[i].replace(path + '/', ''))
    else:
        aux2 = lis.replace('@', '')
        with open(aux2, 'r') as f:
            aux3 = f.readlines()
        newlist = []
        for i in range(len(aux3)):
            newlist.append(aux3[i].rstrip('\n'))
    return newlist

@jit(nopython=True)
def sfcomb(aflux, wei, nit=5, sigma=3):
    nwvl = aflux.shape[1]
    s_mean = np.zeros(nwvl, dtype=np.float64)
    for w in range(nwvl):
        fdata = aflux[:, w].copy()
        weights = wei.copy()
        if np.sum(weights) == 0:
            s_mean[w] = np.mean(fdata)
            continue
        for _ in range(nit):
            total_weight = np.sum(weights)
            if total_weight == 0:
                break
            fmean = np.sum(fdata * weights) / total_weight
            var = np.sum(weights * (fdata - fmean)**2) / total_weight
            fsig = np.sqrt(var)
            if fsig == 0:
                break
            mask = np.abs(fdata - fmean) <= sigma * fsig
            if np.all(mask):
                break
            fdata = fdata[mask]
            weights = weights[mask]
            if len(fdata) < 2:
                break
        if len(weights) > 0 and np.sum(weights) > 0:
            s_mean[w] = np.sum(fdata * weights) / np.sum(weights)
        else:
            s_mean[w] = np.mean(aflux[:, w])
    return s_mean

def setregion(wreg, delta, winf, wsup, amort=0.1):
    reg1 = wreg.split(',')
    reg2 = []
    for i, str1 in enumerate(reg1):
        reg2.append([int(str1.split('-')[0]), int(str1.split('-')[1])])
    reg3 = []
    stat1 = True
    for j, wvx in enumerate(reg2):
        x1 = wvx[0]
        x2 = wvx[1]
        if stat1:
            if x1 >= winf:
                reg3.append(wvx)
                stat1 = False
            elif x1 < winf and x2 > winf:
                wvx[0] = winf
                reg3.append(wvx)
                stat1 = False
            elif x1 < winf and x2 <= winf:
                stat1 = True
        else:
            if x1 > wsup and x2 > wsup:
                break
            elif x1 < wsup and x2 >= wsup:
                wvx[1] = wsup
                reg3.append(wvx)
            elif x1 < wsup and x2 < wsup:
                reg3.append(wvx)
    xmin = reg3[0][0]
    xmax = reg3[-1][1]
    num = int(np.ceil((xmax - xmin) / delta)) + 1
    wvl = np.linspace(xmin, xmin + (num - 2) * delta, num - 1)
    f = np.zeros(len(wvl))
    for k, interv in enumerate(reg3):
        x1 = interv[0]
        x2 = interv[1]
        i1 = np.abs(wvl - x1).argmin(0)
        i2 = np.abs(wvl - x2).argmin(0)
        am2 = amort * (x2 - x1)
        xarr = wvl[i1:i2]
        mask = np.zeros(len(xarr))
        for k, w in enumerate(xarr):
            if w <= (x1 + am2):
                mask[k] = np.sin(np.pi * (w - x1) / (2 * am2))
            elif w > (x1 + am2) and w < (x2 - am2):
                mask[k] = 1
            else:
                mask[k] = np.cos(np.pi * (w - x2 + am2) / (2 * am2))
        f[i1:i2] = mask
    return wvl, f

def continuum(w, f, order=12, type='fit', lo=2, hi=3, nit=10, graph=True):
    from specutils.manipulation import SplineInterpolatedResampler
    spline3 = SplineInterpolatedResampler()
    w_cont = w.copy()
    f_cont = f.copy()
    sigma0 = np.std(f_cont)
    wei = ~np.isnan(f_cont) * 1
    i = 1
    nrej1 = 0
    while i < nit:
        c0 = np.polynomial.chebyshev.Chebyshev.fit(w_cont, f_cont, order, w=wei)(w_cont)
        resid = f_cont - c0
        sigma0 = np.sqrt(np.average((resid)**2, weights=wei))
        wei = 1 * np.logical_and(resid > -lo * sigma0, resid < sigma0 * hi)
        nrej = len(wei) - np.sum(wei)
        if nrej == nrej1:
            break
        nrej1 = nrej
        i = i + 1
    s1 = Spectrum(flux=c0 * u.Jy, spectral_axis=w_cont * 0.1 * u.nm)
    c1 = fit_continuum(s1, model=Chebyshev1D(order), fitter=LinearLSQFitter())
    if type == 'fit':
        fout = c1(w * 0.1 * u.nm).value
    elif type == 'ratio':
        fout = f_cont / c1(w * 0.1 * u.nm).value
    elif type == 'diff':
        fout = f_cont - c1(w * 0.1 * u.nm).value
    if graph:
        fig = plt.figure(figsize=[20, 10])
        ngrid = gridspec.GridSpec(ncols=1, nrows=2, height_ratios=[5, 1])
        ax1 = fig.add_subplot(ngrid[0])
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.ylabel('Flux')
        ax1.plot(w, f, color='gray')
        ax1.plot(w_cont, f_cont, color='blue', linestyle='', marker='.', markersize=2)
        ax1.plot(w_cont, c1(w_cont * 0.1 * u.nm).value, c='r', linestyle='--')
        ax2 = fig.add_subplot(ngrid[1], sharex=ax1, sharey=ax1)
        plt.xlabel('Wavelength [nm]')
        ax2.plot(w, (f - c1(w * 0.1 * u.nm).value), color='gray', linestyle='', marker='.', markersize=1)
        ax2.plot(w_cont, (f_cont - c1(w_cont * 0.1 * u.nm).value), color='blue', linestyle='', marker='.', markersize=1)
        ax2.axhline(y=0, color='red', linestyle='--', linewidth=1)
        plt.tight_layout()
        plt.show()
        plt.close()
    return fout

def on_key(event):
    global yours, nmin, rv1, rv2, gbase, xcent
    yours = event.key
    if event.key == 'y' or event.key == 'Y':
        plt.close()
    elif event.key == 'm' or event.key == 'M':
        if nmin:
            rv1 = event.xdata
            print('RV min: ' + str(round(rv1, 3)) + ' km/s')
            nmin = False
        else:
            rv2 = event.xdata
            print('RV max: ' + str(round(rv2, 3)) + ' km/s')
            nmin = True
            yours = 'r'
            plt.close()
    elif event.key == 'b' or event.key == 'B':
        gbase = event.ydata
        plt.close()
    elif event.key == 'c' or event.key == 'C':
        xcent = event.xdata
        rv1 = xcent - 10
        rv2 = xcent + 10
        plt.close()
    event.canvas.draw()

def fxcor(w, f, wt, ft, mask, fitcont=True, rvcent=None, interac=True):
    from specutils.manipulation import SplineInterpolatedResampler
    spline3 = SplineInterpolatedResampler()
    global yours, nmin, rv1, rv2, gbase, xcent
    nmin = True
    plt.ioff()
    if fitcont:
        fci = continuum(w=w, f=f, type='diff', graph=False)
        fct = continuum(w=wt, f=ft, type='diff', graph=False)
    else:
        fci = f
        fct = ft
    aux_grid = np.log(w)
    dlog = aux_grid[-1] - aux_grid[-2]
    new_log_grid = np.arange(aux_grid[0], aux_grid[-1], dlog)
    aux_mask = Spectrum(flux=mask * u.Jy, spectral_axis=np.log(w) * 0.1 * u.nm)
    aux2_mask = spline3(aux_mask, new_log_grid * 0.1 * u.nm)
    log_mask = splineclean(aux2_mask.flux.value)
    aux_img = Spectrum(flux=fci * u.Jy, spectral_axis=np.log(w) * 0.1 * u.nm)
    aux2_img = spline3(aux_img, new_log_grid * 0.1 * u.nm)
    log_img = splineclean(aux2_img.flux.value)
    fi2 = log_img * log_mask
    aux_sa = Spectrum(flux=fct * u.Jy, spectral_axis=np.log(wt) * 0.1 * u.nm)
    aux2_sa = spline3(aux_sa, new_log_grid * 0.1 * u.nm)
    log_tmp = splineclean(aux2_sa.flux.value)
    ft2 = log_tmp * log_mask
    cc1 = signal.correlate(fi2, ft2, method='fft')
    cc1 = cc1 / (np.sqrt(np.sum(np.power(fi2, 2))) * (np.sqrt(np.sum(np.power(ft2, 2)))))
    i1 = np.where(cc1 == max(cc1))[0][0]
    lamlog1 = new_log_grid - new_log_grid[0]
    lamlog2 = -new_log_grid + new_log_grid[0]
    lamlog2.sort()
    llog = np.concatenate((lamlog2[0:-1], lamlog1), axis=0)
    axisrv = llog * c.value / 1000
    if rvcent is None:
        ccmax = np.argmax(cc1)
        xcent = axisrv[ccmax]
    else:
        xcent = rvcent
    stat = True
    rv1 = xcent - 10
    rv2 = xcent + 10
    gbase = np.mean(cc1)
    nach = 0
    while stat:
        near0 = axisrv.flat[np.abs(axisrv - xcent).argmin()]
        i0 = np.where(axisrv == near0)[0][0]
        near1 = axisrv.flat[np.abs(axisrv - min(rv1, rv2)).argmin()]
        i1 = np.where(axisrv == near1)[0][0]
        near2 = axisrv.flat[np.abs(axisrv - max(rv1, rv2)).argmin()]
        i2 = np.where(axisrv == near2)[0][0]
        try:
            xb = axisrv[i1:i2]
            yb = cc1[i1:i2]
            mb = np.sum(xb * (yb - gbase)) / np.sum(yb - gbase)
            sigb = np.sqrt(np.abs(np.sum((yb - gbase) * (xb - mb)**2) / np.sum(yb - gbase)))
            pb1, pb2 = curve_fit(Gauss, xb, yb - gbase, p0=[np.max(yb - gbase), mb, sigb])
            if nach == 0:
                mod_rv1 = axisrv.flat[np.abs(axisrv - (xcent - pb1[2]*5)).argmin()]
                mod_i1 = np.where(axisrv == mod_rv1)[0][0]
                mod_rv2 = axisrv.flat[np.abs(axisrv - (xcent + pb1[2]*5)).argmin()]
                mod_i2 = np.where(axisrv == mod_rv2)[0][0]
                gbase = min(cc1[mod_i1:mod_i2])
                xb = axisrv[i1:i2]
                yb = cc1[i1:i2]
                mb = np.sum(xb * (yb - gbase)) / np.sum(yb - gbase)
                sigb = np.sqrt(np.abs(np.sum((yb - gbase) * (xb - mb)**2) / np.sum(yb - gbase)))
                pb1, pb2 = curve_fit(Gauss, xb, yb - gbase, p0=[np.max(yb - gbase), mb, sigb])
                nach = 1
            ybgauss = Gauss(xb, *pb1) + gbase
            xrv = pb1[1]
            rverr = np.sqrt(np.diag(pb2))[1]
            sfit = True
            iran = int(2000 * 1000 / (c.value * dlog))
            reg2 = axisrv[i1 - iran:i1 + iran]
            fcc2 = cc1[i1 - iran:i1 + iran]
            ccnew = continuum(w=reg2, f=fcc2, order=4, type='diff', graph=False)
            cc1_neg = ccnew[0:iran]
            cc1_pos = ccnew[iran:iran*2]
            dif_neg = cc1_neg - cc1_pos[::-1]
            dif_pos = cc1_pos - cc1_neg[::-1]
            cc_anti = np.concatenate((dif_neg, dif_pos))
            sig_anti = np.std(cc_anti)
            r = np.max(ccnew) / (np.sqrt(2 * sig_anti))
            wfreq1 = fft(ccnew)
            nmed = int(len(wfreq1.imag) / 4)
            ys2 = np.log(np.sqrt(wfreq1.real[1:nmed]**2 + wfreq1.imag[1:nmed]**2))
            ys2 = ys2 - ys2[-1]
            xs2 = np.arange(1, nmed, 1)
            aux_x = np.concatenate((-1 * xs2[::-1], xs2))
            aux_y = np.concatenate((ys2[::-1], ys2))
            az = np.where(aux_y == 0)
            aux_y[az] = 1e-10
            p1, pc1 = curve_fit(doubleG, aux_x, aux_y)
            B_exp1 = p1[1] * 2.35482 / 2
            err_g1 = len(wfreq1) / (16 * B_exp1 * (1 + r)) + rverr
            serr = str(round(err_g1, 3))
        except Exception:
            xrv = axisrv[np.argmax(cc1)]
            sfit = False
            err_g1 = np.nan
            serr = str(err_g1)
        print('Radial Velocity: ' + str(round(xrv, 3)) + ' +/- ' + serr + ' km/s')
        if interac:
            fig, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 5]})
            plt.setp(a0.get_yticklabels(), visible=False)
            plt.setp(a0.get_xticklabels(), visible=False)
            a0.yaxis.offsetText.set_visible(False)
            a0.plot(axisrv, cc1, color='blue')
            a0.axvline(axisrv[i0], color='red', linestyle='--', linewidth=1)
            a1.set_xlabel("Radial Velocity [km/s]", fontsize=10)
            a1.set_ylabel("Correlation", fontsize=10)
            a1.plot(axisrv, cc1, color='blue')
            a1.set(xlim=(axisrv[i1-10], axisrv[i2+10]))
            ymin2 = min(cc1[i1:i2])
            ymax2 = max(cc1[i1:i2])
            ex2 = (ymax2 - ymin2) * 0.1
            a1.set(ylim=(ymin2 - ex2, ymax2 + ex2))
            if sfit:
                plt.plot(xb, yb, color='black', marker='.', linestyle='')
                plt.plot(xb, ybgauss, color='green', label='fit', linestyle='--')
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
                print('\nProcess completed successfully!\n')
                plt.close()
                break
            elif yours == 'q' or yours == 'Q':
                fsave = False
                plt.close()
                break
        else:
            fsave = True
            print('\nProcess completed successfully!\n')
            break
    plt.close('all')
    return xrv, err_g1, fsave

@jit(nopython=True)
def splineclean(fspl):
    f = fspl.copy()
    n = len(f)
    if n == 0:
        return f
    first_good = 0
    while first_good < n and np.isnan(f[first_good]):
        first_good += 1
    if first_good == n:
        f[:] = 0.0
        return f
    last_good = n - 1
    while last_good >= 0 and np.isnan(f[last_good]):
        last_good -= 1
    if first_good > 0:
        f[0:first_good] = f[first_good]
    if last_good < n - 1:
        f[last_good+1:n] = f[last_good]
    return f

def copyheader(img1, imgout):
    with fits.open(img1, mode='readonly') as hdul:
        with fits.open(imgout, mode='readonly') as hnorm:
            listk = ('CDELT1', 'CTYPE1', 'BUNIT', 'ORIGIN', 'DATE', 'TELESCOP', 'INSTRUME',
                     'OBJECT', 'RA', 'DEC', 'EQUINOX', 'RADECSYS', 'EXPTIME', 'MJD-OBS', 'DATE-OBS', 'UTC', 'LST',
                     'PI-COI', 'CTYPE1', 'CTYPE2', 'ORIGFILE', 'UT', 'ST', 'AIRMASS', 'VRA', 'VRB')
            for i, k in enumerate(listk):
                try:
                    hnorm[0].header[k] = hdul[0].header[k]
                except KeyError:
                    print('Keyword ' + k + ' not found')
            hnorm.flush(output_verify='ignore')

def listmp(ruta=getcwd()):
    return [abspath(arch.path) for arch in scandir(ruta) if arch.is_file()]

def Gauss(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def doubleG(x, a, sigma):
    return a * np.exp(-(x)**2 / (2 * sigma**2))

# ============================================================
# Funciones principales: find2c, hselect, onecomp, qfitg, rvbina,
# rvextract, setrvs, spbina, splot,uniform,vgrid, vexplore
# ============================================================

def find2c(lis, lit, vgamma, spa='A', spb='B', qmin=0.02, qmax=0.5, deltaq=0.01,
           wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900', nproc=None, cpulimit=0.9):

    plt.close('all')
    gc.collect()
    
    from specutils.manipulation import SplineInterpolatedResampler
    spline3 = SplineInterpolatedResampler()
    print('\n\tRunning FIND2C\n')
    VerifyWarning('ignore')
    larch = makelist(lis)
    if lit[len(lit)-1] == '/':
        lit = lit[0:len(lit)-1]
    ltemp = sorted(listmp(lit))
    num_q = int(np.ceil((qmax - qmin) / deltaq)) + 1
    q_array = np.linspace(qmin, qmin + (num_q - 1) * deltaq, num_q)
    path1 = os.getcwd()
    if nproc == None:
        nproc = suggest_nproc(larch,ram_fraction=cpulimit)
    print('Calculating spectra...')
    
    # ========== USAR BACKEND Agg solo para el Pool ==========
    with temp_backend('Agg'):
        with Pool(processes=nproc, maxtasksperchild=1) as pool:
            qa2 = np.array_split(q_array, nproc)
            pres = [pool.apply_async(qparallel, args=[chunk, lis, spa, spb, deltaq, vgamma])
                    for chunk in qa2]
            pres2 = []
            for i, p in enumerate(pres):
                try:
                    result = p.get(timeout=3600)
                    pres2.append(result)
                    gc.collect()
                except mp.TimeoutError:
                    print(f"Timeout en proceso {i}")
                    pool.terminate()
                    pool.join()
                    raise
                except Exception as e:
                    print(f"Error en proceso {i}: {e}")
                    pool.terminate()
                    pool.join()
                    raise
    # ========================================================

    gc.collect()
    plt.close('all')
    print('\nProcess completed successfully!\n')
    print('')

    dlist = []
    with fits.open(larch[0], mode='readonly') as hobs:
        d1 = hobs[0].header['CDELT1']
        dlist.append(d1)
        try:
            obj1 = hobs[0].header['OBJECT']
            obj1 = obj1.replace(' ', '')
        except KeyError:
            obj1 = 'NoObject'
    with fits.open(ltemp[0], mode='readonly') as htmp:
        d2 = htmp[0].header['CDELT1']
        dlist.append(d2)
    waux1, faux1 = pyasl.read1dFitsSpec(larch[0])
    waux2, faux2 = pyasl.read1dFitsSpec(ltemp[0])
    gap = 50
    winf = max(waux1[0], waux2[0]) + gap
    wsup = min(waux1[-1], waux2[-1]) - gap
    new_disp_grid, fmask = setregion(wreg, np.max(dlist), winf, wsup)

    matrix_sq = np.zeros(shape=(len(q_array), len(new_disp_grid)))
    bar1 = ChargingBar('Loading B_q spectra:', max=len(q_array))
    for j, xq in enumerate(q_array):
        aux1 = str(round(xq, len(str(deltaq+1))-2)).replace('.', '')
        wimg, fimg = pyasl.read1dFitsSpec(spb + aux1 + '.fits')
        spec_cont = continuum(wimg, fimg, type='diff', nit=5, lo=2.5, hi=3.5, graph=False)
        aux_img = Spectrum(flux=spec_cont * u.Jy, spectral_axis=wimg * 0.1 * u.nm)
        aux2_img = spline3(aux_img, new_disp_grid * 0.1 * u.nm)
        matrix_sq[j] = splineclean(aux2_img.flux.value) * fmask
        bar1.next()
    bar1.finish()
    print('')

    norma_b = np.max(np.mean(matrix_sq**2, axis=1))
    vector_t = np.zeros(len(ltemp))
    matrix_cc = np.zeros(shape=(len(ltemp), len(q_array)))
    bar2 = ChargingBar('Comparing templates:', max=len(ltemp))
    for k, tmp in enumerate(ltemp):
        with fits.open(tmp, mode='readonly') as htmp:
            teff = htmp[0].header['TEFF']
            vector_t[k] = teff
        wt1, ft1 = pyasl.read1dFitsSpec(tmp)
        temp_cont = continuum(wt1, ft1, order=50, nit=10, type='diff', lo=2, hi=4, graph=False)
        aux_tmp1 = Spectrum(flux=temp_cont * u.Jy, spectral_axis=wt1 * 0.1 * u.nm)
        aux2_tmp1 = spline3(aux_tmp1, new_disp_grid * 0.1 * u.nm)
        template1 = splineclean(aux2_tmp1.flux.value * fmask)
        tt = np.mean(template1**2)
        for l, xq in enumerate(q_array):
            tb = np.mean((template1 * matrix_sq[l]))
            cc = tb / (np.sqrt(norma_b) * np.sqrt(tt))
            matrix_cc[k, l] = cc
        bar2.next()
    bar2.finish()
    print('')

    # Guardar matriz de correlación
    vgindex = str(vgamma)
    auxout = os.path.isfile(obj1 + '_vg_' + vgindex + '.fits')
    if auxout:
        os.remove(obj1 + '_vg_' + vgindex + '.fits')
    fits.writeto(obj1 + '_vg_' + vgindex + '.fits', matrix_cc, overwrite=True)
    with fits.open(obj1 + '_vg_' + vgindex + '.fits', mode='update', verify='ignore') as hcc:
        hcc[0].header['Q0'] = qmin
        hcc[0].header['Q1'] = qmax
        hcc[0].header['DELTA_Q'] = deltaq
        hcc[0].header['T0'] = vector_t[0]
        hcc[0].header['T1'] = vector_t[-1]
        hcc[0].header['DELTA_T'] = vector_t[1] - vector_t[0]
        hcc[0].header['VGAMMA'] = vgamma
        hcc[0].header['WREG'] = '4000-4090,4110-4320,4360-4850,4875-5290,5350-5900'

    # Estimación de q y Teff
    qmed = np.max(matrix_cc, axis=0)
    iq2 = np.where(qmed == np.max(qmed))[0][0]
    if qmed[iq2] != qmed[0] and qmed[iq2] != qmed[-1]:
        best_q = q_array[iq2] - (qmed[iq2+1] - qmed[iq2-1]) * deltaq / 2 / (qmed[iq2-1] + qmed[iq2+1] - 2*qmed[iq2])
    else:
        best_q = q_array[iq2]
    tmed = np.max(matrix_cc, axis=1)
    jt2 = np.argmax(tmed)

    # Imprimir resultados
    print('\t· · · · · · · · · · · · · ·')
    print('\t  Teff=' + str(int(vector_t[jt2])) + ' K\tq = ' + str(round(best_q, 2)) + '  ')
    print('\t  Power=' + str(round(np.max(matrix_cc) / np.mean(matrix_cc) * np.sqrt(len(larch)),2)) + '  ')
    print('\t· · · · · · · · · · · · · ·')

    # Graficar resultados (modo no interactivo)
    plot_find2c_results(obj1, q_array, vector_t, matrix_cc, path1, best_q, jt2, tmed, qmed)
    
def hselect(img, fields):
    VerifyWarning('ignore')
    larch = makelist(img)
    s1 = fields.split(',')
    for i, img in enumerate(larch):
        # CORREGIDO
        with fits.open(img, mode='readonly') as hdul:
            print(img, end='\t')
            for j, ky in enumerate(s1):
                try:
                    print(str(hdul[0].header[ky]), end='\t')
                except KeyError:
                    print(np.nan, end='\t')
            print('\n')

def onecomp(img, lit, wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900'):
    # Crear spline3 localmente
    from specutils.manipulation import SplineInterpolatedResampler
    spline3 = SplineInterpolatedResampler()
    
    plt.ion()
    VerifyWarning('ignore')
    if lit[len(lit)-1] == '/':
        lit = lit[0:len(lit)-1]
    ltemp = sorted(listmp(lit))
    
    # create fading mask
    dlist = []
    # CORREGIDO
    with fits.open(img, mode='readonly') as hobs:
        d1 = hobs[0].header['CDELT1']
        dlist.append(d1)
        try:
            obj1 = hobs[0].header['OBJECT']
            obj1 = obj1.replace(' ', '')
        except KeyError:
            obj1 = 'NoObject'
    
    # CORREGIDO
    with fits.open(ltemp[0], mode='readonly') as htmp:
        d2 = htmp[0].header['CDELT1']
        dlist.append(d2)
    
    waux1, faux1 = pyasl.read1dFitsSpec(img)
    waux2, faux2 = pyasl.read1dFitsSpec(ltemp[0])
    
    # gap expressed wavelength margins in pixels
    gap = 50
    winf = max(waux1[0], waux2[0]) + gap
    wsup = min(waux1[-1], waux2[-1]) - gap
    new_disp_grid, fmask = setregion(wreg, np.max(dlist), winf, wsup)
    
    wimg, fimg = pyasl.read1dFitsSpec(img)
    spec_cont = continuum(wimg, fimg, type='diff', nit=5, lo=2.5, hi=3.5, graph=False)
    aux_img = Spectrum(flux=spec_cont*u.Jy, spectral_axis=wimg*0.1*u.nm)
    aux2_img = spline3(aux_img, new_disp_grid*0.1*u.nm)
    matrix_sq = splineclean(aux2_img.flux.value) * fmask
    
    vector_t = np.zeros(len(ltemp))
    matrix_cc = np.zeros(len(ltemp))
    bar2 = ChargingBar('Comparing templates:', max=len(ltemp))
    for k, tmp in enumerate(ltemp):
        # CORREGIDO
        with fits.open(tmp, mode='readonly') as htmp:
            teff = htmp[0].header['TEFF']
            vector_t[k] = teff
        
        wt1, ft1 = pyasl.read1dFitsSpec(tmp)
        temp_cont = continuum(wt1, ft1, order=50, nit=10, type='diff', lo=2, hi=4, graph=False)
        aux_tmp1 = Spectrum(flux=temp_cont*u.Jy, spectral_axis=wt1*0.1*u.nm)
        aux2_tmp1 = spline3(aux_tmp1, new_disp_grid*0.1*u.nm)
        template1 = splineclean(aux2_tmp1.flux.value * fmask)
        tt = np.mean(template1**2)
        tb = np.mean((template1 * matrix_sq))
        matrix_cc[k] = tb/(np.sqrt(np.mean(matrix_sq**2)) * np.sqrt(tt))
        bar2.next()
    bar2.finish()
    
    plt.figure()
    plt.rc('xtick', labelsize=9)
    plt.rc('ytick', labelsize=9)
    plt.ylabel("Correlation", fontsize=9)
    plt.xlabel("Temp. [1000 K]", fontsize=9)
    plt.plot(vector_t, matrix_cc, marker='', ls='-', color='red')# Usar archivos temporales para evitar conflictos

def qfitg(lista, ordcon=1, output_csv="output.csv", graph=True, save=True):
    """
    Procesa archivos FITS, ajusta gaussiana al perfil qmed (máximo a lo largo del eje Y)
    usando continuum para el fondo y Gauss para el pico.
    Guarda resultados en CSV y muestra gráficas interactivas.
    """
    resultados = []
    larch = makelist(lista)
    for i, img in enumerate(larch):
        print(f"Procesando {i+1}/{len(larch)}: {img}")
        with fits.open(img, mode='update', verify='ignore') as hdul:
            if i == 0:
                qmin = hdul[0].header['Q0']
                qmax = hdul[0].header['Q1']
                deltaq = hdul[0].header['DELTA_Q']
                num_q = int(np.ceil((qmax - qmin) / deltaq)) + 1
                q_array = np.linspace(qmin, qmin + (num_q - 1) * deltaq, num_q)
            matrix_cc = hdul[0].data
            qmed = np.max(matrix_cc, axis=0)  # perfil

        # 1. Estimar el continuo (fondo)
        fondo = continuum(q_array, qmed, order=ordcon, type='fit', lo=3, hi=3, nit=5, graph=False)

        # 2. Restar fondo para obtener solo el pico
        pico = qmed - fondo

        # 3. Ajuste de gaussiana sin offset
        idx_max = np.argmax(pico)
        max_val = pico[idx_max]
        # Estimar sigma inicial a partir de la anchura a media altura
        half_max = max_val / 2
        left = idx_max
        while left > 0 and pico[left] > half_max:
            left -= 1
        right = idx_max
        while right < len(pico) - 1 and pico[right] > half_max:
            right += 1
        fwhm = right - left
        sigma_guess = fwhm / 2.355 if fwhm > 1 else 5

        # Ventana de ajuste (5 sigma alrededor del pico)
        radio = int(5 * sigma_guess)
        radio = min(radio, idx_max, len(pico) - idx_max - 1)
        left_fit = max(0, idx_max - radio)
        right_fit = min(len(pico), idx_max + radio + 1)
        x_fit = q_array[left_fit:right_fit]
        y_fit = pico[left_fit:right_fit]

        try:
            popt, _ = curve_fit(Gauss, x_fit, y_fit, p0=[max_val, q_array[idx_max], sigma_guess])
            amp, mu, sigma = popt
        except Exception as e:
            print(f"Error en ajuste gaussiano para {img}: {e}. Usando estimaciones.")
            amp = max_val
            mu = q_array[idx_max]
            sigma = sigma_guess
        
        sigma = abs(sigma)
        # Evaluar gaussiana ajustada en todo el rango
        gauss_ajustada = Gauss(q_array, amp, mu, sigma)

        # 4. Estimar el ruido usando puntos fuera del pico (fuera de 3 sigma)
        mask = np.abs(q_array - mu) > 3 * sigma
        puntos_ruido = qmed[mask]  # datos originales sin restar fondo
        if len(puntos_ruido) < 10:
            # Si no hay suficientes, tomar colas del array
            n = len(qmed)
            n_noise = max(10, int(0.2 * n))
            cola_izq = qmed[:n_noise//2]
            cola_der = qmed[-n_noise//2:]
            puntos_ruido = np.concatenate([cola_izq, cola_der])
        ruido_std = np.std(puntos_ruido)

        # 5. Guardar resultados
        idx_mu = np.argmin(np.abs(q_array - mu))
        offset = fondo[idx_mu]

        resultados.append({
            'archivo': img,
            'peak': max_val,
            'amp_gauss': amp,
            'sigma': sigma,
            'offset': offset,
            'ruido': ruido_std
        })

        # 6. Graficar (modo mejorado)
        if graph:
            plt.figure(figsize=(12, 6))
            plt.plot(q_array, qmed, 'b-', alpha=0.7, label='Calculated profile')
            plt.plot(q_array, fondo, 'g-', linewidth=1.5, label='Noise level')
            plt.plot(q_array, gauss_ajustada + fondo, 'r--', linewidth=2, label='Gauss fitting')
            plt.xlabel('q')
            plt.ylabel('Power')
            plt.title(f'File: {img}\nPeak={max_val:.4f}, amp_gauss={amp:.4f}, sigma={sigma:.4f}, noise={ruido_std:.4f}')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show(block=False)
            print("\nPress Enter to continue to next file...")
            input()  # Espera a que el usuario presione Enter
            plt.close()

    if save:
        # 7. Escribir CSV
        with open(output_csv, 'w', newline='') as csvfile:
            fieldnames = ['archivo', 'peak', 'amp_gauss', 'sigma', 'offset', 'ruido']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for row in resultados:
                writer.writerow(row)

        print(f"Save file")
    return resultados


def spbina(lis, spa='A', spb='B', nit=5, frat=0.01, cord=10, reject=True, q=None, vgamma=None, obspha=False, showtit=True):

    if showtit:
        print('')
        print('\t  Running SPBINA')
        print('')
    
    from specutils.manipulation import SplineInterpolatedResampler
    spline3 = SplineInterpolatedResampler()
    VerifyWarning('ignore')
    
    larch = makelist(lis)
    nimg = len(larch)
    
    with fits.open(larch[0], mode='readonly') as haux:
        delta = haux[0].header['CDELT1']
    
    xwmin = []
    xwmax = []
    for i in range(nimg):
        wbusq, dumm2 = pyasl.read1dFitsSpec(larch[i])
        xwmin.append(wbusq[0])
        xwmax.append(wbusq[-1])
    xmin = np.max(xwmin)
    xmax = np.min(xwmax)
    num = int(np.ceil((xmax - xmin) / delta)) + 1
    new_disp_grid = np.linspace(xmin, xmin + (num - 2) * delta, num - 1)
    nwvl = len(new_disp_grid)
    
    baux1 = os.path.isfile(spb + '.fits')
    
    if not baux1:
        fa = 1.0 / (1.0 + frat)
        fb = frat * fa
        tmp_matrix = np.zeros(shape=(nimg, nwvl))
        cont = 0
        for img in larch:
            wimg, fimg = pyasl.read1dFitsSpec(img)
            aux_img = Spectrum(flux=fimg * u.Jy, spectral_axis=wimg * 0.1 * u.nm)
            aux2_img = spline3(aux_img, new_disp_grid * 0.1 * u.nm)
            tmp = aux2_img.flux.value
            tmp_matrix[cont] = splineclean(tmp)
            cont += 1
        f_mean = np.zeros(nwvl)
        for i in range(nwvl):
            f_mean[i] = np.mean(tmp_matrix[:, i])
        f_cont = continuum(new_disp_grid, f_mean, order=cord, nit=10, type='fit', graph=False)
        f_cont = splineclean(f_cont)
        B = f_cont * fb
        del tmp_matrix
    else:
        waux, faux = pyasl.read1dFitsSpec(spb + '.fits')
        aux_B = Spectrum(flux=faux * u.Jy, spectral_axis=waux * 0.1 * u.nm)
        aux2_B = spline3(aux_B, new_disp_grid * 0.1 * u.nm)
        tmp2 = aux2_B.flux.value
        B = splineclean(tmp2)
    
    obs_matrix = np.zeros(shape=(nimg, nwvl))
    dsA_matrix = np.zeros(shape=(nimg, nwvl))
    dsB_matrix = np.zeros(shape=(nimg, nwvl))
    za_matrix = np.zeros(shape=(nimg, nwvl))
    zb_matrix = np.zeros(shape=(nimg, nwvl))
    cont = 0
    vra_array = np.zeros(nimg)
    vrb_array = np.zeros(nimg)
    
    for idx_img, img in enumerate(larch):
        with fits.open(img, mode='readonly') as hdul:
            vra = hdul[0].header['VRA']
            if q is None:
                vrb = hdul[0].header['VRB']
            elif q > 0:
                if vgamma is None:
                    vgamma = 0
                vrb = vgamma - (vra - vgamma) / q
        vra_array[cont] = vra
        vrb_array[cont] = vrb
        wimg, fimg = pyasl.read1dFitsSpec(img)
        
        wlprime_B = new_disp_grid * np.sqrt((1. + vrb / 299792.458) / (1. - vrb / 299792.458))
        if not np.all(np.diff(wlprime_B) > 0):
            idx = np.argsort(wlprime_B)
            wlprime_B = wlprime_B[idx]
        aux_sb = Spectrum(flux=B * u.Jy, spectral_axis=wlprime_B * 0.1 * u.nm)
        aux2_sb = spline3(aux_sb, new_disp_grid * 0.1 * u.nm)
        fb_dop = aux2_sb.flux.value
        
        aux_img = Spectrum(flux=fimg * u.Jy, spectral_axis=wimg * 0.1 * u.nm)
        aux2_img = spline3(aux_img, new_disp_grid * 0.1 * u.nm)
        dsB = aux2_img.flux.value - fb_dop
        dsB_matrix[cont] = splineclean(dsB)
        obs_matrix[cont] = splineclean(aux2_img.flux.value)
        cont += 1
    
    wei = np.mean(obs_matrix, axis=1) / np.max(np.mean(obs_matrix, axis=1))

    for i in range(nit):
        if showtit and i == 0:
            bar3 = ChargingBar('Calculating spectra:', max=nit)
        for j in range(nimg):
            wlprime_A = new_disp_grid * np.sqrt((1. - vra_array[j] / 299792.458) / (1. + vra_array[j] / 299792.458))
            if not np.all(np.diff(wlprime_A) > 0):
                idx = np.argsort(wlprime_A)
                wlprime_A = wlprime_A[idx]
            aux_sa = Spectrum(flux=dsB_matrix[j] * u.Jy, spectral_axis=wlprime_A * 0.1 * u.nm)
            aux2_sa = spline3(aux_sa, new_disp_grid * 0.1 * u.nm)
            fa_dop = aux2_sa.flux.value
            za_matrix[j] = splineclean(fa_dop)
        
        if reject:
            A = sfcomb(za_matrix, wei)
        else:
            A = np.average(za_matrix, axis=0, weights=wei)
        
        for j in range(nimg):
            wlprime_A = new_disp_grid * np.sqrt((1. + vra_array[j] / 299792.458) / (1. - vra_array[j] / 299792.458))
            if not np.all(np.diff(wlprime_A) > 0):
                idx = np.argsort(wlprime_A)
                wlprime_A = wlprime_A[idx]
            aux_sa = Spectrum(flux=A * u.Jy, spectral_axis=wlprime_A * 0.1 * u.nm)
            aux2_sa = spline3(aux_sa, new_disp_grid * 0.1 * u.nm)
            fa_dop = aux2_sa.flux.value
            fa_dop = splineclean(fa_dop)
            dsA_matrix[j] = obs_matrix[j] - fa_dop
            
            wlprime_B = new_disp_grid * np.sqrt((1. - vrb_array[j] / 299792.458) / (1. + vrb_array[j] / 299792.458))
            if not np.all(np.diff(wlprime_B) > 0):
                idx = np.argsort(wlprime_B)
                wlprime_B = wlprime_B[idx]
            aux_sb = Spectrum(flux=dsA_matrix[j] * u.Jy, spectral_axis=wlprime_B * 0.1 * u.nm)
            aux2_sb = spline3(aux_sb, new_disp_grid * 0.1 * u.nm)
            fb_dop = aux2_sb.flux.value
            zb_matrix[j] = splineclean(fb_dop)
        
        if reject:
            B = sfcomb(zb_matrix, wei)
        else:
            B = np.average(zb_matrix, axis=0, weights=wei)
        

        for j in range(nimg):
            wlprime_B = new_disp_grid * np.sqrt((1. + vrb_array[j] / 299792.458) / (1. - vrb_array[j] / 299792.458))
            if not np.all(np.diff(wlprime_B) > 0):
                idx = np.argsort(wlprime_B)
                wlprime_B = wlprime_B[idx]
            aux_sb = Spectrum(flux=B * u.Jy, spectral_axis=wlprime_B * 0.1 * u.nm)
            aux2_sb = spline3(aux_sb, new_disp_grid * 0.1 * u.nm)
            fb_dop = aux2_sb.flux.value
            dsB_matrix[j] = obs_matrix[j] - splineclean(fb_dop)
        
        if showtit:
            bar3.next()
    
    if showtit:
        bar3.finish()

    # Usar archivos temporales para evitar conflictos
    temp_spa = spa + '.tmp'
    temp_spb = spb + '.tmp'
    pyasl.write1dFitsSpec(temp_spa, A, wvl=new_disp_grid, clobber=True)
    pyasl.write1dFitsSpec(temp_spb, B, wvl=new_disp_grid, clobber=True)
    os.rename(temp_spa, spa + '.fits')
    os.rename(temp_spb, spb + '.fits')
    
    if obspha:
        for i, img in enumerate(larch):
            pyasl.write1dFitsSpec('ds-B_' + img, dsB_matrix[i], wvl=new_disp_grid, clobber=True)
            copyheader(img, 'ds-B_' + img)
            pyasl.write1dFitsSpec('ds-A_' + img, dsA_matrix[i], wvl=new_disp_grid, clobber=True)
            copyheader(img, 'ds-A_' + img)
    
    # Liberar grandes arrays
    for var_name in ['obs_matrix', 'dsA_matrix', 'dsB_matrix', 'za_matrix', 'zb_matrix']:
        try:
            del locals()[var_name]
        except:
            pass
    try:
        del spline3
    except:
        pass
    gc.collect()

def splot(file, xmin=None, xmax=None, ymin=None, ymax=None, scale=1., markpix=False, newfig=True, color='r'):
    plt.ion()
    w, f = pyasl.read1dFitsSpec(file)
    if newfig:
        plt.figure(figsize=[20, 10])
    if xmin is None:
        x1 = np.min(w)
    else:
        x1 = xmin
    if xmax is None:
        x2 = np.max(w)
    else:
        x2 = xmax
    if ymin is None:
        y1 = np.min(f) * scale
    else:
        y1 = ymin
    if ymax is None:
        y2 = np.max(f) * scale
    else:
        y2 = ymax * scale
    plt.axis([x1, x2, y1, y2])
    plt.ylabel('Flux')
    plt.xlabel('Wavelength')
    plt.title(file)
    plt.plot(w, f * scale, marker='', color=color, linewidth=1)
    if markpix:
        plt.plot(w, f * scale, marker='.', markersize=2, color='black', linestyle='')
    plt.tight_layout()
    plt.show()
    return w, f

def setrvs(lis, ta='templateA', tb=None, wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900',
           fitcont=True, interac=True):
    spline3 = SplineInterpolatedResampler()
    plt.ioff()
    VerifyWarning('ignore')
    larch = makelist(lis)
    wmins = []
    wmaxs = []
    deltamin = []
    wta, fta = None, None
    if ta is not None:
        ta = ta.replace('.fits', '')
        if os.path.isfile(ta + '.fits'):
            wta, fta = pyasl.read1dFitsSpec(ta + '.fits')
            wmins.append(wta[0])
            wmaxs.append(wta[-1])
            with fits.open(ta + '.fits', mode='readonly') as tadul:
                de1 = tadul[0].header['CDELT1']
                deltamin.append(de1)
    wtb, ftb = None, None
    if tb is not None:
        tb = tb.replace('.fits', '')
        if os.path.isfile(tb + '.fits'):
            wtb, ftb = pyasl.read1dFitsSpec(tb + '.fits')
            wmins.append(wtb[0])
            wmaxs.append(wtb[-1])
            with fits.open(tb + '.fits', mode='readonly') as tbdul:
                de2 = tbdul[0].header['CDELT1']
                deltamin.append(de2)
    for k, img in enumerate(larch):
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
        print('Set RVs for ' + img + '\n')
        wimg, fimg = pyasl.read1dFitsSpec(img)
        with fits.open(img, mode='update', verify='ignore') as hdul:
            if k == 0:
                de3 = hdul[0].header['CDELT1']
                deltamin.append(de3)
                wmins.append(wimg[0])
                wmaxs.append(wimg[-1])
                winf = max(wmins)
                wsup = min(wmaxs)
                delta = np.max(deltamin)
                new_disp_grid, fmask = setregion(wreg, delta, winf, wsup)
            aux_img = Spectrum(flux=fimg * u.Jy, spectral_axis=wimg * 0.1 * u.nm)
            aux2_img = spline3(aux_img, new_disp_grid * 0.1 * u.nm)
            fnew = aux2_img.flux.value
            if ta is not None and wta is not None:
                print('Cross-correlation for primary component\n')
                best_vra, aerr, asave = fxcor(new_disp_grid, fnew, wta, fta, fmask, fitcont=fitcont, interac=interac)
                if asave:
                    hdul[0].header['VRA'] = round(best_vra, 6)
            if tb is not None and wtb is not None:
                print('\nCross-correlation for secondary component\n')
                best_vrb, berr, bsave = fxcor(new_disp_grid, fnew, wtb, ftb, fmask, fitcont=fitcont, interac=interac)
                if bsave:
                    hdul[0].header['VRB'] = round(best_vrb, 6)
    plt.close('all')

def rvbina(lis, spa='A', spb='B', ta='templateA', tb='templateB',
           wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900',
           aconv=0.5, keyjd='MJD-OBS', fitcont=True, interac=True):
    spline3 = SplineInterpolatedResampler()
    plt.ioff()
    VerifyWarning('ignore')
    larch = makelist(lis)
    ta = ta.replace('.fits', '')
    tb = tb.replace('.fits', '')
    aaux1 = os.path.isfile(ta + '.fits')
    baux1 = os.path.isfile(tb + '.fits')
    if not aaux1:
        print('Can not access to primary template spectrum')
        print('END')
        return
    if not baux1:
        print('Can not access to secondary template spectrum')
        print('END')
        return
    aaux2 = os.path.isfile(spa + '.fits')
    baux2 = os.path.isfile(spb + '.fits')
    if not aaux2:
        print('Can not access to ' + spa + '.fits')
        print('END')
        return
    if not baux2:
        print('Can not access to ' + spb + '.fits')
        print('END')
        return
    if aaux1 and baux1 and aaux2 and baux2:
        for k, img in enumerate(larch):
            print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::')
            print('Processing ' + img + '...\n')
            wimg, fimg = pyasl.read1dFitsSpec(img)
            with fits.open(img, mode='update', verify='ignore') as hdul:
                if k == 0:
                    wa, fa = pyasl.read1dFitsSpec(spa + '.fits')
                    wb, fb = pyasl.read1dFitsSpec(spb + '.fits')
                    wta, fta = pyasl.read1dFitsSpec(ta + '.fits')
                    with fits.open(ta + '.fits', mode='readonly') as tadul:
                        del_a = tadul[0].header['CDELT1']
                    wtb, ftb = pyasl.read1dFitsSpec(tb + '.fits')
                    with fits.open(tb + '.fits', mode='readonly') as tbdul:
                        del_b = tbdul[0].header['CDELT1']
                    del_img = hdul[0].header['CDELT1']
                    delta = np.max([del_a, del_b, del_img])
                    winf = max(wimg[0], wta[0], wtb[0])
                    wsup = min(wimg[-1], wta[-1], wtb[-1])
                    new_disp_grid, fmask = setregion(wreg, delta, winf, wsup)
                    try:
                        obj1 = hdul[0].header['OBJECT']
                    except KeyError:
                        obj1 = 'NoObject'
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
            if vra is None:
                print('RV for primary component not found.')
                setrvs(img, ta=ta, wreg=wreg, fitcont=fitcont)
                with fits.open(img, mode='readonly') as hdul:
                    vra = hdul[0].header['VRA']
            if vrb is None:
                print('RV for secondary component not found.')
                setrvs(img, tb=tb, wreg=wreg, fitcont=fitcont)
                with fits.open(img, mode='readonly') as hdul:
                    vrb = hdul[0].header['VRB']
            wlprime_A = wa * np.sqrt((1. + vra / 299792.458) / (1. - vra / 299792.458))
            if not np.all(np.diff(wlprime_A) > 0):
                idx = np.argsort(wlprime_A)
                wlprime_A = wlprime_A[idx]
            wlprime_B = wb * np.sqrt((1. + vrb / 299792.458) / (1. - vrb / 299792.458))
            if not np.all(np.diff(wlprime_B) > 0):
                idx = np.argsort(wlprime_B)
                wlprime_B = wlprime_B[idx]
            aux_img = Spectrum(flux=fimg * u.Jy, spectral_axis=wimg * 0.1 * u.nm)
            aux_sa = Spectrum(flux=fa * u.Jy, spectral_axis=wlprime_A * 0.1 * u.nm)
            aux_sb = Spectrum(flux=fb * u.Jy, spectral_axis=wlprime_B * 0.1 * u.nm)
            aux2_img = spline3(aux_img, new_disp_grid * 0.1 * u.nm)
            aux2_sa = spline3(aux_sa, new_disp_grid * 0.1 * u.nm)
            aux2_sb = spline3(aux_sb, new_disp_grid * 0.1 * u.nm)
            dsA = aux2_img.flux.value - aux2_sa.flux.value
            dsA = splineclean(dsA)
            dsB = aux2_img.flux.value - aux2_sb.flux.value
            dsB = splineclean(dsB)
            print('Cross-correlation for primary component\n')
            best_vra, aerr, asave = fxcor(new_disp_grid, dsB, wta, fta, fmask, fitcont=fitcont, interac=interac)
            print('Cross-correlation for secondary component\n')
            best_vrb, berr, bsave = fxcor(new_disp_grid, dsA, wtb, ftb, fmask, fitcont=fitcont, interac=interac)
            with fits.open(img, mode='update', verify='ignore') as hdul:
                if asave:
                    new_vra = best_vra * aconv + vra * (1 - aconv)
                    err_a = np.abs(new_vra - vra)
                    hdul[0].header['VRA'] = round(new_vra, 6)
                    print(img + ', VRA: ' + str(round(vra, 3)) + ' km/s --> ' + str(round(new_vra, 3)) + ' km/s')
                else:
                    new_vra = vra
                    err_a = 0
                if bsave:
                    new_vrb = best_vrb * aconv + vrb * (1 - aconv)
                    err_b = np.abs(new_vrb - vrb)
                    hdul[0].header['VRB'] = round(new_vrb, 6)
                    print(img + ', VRB: ' + str(round(vrb, 3)) + ' km/s --> ' + str(round(new_vrb, 3)) + ' km/s')
                else:
                    new_vrb = vrb
                    err_b = 0
            name = img.replace('.fits', '')
            if not os.path.isfile(name + '.log'):
                with open(name + '.log', 'w') as flog:
                    flog.write('#RV_A\te_it_A\tRV_B\te_it_B\n')
            with open(name + '.log', 'a') as flog:
                flog.write(str(round(new_vra, 6)) + '\t' + str(round(err_a, 6)) + '\t' +
                          str(round(new_vrb, 6)) + '\t' + str(round(err_b, 6)) + '\n')
            if not os.path.isfile(obj1 + '_RV.txt'):
                with open(obj1 + '_RV.txt', 'w') as frv:
                    frv.write('#JD\tRV_A\te_A\tRV_B\te_B\n')
            with open(obj1 + '_RV.txt', 'a') as frv:
                frv.write(str(xjd) + '\t' + str(round(new_vra, 6)) + '\t' +
                         str(round(err_a, 6)) + '\t' + str(round(new_vrb, 6)) + '\t' + str(round(err_b, 6)) + '\n')
    plt.close('all')

def uniform(lis, wmin=None, wmax=None, keyrv='VRA', cmcrem=True, nit=5, lo=5, hi=2, interac=True):
    spline3 = SplineInterpolatedResampler()
    global yours
    plt.close()
    plt.ioff()
    VerifyWarning('ignore')
    larch = makelist(lis)
    with fits.open(larch[0], mode='readonly') as hobs:
        dx = hobs[0].header['CDELT1']
    xwmin = []
    xwmax = []
    for i, imin in enumerate(larch):
        wbusq, dumm2 = pyasl.read1dFitsSpec(imin)
        xwmin.append(wbusq[0])
        xwmax.append(wbusq[-1])
    xmin = np.max(xwmin)
    xmax = np.min(xwmax)
    num = int(np.ceil((xmax - xmin) / dx)) + 1
    new_disp_grid = np.linspace(xmin, xmin + (num - 2) * dx, num - 1)
    nwvl = len(new_disp_grid)
    obs_matrix = np.zeros(shape=(len(larch), nwvl))
    ldel = []
    print('\n\tProcessing spectra list...\n')
    larch2 = larch.copy()
    for k, img in enumerate(larch):
        wimg, fimg = pyasl.read1dFitsSpec(img)
        rp = np.where(fimg == 0)
        fimg[rp] = np.mean(fimg)
        with fits.open(img, mode='readonly') as hdul:
            vra = hdul[0].header[keyrv]
        w2 = wimg * np.sqrt((1. - vra / 299792.458) / (1. + vra / 299792.458))
        aux_img = Spectrum(flux=fimg * u.Jy, spectral_axis=w2 * 0.1 * u.nm)
        aux_img2 = spline3(aux_img, new_disp_grid * 0.1 * u.nm)
        if interac:
            fig = plt.figure(figsize=[15, 8])
            plt.title('Spectrum ' + img)
            plt.xlabel('Wavelength [A]')
            plt.ylabel('Flux')
            plt.plot(wimg, fimg, color='red', ls='-')
            plt.plot(new_disp_grid, aux_img2.flux.value, ls='--', marker='', color='blue')
            plt.legend(('Original', 'RV corrected'))
            plt.tight_layout()
            fig.canvas.mpl_connect('key_press_event', on_key)
            print('\t............................................')
            print('\t:      Press Y to save and continue        :')
            print('\t:  Press any button to discard spectrum    :')
            print('\t············································')
            plt.show()
            if yours == 'y' or yours == 'Y':
                obs_matrix[k] = splineclean(aux_img2.flux.value)
            else:
                ldel.append(k)
                larch2.remove(img)
        else:
            obs_matrix[k] = splineclean(aux_img2.flux.value)
    obs_matrix = np.delete(obs_matrix, ldel, axis=0)
    wei = np.mean(obs_matrix, axis=1) / np.max(np.mean(obs_matrix, axis=1))
    smean1 = sfcomb(obs_matrix, wei)
    if cmcrem:
        print('Removing cosmic rays...')
        for i in range(nit):
            cont = 0
            print('Iteration #' + str(i))
            for j in np.arange(len(new_disp_grid)):
                fdev = np.std(obs_matrix[:, j])
                for k in np.arange(len(obs_matrix)):
                    fx = obs_matrix[k, j]
                    fcomp = np.mean(obs_matrix[:, j])
                    if (fx < fcomp - fdev * lo) or (fx > fcomp + fdev * hi):
                        obs_matrix[k, j] = fcomp
                        cont += 1
            print('\tReplaced ' + str(cont) + ' pixels')
            smean1 = sfcomb(obs_matrix, wei)
    f_cont = continuum(new_disp_grid, smean1, type='fit', graph=False)
    for i, img in enumerate(larch2):
        with fits.open(img, mode='readonly') as hobs:
            vra = hobs[0].header[keyrv]
        grid_aux = new_disp_grid * np.sqrt((1. + vra / 299792.458) / (1. - vra / 299792.458))
        f_norm = obs_matrix[i] / f_cont
        if wmin is not None:
            i_min = np.searchsorted(new_disp_grid, wmin, side='left')
        else:
            i_min = 0
        if wmax is not None:
            i_max = np.searchsorted(new_disp_grid, wmax, side='right') - 1
        else:
            i_max = len(new_disp_grid) - 1
        naux = img.replace('.fits', '')
        print('')
        print('\tSaving ' + img + '...')
        header_nuevo = fits.Header()
        keywords_esenciales = [
            'SIMPLE', 'BITPIX', 'NAXIS', 'NAXIS1', 'EXTEND', 'OBJECT', 'CRPIX1',
            'CRVAL1', 'CDELT1', 'CTYPE1', 'BUNIT', 'CUNIT1', 'RA', 'DEC', 'EQUINOX',
            'DATE-OBS', 'MJD-OBS', 'EXPTIME', 'TELESCOP', 'INSTRUME',
            'RADECSYS', 'ORIGIN', 'DATE', 'AIRMASS', 'UTC', 'LST', 'HJD', 'CD1_1', keyrv]
        for kw in keywords_esenciales:
            if kw in hobs[0].header:
                try:
                    value = hobs[0].header[kw]
                    if isinstance(value, str):
                        value = ''.join(char for char in value if ord(char) < 128)
                    header_nuevo[kw] = value
                except:
                    pass
        header_nuevo['NAXIS'] = 1
        header_nuevo['NAXIS1'] = len(new_disp_grid[i_min:i_max])
        header_nuevo['CRVAL1'] = new_disp_grid[i_min]
        header_nuevo['CRPIX1'] = 1.0
        header_nuevo['CTYPE1'] = 'WAVELEN'
        header_nuevo['CUNIT1'] = 'Angstrom'
        pyasl.write1dFitsSpec(naux + '_CMC.fits', f_norm[i_min:i_max], grid_aux[i_min:i_max], clobber=True, header=header_nuevo)
    plt.close('all')


def vgrid(lis, lit, svmin=-1, svmax=1, step=0.1, qmin=0.02, qmax=0.5, deltaq=0.01,
          wreg='4000-4090,4110-4320,4360-4850,4875-5290,5350-5900', nproc=None, cpulimit=0.9):
    plt.close('all')
    gc.collect()
    
    from specutils.manipulation import SplineInterpolatedResampler
    spline3 = SplineInterpolatedResampler()
    istr = '00'
    if not os.path.isdir('output_00'):
        os.mkdir('output_00')
    else:
        it = 1
        while os.path.isdir('output_' + istr):
            if it < 10:
                istr = '0' + str(it)
            else:
                istr = str(it)
            it += 1
        os.mkdir('output_' + istr)
    outfolder = 'output_' + istr
    print('Output folder: ' + outfolder)
    svrange = np.arange(svmin, svmax + step, step)
    if math.modf(step)[0] == 0:
        nrd = 0
    else:
        nrd = len(str(step)) - str(step).find('.') - 1
    plt.ion()
    print('\n\t  Running VGRID\n')
    VerifyWarning('ignore')
    larch = makelist(lis)
    if nproc == None:
        nproc = suggest_nproc(larch,ram_fraction=cpulimit)
    if lit[len(lit)-1] == '/':
        lit = lit[0:len(lit)-1]
    ltemp = sorted(listmp(lit))
    num_q = int(np.ceil((qmax - qmin) / deltaq)) + 1
    q_array = np.linspace(qmin, qmin + (num_q - 1) * deltaq, num_q)
    dlist = []
    with fits.open(larch[0], mode='readonly') as hobs:
        d1 = hobs[0].header['CDELT1']
        dlist.append(d1)
        try:
            obj1 = hobs[0].header['OBJECT']
            obj1 = obj1.replace(' ', '')
        except KeyError:
            obj1 = 'NoObject'
    with fits.open(ltemp[0], mode='readonly') as htmp:
        d2 = htmp[0].header['CDELT1']
        dlist.append(d2)
    waux1, faux1 = pyasl.read1dFitsSpec(larch[0])
    waux2, faux2 = pyasl.read1dFitsSpec(ltemp[0])
    gap = 50
    winf = max(waux1[0], waux2[0]) + gap
    wsup = min(waux1[-1], waux2[-1]) - gap
    new_disp_grid, fmask = setregion(wreg, np.max(dlist), winf, wsup)
    vector_t = np.zeros(len(ltemp))
    tt_array = np.zeros(len(ltemp))
    matrix_tmp = np.zeros(shape=(len(ltemp), len(new_disp_grid)))
    bar2 = ChargingBar('Loading templates:', max=len(ltemp))
    for k, tmp in enumerate(ltemp):
        with fits.open(tmp, mode='readonly') as htmp:
            teff = htmp[0].header['TEFF']
            vector_t[k] = teff
        wt1, ft1 = pyasl.read1dFitsSpec(tmp)
        temp_cont = continuum(wt1, ft1, order=50, nit=10, type='diff', lo=2, hi=4, graph=False)
        aux_tmp1 = Spectrum(flux=temp_cont * u.Jy, spectral_axis=wt1 * 0.1 * u.nm)
        aux2_tmp1 = spline3(aux_tmp1, new_disp_grid * 0.1 * u.nm)
        template1 = splineclean(aux2_tmp1.flux.value * fmask)
        matrix_tmp[k] = template1
        tt_array[k] = np.mean(template1**2)
        bar2.next()
    bar2.finish()
    print('')
    bar0 = ChargingBar('Calculating syst. RV:', max=len(svrange))
    
    for vgamma in svrange:
        matrix_cc = np.zeros(shape=(len(ltemp), len(q_array)))
        vgindex = str(round(vgamma, nrd))
        for xq in q_array:
            aux1 = str(round(xq, len(str(deltaq+1))-2)).replace('.', '')
            if os.path.isfile('B' + aux1 + '.fits'):
                os.remove('B' + aux1 + '.fits')
            if os.path.isfile('A' + aux1 + '.fits'):
                os.remove('A' + aux1 + '.fits')
        
        # ========== USAR BACKEND Agg solo para el Pool ==========
        with temp_backend('Agg'):
            with Pool(processes=nproc, maxtasksperchild=1) as pool:
                qa2 = np.array_split(q_array, nproc)
                pres = [pool.apply_async(qparallel, args=[chunk, lis, 'A', 'B', deltaq, vgamma])
                        for chunk in qa2]
                for i, p in enumerate(pres):
                    try:
                        p.get(timeout=3600)
                    except mp.TimeoutError:
                        print(f"Timeout en proceso {i}")
                        pool.terminate()
                        pool.join()
                        raise
                    except Exception as e:
                        print(f"Error en proceso {i}: {e}")
                        pool.terminate()
                        pool.join()
                        raise
        # ========================================================
        
        matrix_sq = np.zeros(shape=(len(q_array), len(new_disp_grid)))
        for j, xq in enumerate(q_array):
            aux1 = str(round(xq, len(str(deltaq+1))-2)).replace('.', '')
            wimg, fimg = pyasl.read1dFitsSpec('B' + aux1 + '.fits')
            spec_cont = continuum(wimg, fimg, type='diff', lo=2.5, hi=3.5, graph=False)
            aux_img = Spectrum(flux=spec_cont * u.Jy, spectral_axis=wimg * 0.1 * u.nm)
            aux2_img = spline3(aux_img, new_disp_grid * 0.1 * u.nm)
            matrix_sq[j] = splineclean(aux2_img.flux.value) * fmask
        norma_b = np.max(np.mean(matrix_sq**2, axis=1))
        for k, tt in enumerate(tt_array):
            for l, xq in enumerate(q_array):
                tb = np.mean(matrix_tmp[k] * matrix_sq[l])
                cc = tb / (np.sqrt(norma_b) * np.sqrt(tt))
                matrix_cc[k, l] = cc
        bar0.next()
        fits.writeto(outfolder + '/' + obj1 + '_vg_' + vgindex + '.fits', matrix_cc, overwrite=True)
        with fits.open(outfolder + '/' + obj1 + '_vg_' + vgindex + '.fits', mode='update', verify='ignore') as hcc:
            hcc[0].header['Q0'] = qmin
            hcc[0].header['Q1'] = qmax
            hcc[0].header['DELTA_Q'] = deltaq
            hcc[0].header['T0'] = vector_t[0]
            hcc[0].header['T1'] = vector_t[-1]
            hcc[0].header['DELTA_T'] = vector_t[1] - vector_t[0]
            hcc[0].header['VGAMMA'] = float(vgindex)
            hcc[0].header['WREG'] = '4000-4090,4110-4320,4360-4850,4875-5290,5350-5900'
        for xq in q_array:
            aux1 = str(round(xq, len(str(deltaq+1))-2)).replace('.', '')
            if os.path.isfile('B' + aux1 + '.fits'):
                os.remove('B' + aux1 + '.fits')
            if os.path.isfile('A' + aux1 + '.fits'):
                os.remove('A' + aux1 + '.fits')
        gc.collect()
        plt.close('all')
    
    vexplore(outfolder)
    bar0.finish()
    gc.collect()
    plt.close('all')

def vexplore(folder):
    plt.ioff()
    folder = folder.replace('/', '')
    nfiles = os.listdir(folder)
    vgdicc = {}
    for xf in nfiles:
        if os.path.isfile(os.path.join(folder, xf)) and xf.endswith('.fits'):
            with fits.open(folder + '/' + xf, mode='readonly') as hdul:
                xvg = hdul[0].header['VGAMMA']
                vgdicc[xf] = xvg
    vgd_sort = sorted(vgdicc.items(), key=operator.itemgetter(1), reverse=False)
    n = len(vgd_sort)
    qall = []
    tall = []
    ccall = []
    for i, img in enumerate(vgd_sort):
        with fits.open(folder + '/' + img[0], mode='readonly') as hdul:
            qmin = hdul[0].header['Q0']
            qmax = hdul[0].header['Q1']
            deltaq = hdul[0].header['DELTA_Q']
            num_q = int(np.ceil((qmax - qmin) / deltaq)) + 1
            q_array = np.linspace(qmin, qmin + (num_q - 1) * deltaq, num_q)
            t0 = hdul[0].header['T0']
            t1 = hdul[0].header['T1']
            deltat = hdul[0].header['DELTA_T']
            vector_t = np.arange(t0, t1 + deltat, deltat)
            ccx = hdul[0].data
            qall.append(q_array)
            tall.append(vector_t)
            ccall.append(ccx)
    fig = plt.figure(figsize=[6, 6])
    ax = fig.add_subplot(111, projection='3d')
    ax.set_autoscalez_on(False)
    plt.rc('xtick', labelsize=9)
    plt.rc('ytick', labelsize=9)
    ax.set_zlim3d(bottom=np.min(ccall), top=np.max(ccall))
    ax.set_xlabel("Mass ratio ($q$)", fontsize=9)
    ax.set_ylabel("Temp. [1000 K]", fontsize=9)
    ax.set_title('Systemic RV = ' + str(vgd_sort[0][1]) + ' km/s', y=1)
    X, Y = np.meshgrid(qall[0], tall[0] / 1000)
    ax.plot_surface(X, Y, ccall[0], rstride=1, cstride=1, cmap='jet', vmin=np.min(ccall), vmax=np.max(ccall))
    fig.subplots_adjust(top=1, bottom=0.05, left=0.05, right=0.92, hspace=0.2, wspace=0.02)
    class Index:
        ind = 0
        def next(self, event):
            self.ind += 1
            i = self.ind % n
            ax.cla()
            ax.set_zlim3d(bottom=np.min(ccall), top=np.max(ccall))
            ax.plot_surface(X, Y, ccall[i], rstride=1, cstride=1, cmap='jet', vmin=np.min(ccall), vmax=np.max(ccall))
            ax.set_title('Systemic RV = ' + str(vgd_sort[i][1]) + ' km/s', y=1)
            ax.set_xlabel("Mass ratio ($q$)", fontsize=9)
            ax.set_ylabel("Temp. [1000 K]", fontsize=9)
            plt.draw()
        def prev(self, event):
            self.ind = (self.ind - 1) % n
            i = self.ind % n
            ax.cla()
            ax.set_zlim3d(bottom=np.min(ccall), top=np.max(ccall))
            ax.plot_surface(X, Y, ccall[i], rstride=1, cstride=1, cmap='jet', vmin=np.min(ccall), vmax=np.max(ccall))
            ax.set_title('Systemic RV = ' + str(vgd_sort[i][1]) + ' km/s', y=1)
            ax.set_xlabel("Mass ratio ($q$)", fontsize=9)
            ax.set_ylabel("Temperature [x1000 K]", fontsize=9)
            plt.draw()
    callback = Index()
    axprev = plt.axes([0.05, 0.88, 0.1, 0.075])
    axnext = plt.axes([0.85, 0.88, 0.1, 0.075])
    bnext = Button(axnext, 'Next')
    bnext.on_clicked(callback.next)
    bprev = Button(axprev, 'Previous')
    bprev.on_clicked(callback.prev)
    plt.show()
    plt.close('all')

# ============================================================
# Función qparallel (para ser usada por Pool)
# ============================================================

def qparallel(q_array, lis, spa, spb, deltaq, vgamma):
    import os, sys, time, gc
    import numpy as np
    from astropy.io import fits
    from astropy import units as u
    from specutils import Spectrum
    from specutils.manipulation import SplineInterpolatedResampler
    from PyAstronomy import pyasl
    from astropy.io.fits.verify import VerifyWarning
    import matplotlib
    matplotlib.use('Agg')
    VerifyWarning('ignore')
    
    # Calentamiento de Numba
    try:
        _test_arr = np.array([1.0, 2.0, np.nan, 3.0, 4.0], dtype=np.float64)
        splineclean(_test_arr)
        del _test_arr
    except:
        pass
    
    for idx, xq in enumerate(q_array):
        aux1 = str(round(xq, len(str(deltaq+1))-2)).replace('.', '')
        if not os.path.isfile(spb + aux1 + '.fits'):
            start = time.time()
            try:
                spbina(lis, spa=spa+aux1, spb=spb+aux1, nit=1, q=xq, vgamma=vgamma, showtit=False)
            except Exception as e:
                raise
            elapsed = time.time() - start
            gc.collect()
    gc.collect()
    return None



