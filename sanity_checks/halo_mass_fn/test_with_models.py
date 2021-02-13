import numpy as np
import h5py
import tables
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 1
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath,xfrac}']
# mpl.rcParams['text.latex.unicode'] = True
# mpl.rcParams['text.usetex'] = True


import argparse, os, sys

import hmf
from hmf import MassFunction

# from hmf.halos import mass_definitions as md

from astropy.cosmology import Planck18
from scipy.interpolate import interp1d

import vrpy_tools

# plt.style.use('dark_background')
# plt.style.use('default')
# plt.rc('text', usetex=True)
# plt.rc(r'text.latex', preamble=r'\usepackage{xfrac}')

parser = argparse.ArgumentParser()

parser.add_argument('--simnm', type=str)
parser.add_argument('--rund', type=str)
parser.add_argument('--i', type=int)

args = parser.parse_args()

L = int(args.simnm.split('_')[0][1:])
# i = args.i
rund = args.rund

dir_simdata_base = "/scratch/cprem/sims/"

dir_simdata_allruns = os.path.join(dir_simdata_base, args.simnm)

dir_simdata = os.path.join(dir_simdata_allruns, rund)
dir_simdata_vr = os.path.join(dir_simdata, 'halos_vr_6d')
dir_simdata_rs = os.path.join(dir_simdata, 'halos_rs')

siminfo = vrpy_tools.ReadSimInfo(os.path.join(dir_simdata_vr, f'out_{args.i:03d}'))

z = 1/siminfo['ScaleFactor'] - 1

tnk_inp_fn = {}
for model, massdef in zip(["SOVirial", 'SOMean', 'SOCritical'], ['mvir', 'm200m', 'm200c']):
    hal_mass_fn = MassFunction()
    hal_mass_fn.update(cosmo_model=Planck18)
    hal_mass_fn.update(z=z)
    hal_mass_fn.update(hmf_model="Tinker08")
    hal_mass_fn.update(Mmin  = 8, Mmax = 15)
    hal_mass_fn.update(mdef_model  = model,)
    tnk_inp_fn[massdef] = interp1d(np.log10(hal_mass_fn.m), hal_mass_fn.dndlog10m)

# hal_mass_fn.update(mdef_model  = model, mdef_params = {"overdensity":200,})

bw=0.04
bins_edges = np.arange(8,15,bw)
bins_cen = (bins_edges[:-1] + bins_edges[1:]) / 2

numbins = bins_cen.shape[0]

hist_rsvir = np.zeros(numbins)
hist_rs200c = np.zeros(numbins)
hist_rs200m = np.zeros(numbins)

hist_vr200c = np.zeros(numbins)
hist_vr200m = np.zeros(numbins)
hist_vrvir = np.zeros(numbins)
hist_vrfof = np.zeros(numbins)

for rund in ('r4', 'r5', 'r6', 'r7'):
    dir_simdata = os.path.join(dir_simdata_allruns, rund)
    dir_simdata_vr = os.path.join(dir_simdata, 'halos_vr_6d')
    dir_simdata_rs = os.path.join(dir_simdata, 'halos_rs')

    hal_rs = pd.read_csv(os.path.join(dir_simdata_rs, f'out_wp_{args.i:d}.list'), sep=r'\s+', header=0, skiprows=list(range(1,16)), engine='c')
    hal_rs = hal_rs[hal_rs['PID']==-1]

    hal_vr = tables.open_file(os.path.join(dir_simdata_vr, f'out_{args.i:03d}.properties.0'), 'r')
    select_hal = np.where(hal_vr.root.hostHaloID[:]==-1)

    num_hal_vr = select_hal[0].shape[0]
    num_hal_rs = hal_rs['#ID'].shape[0]

    hist_rsvir += np.histogram(np.log10(hal_rs['Mvir']), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, )[0]
    hist_rs200c += np.histogram(np.log10(hal_rs['M200c']), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, )[0]
    hist_rs200m += np.histogram(np.log10(hal_rs['M200b']), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, )[0]

    hist_vr200c += np.histogram(np.log10(hal_vr.root.Mass_200crit[select_hal]) + 10-0.1, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]
    hist_vr200m += np.histogram(np.log10(hal_vr.root.Mass_200mean[select_hal]) + 10-0.1, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]
    hist_vrvir += np.histogram(np.log10(hal_vr.root.Mvir[select_hal]) + 10, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]
    hist_vrfof += np.histogram(np.log10(hal_vr.root.Mass_FOF[select_hal]) + 10, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]


    hal_vr.close()  
    del hal_rs, hal_vr

hist_rsvir /= 4
hist_rs200c /= 4
hist_rs200m /= 4

hist_vr200c /= 4
hist_vr200m /= 4
hist_vrvir /= 4
hist_vrfof /= 4


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

hist_rsvir = smooth(hist_rsvir, 3)
hist_rs200c = smooth(hist_rs200c, 3)
hist_rs200m = smooth(hist_rs200m, 3)

hist_vr200c = smooth(hist_vr200c, 3)
hist_vr200m = smooth(hist_vr200m, 3)
hist_vrvir = smooth(hist_vrvir, 3)
hist_vrfof = smooth(hist_vrfof, 3)




fig, (ax1,ax2) = plt.subplots(2, figsize=(14,9), dpi=250, sharex=True)
plt.subplots_adjust(hspace=.1)

color=next(ax1._get_lines.prop_cycler)['color']

ax1.plot(bins_cen, hist_rsvir, label="RS_Mvir")
ax1.plot(bins_cen, hist_rs200c, label="RS M200c")
ax1.plot(bins_cen, hist_rs200m, label="RS M200m")
ax1.plot(bins_cen, hist_vrfof, label="VR MFOF")
ax1.plot(bins_cen, hist_vrvir, label="VR Mvir")
ax1.plot(bins_cen, hist_vr200c, label="VR M200c")
ax1.plot(bins_cen, hist_vr200m, label="VR M200m")

 
# rsvir = ax1.hist(np.log10(hal_rs['Mvir']), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, histtype='step', log=True,label="RS_Mvir")
# rs200c = ax1.hist(np.log10(hal_rs['M200c']), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, histtype='step', log=True,label="RS_M200c")

# vr200c= ax1.hist(np.log10(hal_vr.root.Mass_200crit[select_hal]) + 10-0.1, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3, histtype='step', log=True, label="VR_M200c")
# vrvir= ax1.hist(np.log10(hal_vr.root.Mvir[select_hal]) + 10, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3, histtype='step', log=True, label="VR_Mvir")
# vrfof= ax1.hist(np.log10(hal_vr.root.Mass_FOF[select_hal]) + 10, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3, histtype='step', log=True, label="VR_MFOF")

# tinkerplot = ax1.plot(np.log10(hal_mass_fn.m), hal_mass_fn.dndlog10m, label="Tinker-08", color='black')

# ax1.plot(bins_cen, tinker_m200c(bins_cen), label="Tinker-08", color='black')
ax1.plot(bins_cen, tnk_inp_fn['mvir'](bins_cen), label="Tinker Mvir")
ax1.plot(bins_cen, tnk_inp_fn['m200c'](bins_cen), label="Tinker M200c")
ax1.plot(bins_cen, tnk_inp_fn['m200m'](bins_cen), label="Tinker  M200m")

ax2.set_xlabel('log(M) where mass is in $h^{-1}~M_{\odot}$')
ax1.set_ylabel(r'$\frac{dn}{d (\log ~M)}$ in $h^{3}Mpc^{-3}$')
ax1.set_title(f'Halo mass function at redshift z = {z:.1g}')
ax1.set_xlim(8, 15)
# ax1.set_ylim(1e-6,6e-2)
ax1.set_yscale('log')
ax1.legend()
ax1.grid()

ax2.plot(bins_cen, tnk_inp_fn['mvir'](bins_cen)/hist_vrfof, label="Tinker/ VR MFOF")
ax2.plot(bins_cen, tnk_inp_fn['m200m'](bins_cen)/hist_vrvir, label="Tinker/ VR Mvir")
ax2.plot(bins_cen, tnk_inp_fn['m200c'](bins_cen)/hist_vr200c, label="Tinker/ VR M200c")
ax2.plot(bins_cen, tnk_inp_fn['m200m'](bins_cen)/hist_vr200m, label="Tinker/ VR M200m")
ax2.plot(bins_cen, tnk_inp_fn['m200c'](bins_cen)/hist_rs200c, label="Tinker/ RS M200c")
ax2.set_ylim(0.5,1.4)
# ax2.set_yscale('log')
ax2.legend()
ax2.grid()

 
fig.savefig("plots/hmf.png")
