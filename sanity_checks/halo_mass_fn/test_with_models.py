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

from astropy.cosmology import FlatLambdaCDM, Planck18
from scipy.interpolate import interp1d

import vrpy_tools

# plt.style.use('dark_background')
# plt.style.use('default')
# plt.rc('text', usetex=True)
# plt.rc(r'text.latex', preamble=r'\usepackage{xfrac}')

parser = argparse.ArgumentParser()

parser.add_argument('--simnm', type=str)
parser.add_argument('--runds', type=str)
parser.add_argument('--i', type=int)
parser.add_argument('--suffix', type=str)

args = parser.parse_args()

# L = 200
L = int(args.simnm.split('_')[0][1:])
# i = args.i
runds = args.runds.split(' ')
rund = runds[0]

dir_simdata_base = "/scratch/cprem/sims/"

dir_simdata_allruns = os.path.join(dir_simdata_base, args.simnm)

dir_simdata = os.path.join(dir_simdata_allruns, rund)

halos_vr_dirname = 'halos_vr_6d'
if args.suffix is not None: halos_vr_dirname += '_'+ args.suffix

dir_simdata_vr = os.path.join(dir_simdata, halos_vr_dirname)
dir_simdata_rs = os.path.join(dir_simdata, 'halos_rs')

siminfo = vrpy_tools.ReadSimInfo(os.path.join(dir_simdata_vr, f'out_{args.i:03d}'))

z = 1/siminfo['ScaleFactor'] - 1
Om_m = siminfo['Omega_m']

p18py = Planck18.clone(name='Planck18 modified', H0=siminfo['h_val']*siminfo['Hubble_unit'], Om0=siminfo['Omega_m'], Ob0=0.0484103)
# p18py = Planck18
print(p18py)

tnk_inp_fn = {}
for model, massdef, overdensity in zip(["SOVirial", 'SOMean', 'SOCritical', 'SOCritical'], ['mvir', 'm200m', 'm200c', 'm500c'], [None, 200, 200, 500]):
    hal_mass_fn = MassFunction()
    hal_mass_fn.update(cosmo_model=p18py)
    hal_mass_fn.update(z=z)
    hal_mass_fn.update(hmf_model="Tinker08")
    hal_mass_fn.update(Mmin  = 8, Mmax = 15)
    if overdensity is None:
        hal_mass_fn.update(mdef_model  = model)
    else:
        # overdensity *= p18py.Om(z)
        hal_mass_fn.update(mdef_model  = model, mdef_params = {"overdensity": overdensity})
    tnk_inp_fn[massdef] = interp1d(np.log10(hal_mass_fn.m), hal_mass_fn.dndlog10m)

# hal_mass_fn.update(mdef_model  = model, mdef_params = {"overdensity":200,})

bw=0.04
bins_edges = np.arange(8,15,bw)
bins_cen = (bins_edges[:-1] + bins_edges[1:]) / 2

numbins = bins_cen.shape[0]

hist_rsvir = np.zeros(numbins)
hist_rs200m = np.zeros(numbins)
hist_rs200c = np.zeros(numbins)
hist_rs500c = np.zeros(numbins)

hist_vrvir = np.zeros(numbins)
hist_vr200m = np.zeros(numbins)
hist_vr200c = np.zeros(numbins)
hist_vr500c = np.zeros(numbins)
hist_vrfof = np.zeros(numbins)

# runds = ('r4', 'r5', 'r6', 'r7')

for rund in runds:
    dir_simdata = os.path.join(dir_simdata_allruns, rund)
    dir_simdata_vr = os.path.join(dir_simdata, halos_vr_dirname)
    dir_simdata_rs = os.path.join(dir_simdata, 'halos_rs')
    # dir_simdata_rs = '/scratch/cprem/sims/L200_N512_Cp18/r1/halos_rs/'

    hal_rs_all = pd.read_csv(os.path.join(dir_simdata_rs, f'out_wp_{args.i:d}.list'), sep=r'\s+', header=0, skiprows=list(range(1,16)), engine='c')
    hal_rs_all = hal_rs_all.loc[hal_rs_all.PID==-1]
    hal_rs_all['kin_rat'] = 2*hal_rs_all['T/|U|']

    # hal_vr = tables.open_file(os.path.join(dir_simdata_vr, f'out_{args.i:03d}.properties.0'), 'r')
    hal_vr_dict, N_hal_tot_vr = vrpy_tools.ReadPropertyFile(os.path.join(dir_simdata_vr, f'out_{args.i:03d}'), ibinary=2, isiminfo=False, iunitinfo=False, iconfiginfo=False)
    hal_vr_all = pd.DataFrame.from_dict(hal_vr_dict)
    hal_vr_all = hal_vr_all.loc[hal_vr_all.hostHaloID==-1]
    hal_vr_all['kin_rat'] = 2*hal_vr_all.Ekin/np.abs(hal_vr_all.Epot)
    

    hal_rs = hal_rs_all.loc[hal_rs_all.kin_rat.between(0.5,1.8)]
    hal_vr = hal_vr_all.loc[hal_vr_all.kin_rat.between(0.5,1.8)]

    num_hal_vr = hal_vr.shape[0]
    num_hal_rs = hal_rs.shape[0]

    print(f"Number of host halos found by ROCKSTAR is {num_hal_rs}.")
    print(f"Number of host halos found by VELOCIraptor is {num_hal_vr}.")

    offset_Mfactor = 0 #.1687
    hist_rsvir += np.histogram(np.log10(hal_rs.Mvir), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, )[0]
    hist_rs200m += np.histogram(np.log10(hal_rs.M200b), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, )[0]
    hist_rs200c += np.histogram(np.log10(hal_rs.M200c), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, )[0]
    hist_rs500c += np.histogram(np.log10(hal_rs.M500c), bins=bins_edges, weights=1*np.ones(num_hal_rs)/bw/L**3, )[0]

    hist_vrvir += np.histogram(np.log10(hal_vr.Mass_BN98)+offset_Mfactor, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]
    hist_vr200m += np.histogram(np.log10(hal_vr.Mass_200mean)+offset_Mfactor, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]
    hist_vr200c += np.histogram(np.log10(hal_vr.Mass_200crit)+offset_Mfactor, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]
    hist_vr500c += np.histogram(np.log10(hal_vr.SO_Mass_500_rhocrit)+offset_Mfactor, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]
    hist_vrfof += np.histogram(np.log10(hal_vr.Mass_FOF)+offset_Mfactor, bins=bins_edges, weights=1*np.ones(num_hal_vr)/bw/L**3)[0]
    
    hal_vr.close() 
    del hal_rs, hal_vr

num_runs = len(runds)

hist_rsvir /= num_runs
hist_rs200m /= num_runs
hist_rs200c /= num_runs
hist_rs500c /= num_runs

hist_vrvir /= num_runs
hist_vr200m /= num_runs
hist_vr200c /= num_runs
hist_vr500c /= num_runs
hist_vrfof /= num_runs


def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

hist_rsvir = smooth(hist_rsvir, 4)
hist_rs200m = smooth(hist_rs200m, 4)
hist_rs200c = smooth(hist_rs200c, 4)
hist_rs500c = smooth(hist_rs500c, 4)

hist_vrvir = smooth(hist_vrvir, 4)
hist_vr200m = smooth(hist_vr200m, 4)
hist_vr200c = smooth(hist_vr200c, 4)
hist_vr500c = smooth(hist_vr500c, 4)
hist_vrfof = smooth(hist_vrfof, 4)




colors = ['blue', 'green', 'brown', 'purple', 'orange' ]
mpl.rcParams['lines.linewidth'] = 1

fig, (ax1,ax2) = plt.subplots(2, figsize=(10,7), dpi=200, sharex=True)
plt.subplots_adjust(hspace=.1)

# color=next(ax1._get_lines.prop_cycler)['color']

mpl.rcParams['lines.linestyle'] = 'dotted'
pltrsmvir, = ax1.plot(bins_cen, hist_rsvir, color=colors[0])
pltrsm200m, = ax1.plot(bins_cen, hist_rs200m, color=colors[1])
pltrsm200c, = ax1.plot(bins_cen, hist_rs200c, color=colors[2])
pltrsm500c, = ax1.plot(bins_cen, hist_rs500c, color=colors[3])

mpl.rcParams['lines.linestyle'] = '--'

pltvrmvir, = ax1.plot(bins_cen, hist_vrvir, color=colors[0])
pltvrm200m, = ax1.plot(bins_cen, hist_vr200m, color=colors[1])
pltvrm200c, = ax1.plot(bins_cen, hist_vr200c, color=colors[2])
pltvrm500c, = ax1.plot(bins_cen, hist_vr500c, color=colors[3])
pltvrmfof, = ax1.plot(bins_cen, hist_vrfof, color=colors[4])


mpl.rcParams['lines.linestyle'] = 'solid'
# ax1.plot(bins_cen, tinker_m200c(bins_cen), label="Tinker-08", color='black')
ax1.plot(bins_cen, tnk_inp_fn['mvir'](bins_cen), ls='-', color=colors[0])
ax1.plot(bins_cen, tnk_inp_fn['m200m'](bins_cen), ls='-', color=colors[1])
plttkm200c, = ax1.plot(bins_cen, tnk_inp_fn['m200c'](bins_cen), ls='-', color=colors[2])
ax1.plot(bins_cen, tnk_inp_fn['m500c'](bins_cen), ls='-', color=colors[3])

ax2.set_xlabel(r'log(M) where mass is in $h^{-1}~M_{\odot}$')
ax1.set_ylabel(r'$\frac{dn}{d (\log ~M)}$ in $h^{3}Mpc^{-3}$')
ax1.set_title(f'Halo mass function at redshift z = {z:.3g}')
ax1.set_xlim(10, 15)
ax1.set_ylim(1e-6,1e0)
ax1.set_yscale('log')

lines = ax1.get_lines()
legend1 = ax1.legend([pltrsm200c, pltvrm200c, plttkm200c], ["Rockstar", "Velociraptor", "Tinker 2008"], loc=1)
legend2 = ax1.legend([pltvrmvir, pltvrm200m, pltvrm200c, pltvrm500c, pltvrmfof], [r'$M_{\rm{vir}}$', r'$M_{\rm{200m}}$', r'$M_{\rm{200c}}$', r'$M_{\rm{500c}}$', r'$M_{\rm{FOF}}$'], loc=3)
ax1.add_artist(legend1)
ax1.add_artist(legend2)

# ax1.legend()
ax1.grid()

mpl.rcParams['lines.linestyle'] = 'dotted'
ax2.plot(bins_cen, hist_rsvir/tnk_inp_fn['mvir'](bins_cen), color=colors[0])
ax2.plot(bins_cen, hist_rs200m/tnk_inp_fn['m200m'](bins_cen), color=colors[1])
ax2.plot(bins_cen, hist_rs200c/tnk_inp_fn['m200c'](bins_cen), color=colors[2])
ax2.plot(bins_cen, hist_rs500c/tnk_inp_fn['m500c'](bins_cen), color=colors[3])
ax2.plot([],[], color=colors[0], label='Rockstar/Tinker')

mpl.rcParams['lines.linestyle'] = (0,(5,1,3,1))
ax2.plot(bins_cen, hist_vrvir/tnk_inp_fn['mvir'](bins_cen), color=colors[0])
ax2.plot(bins_cen, hist_vr200m/tnk_inp_fn['m200m'](bins_cen), color=colors[1])
ax2.plot(bins_cen, hist_vr200c/tnk_inp_fn['m200c'](bins_cen), color=colors[2])
ax2.plot(bins_cen, hist_vr500c/tnk_inp_fn['m500c'](bins_cen), color=colors[3])
# ax2.plot(bins_cen, hist_vrfof/tnk_inp_fn['mvir'](bins_cen), color=colors[4])
ax2.plot([],[], color=colors[0], label='Velociraptor/Tinker')

ax2.set_ylim(0.55,1.6)
# ax2.set_yscale('log')

mpl.rcParams['lines.linestyle'] = 'solid'

ax2.legend()

major_ticks = np.arange(0.6, 1.41, 0.2)
minor_ticks = np.arange(0.5, 1.51, 0.1)

ax2.set_yticks(major_ticks)
ax2.set_yticks(minor_ticks, minor=True)

ax2.grid(b=True, which='major')
ax2.grid(b=True, which='minor', ls='dashdot')

# plt.show()

fig.tight_layout()
plot_path = f"plots/{args.simnm}/hmf/"
os.makedirs(plot_path, exist_ok = True)
plot_suffix = f'-{args.suffix:s}' if args.suffix is not None else ''
fig.savefig(f"{plot_path}{''.join(runds)}_{args.i:03d}{plot_suffix:s}.pdf")
# fig.savefig(f"{plot_path}{''.join(runds)}_{args.i:03d}{plot_suffix:s}.png", dpi=200)