import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.animation
import pandas as pd
import argparse
import pdb
from munch import Munch
import tables

import camb
from camb import model, initialpower

from gadget_tools import Snapshot
from field_tools import compute_power_spec
from fitting_fns import halofit


parser = argparse.ArgumentParser(
    description='Density field evolution from Gadget simulation.',
    usage= 'python ./visualize_simulation.py')

parser.add_argument('--simdir', default='/scratch/cprem/sims/', type=str, help='Directory path for all simulations')

parser.add_argument('--simname', type=str, default='bdm_cdm1024', help='Directory name containing the saved data')
parser.add_argument('--cosmo', type=str, default='P18', help='cosmology parameters data from')
parser.add_argument('--rundirs', type=str, default='r1',
                help='Directory name containing the snapshot binaries')

parser.add_argument('--plots_into', type=str, default='/mnt/home/student/cprem/cosmo-sims/sanity_checks/plots/')

parser.add_argument('--snap_i_list', type=str, default=list(range(50,201,50)), help='Snapshot index number')

parser.add_argument('--light_snaps', type=int, default=1, help='save white bg images for pdf notes')

args = parser.parse_args()

grid_size = 512
scheme = 'TSC'
# rundir = 'r1'
interlaced = True

inlcd_str = '_interlaced' if interlaced==True else ''

if args.light_snaps:
    theme = '_light'
    plt.style.use('default')
    # plt.set_cmap('nipy_spectral')
else:
    theme = ''
    plt.style.use('dark_background')
    # plt.set_cmap('nipy_spectral')

# simname = 'bdm_cdm1024' if args.simname is None else args.simname

schemes = ['NGP', 'CIC', 'TSC']
p = schemes.index(scheme) + 1

rundirs = args.rundirs.split(' ')

rundir = rundirs[0]
rundir_str = rundir.replace('/', '_') + '-' + rundirs[-1].split('/')[-1] if len(rundirs)>1 else rundir.replace('/', '_')


plotsdir = os.path.join(args.plots_into, f'{args.simname:s}', f'matter_pk')
os.makedirs(plotsdir, exist_ok=True)

snapdir = os.path.join(args.simdir, args.simname, rundir, 'snaps_sw')

def snapfilen_prefix(snapdirectory, snap_i):
    if os.path.exists(os.path.join(snapdir, f'snapdir')):
        return os.path.join(snapdir, 'snapdir/snapshot_{0:04d}'.format(snap_i))
    else:
        return os.path.join(snapdir, 'snapshot_{0:04d}'.format(snap_i))

def snapfilen(snapdirectory, snap_i):
    snapfilen_prefix_i = snapfilen_prefix(snapdirectory, snap_i)
    if os.path.exists(snapfilen_prefix_i) or os.path.exists(f'{snapfilen_prefix_i:s}.hdf5'):
        return snapfilen_prefix_i
    else:
        return snapfilen_prefix_i + '.0'



# cosmology = 'P18' if args.simname=='bdm_cdm1024' else 'WMAP7'
if args.cosmo =='P18':
    cos_par_vals = (0.3063375, 0.6936625, 0.0484103, 0.6781, 0.9677, 0.815)
elif args.cosmo=='WMAP7':
    cos_par_vals = (0.276, 0.724, 0.045, 0.7, 0.961, 0.811)

cos_pars = Munch()
cos_pars.Om0, cos_pars.Ode0, cos_pars.Ob0, cos_pars.h, cos_pars.ns, cos_pars.sig8 = cos_par_vals
cos_pars.Ombh2 = (cos_pars.Ob0)*cos_pars.h**2
cos_pars.Omch2 = (cos_pars.Om0-cos_pars.Ob0)*cos_pars.h**2

print(args.cosmo, vars(cos_pars))

def Omega(z, Om0):
    E = Om0 * (1+z)**3 + (1-Om0)
    return Om0 * (1+z)**3 / E

def D1(z, Om0):
    Om_m = Omega(z, Om0)
    Om_L = 1 - Om_m
    return 5/2* 1/(1+z) * Om_m / (Om_m**(4/7) - Om_L + (1+Om_m/2)*(1+Om_L/70))

def adjust_lightness(color, amount=0.5):
    # This function from stack for creating lighter or darker versions of matplotlib color
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

def lighter(color): return adjust_lightness(color, 1.5)

def darker(color): return adjust_lightness(color, 0.7)


# i_list = list(range(50,201,50))
# i_list = [0,1]

i_list = [int(x) for x in args.snap_i_list.split(',')]
i_list.sort()

snap = Snapshot(snapfilen(snapdir, 0), snapfrmt='swift')

box_size = snap.box_size
k_nyq = np.pi * grid_size / box_size
k_start = 2* np.pi / box_size

redshifts=[]
for i in i_list:
    snap = Snapshot(snapfilen(snapdir, i), snapfrmt='swift')
    redshifts.append(snap.redshift)


pars_camb = camb.CAMBparams()
pars_camb.set_cosmology(H0=cos_pars.h*100, ombh2=cos_pars.Ombh2, omch2=cos_pars.Omch2)
pars_camb.InitPower.set_params(ns=cos_pars.ns)
#Note non-linear corrections couples to smaller scales than we want
pars_camb.set_matter_power(redshifts=redshifts, kmax=10*k_nyq)

#Linear spectra
pars_camb.NonLinear = model.NonLinear_none
pars_camb.NonLinearModel.set_params(halofit_version='takahashi')
results_camb = camb.get_results(pars_camb)
kh_camb, z_camb, pk_camb = results_camb.get_matter_power_spectrum(minkh=2e-7, maxkh=100, npoints = 5000)
# s8 = np.array(results.get_sigma8())

#Non-Linear spectra (Halofit)
pars_camb.NonLinear = model.NonLinear_both
results_camb.calc_power_spectra(pars_camb)
kh_camb_nonlin, z_camb_nonlin, pk_camb_nonlin = results_camb.get_matter_power_spectrum(minkh=k_start*0.8, maxkh=1.2*k_nyq, npoints = 200)




fig1, ax2 = plt.subplots(1, figsize=(7.5,7), dpi=150)
plt.rcParams['lines.linewidth'] = 1

# ax2.plot([],[], ' ', label=f"Scheme-{scheme}, Grid-size: {grid_size:d}")

for index, i in enumerate(i_list[::-1]):
    color=next(ax2._get_lines.prop_cycler)['color']

    k_full = kh_camb
    pk_lin = pk_camb[index]


    ax2.plot(k_full, pk_lin, color=color, linestyle=(0, (1, 1)), zorder=3)
    ax2.set_xscale('log')
    ax2.set_yscale('log')


    ax2.plot(kh_camb_nonlin, pk_camb_nonlin[index], linestyle='solid', color=darker(color), zorder=1)

    pk_fit = halofit.NonLinPowerSpecCDM(Omega(snap.redshift, snap.Omega_m_0))
    pk_fit.set_Del2L_interpolate(k_full, pk_lin)
    pk_fit.compute_params(param_from='smith')
    print(vars(pk_fit))
    # ax2.plot(kh_camb_nonlin, pk_fit.P(kh_camb_nonlin), linestyle='solid', color=darker(color), zorder=2)

    lin_bin = np.linspace(np.sqrt(k_start*k_nyq),k_nyq, 30)
    log_bin = np.logspace(np.log10(k_start),np.log10(k_nyq), 25)
    # merge_bin = np.concatenate([lin_bin,log_bin])
    # merge_bin.sort()power_spec['k'].iloc[1]
    merge_bin = np.union1d(lin_bin, log_bin)


    # power_spec_grouped1 = power_spec.groupby(pd.cut(power_spec['k'], bins=lin_bin)).mean()
    # ax2.plot(power_spec_grouped1['k'],power_spec_grouped1['Pk'], color=color, linestyle='dashed', label=f"{scheme:s} scheme without interlacing")[0]
    
    # if interlaced:
    
    power_spec_all = None

    for rundir in rundirs:
        snapdir = os.path.join(args.simdir, args.simname, rundir, 'snaps_sw')
        snap = Snapshot(snapfilen(snapdir, i), snapfrmt='swift')
        savesdir = os.path.join(args.simdir, args.simname, rundir)
        print(savesdir)
        dens_griddir = os.path.join(savesdir,'meshgrid')
        Pkdir = os.path.join(savesdir,'power_spectrum')

        os.makedirs(Pkdir, exist_ok=True)

        pk_filepath = os.path.join(Pkdir, f'{scheme:s}_{grid_size:d}{inlcd_str}_{i:03d}.csv' )
        
        h5file = tables.open_file( os.path.join(dens_griddir, f'{scheme:s}_{grid_size:d}_{i:03d}.hdf5'), 'r' )

        if interlaced:
            power_spec_this = compute_power_spec(h5file.root.PartType1.density,snap.box_size, interlace_with_FX=h5file.root.PartType1.density_shifted, Win_correct_scheme=scheme, grid_size=grid_size)

        else:
            power_spec_this = compute_power_spec(h5file.root.PartType1.density,snap.box_size, interlace_with_FX=None, Win_correct_scheme=scheme, grid_size=grid_size)

        power_spec_this_grp = power_spec_this.groupby(pd.cut(power_spec_this['k'], bins=merge_bin)).mean()
        power_spec_this_grp.to_csv(pk_filepath, sep='\t', index=False, float_format='%.8e', header=['k (h/Mpc)', 'P(k) (Mpc/h)^3'])

        h5file.close()

        power_spec_this_grp = pd.read_csv(pk_filepath, sep='\t', dtype='float64', names=['k', 'Pk'], header=0)
        # power_spec_folding_this = pd.read_csv(os.path.join(snapdir,f"powerspecs/powerspec_{i:03d}.txt"), sep='\s+', usecols=[0,1], names=['k', 'Delk'], skiprows=5)

        if power_spec_all is None:
            power_spec_all = power_spec_this_grp.copy()
            # power_spec_folding_all = power_spec_folding_this.copy()
        else:
            power_spec_all.append( power_spec_this_grp )
            # power_spec_folding_all.append( power_spec_folding_this )

        # print(power_spec_folding_this)

    power_spec_grp1 = power_spec_all.groupby(pd.cut(power_spec_all['k'], bins=merge_bin)).mean()

    # power_spec_folding_all.sort_values('k', inplace=True)
    # power_spec_folding = power_spec_folding_all[power_spec_folding_all['k'].between(1e-3,1e2)]
    # power_spec_folding_grp1 = power_spec_folding.groupby(pd.cut(power_spec_folding['k'], bins=merge_bin)).mean()
    # power_spec_folding_grp1['pk'] = power_spec_folding_grp1['Delk']*power_spec_folding_grp1['k']**-3*2*np.pi**2

    ax2.plot(power_spec_grp1['k'],power_spec_grp1['Pk'], color=color, linestyle='dashed', label=f"z={f'{round(snap.redshift,1):f}'.rstrip('0').rstrip('.'):s}")[0]
    ax2.scatter(power_spec_grp1['k'],power_spec_grp1['Pk'], color=color, s=4)

    # power_spec_folding_grp1.plot('k', 'pk', loglog=True, color=lighter(color), linestyle='dashdot', lw=0.8, ax=ax2, label='', legend=False)

ax2.plot([],[], ' ', label=f"GADGET-4 simulation")
ax2.plot([],[], linestyle='dashed', color='gray', label=f"  {scheme}-{grid_size:d} grid")
ax2.plot([],[], linestyle='dashdot', color=lighter('gray'), label=f"  folding technique")
# ax2.plot([],[], linestyle='dashed', color=lighter(lighter('gray')), label=f"  for reference")
ax2.plot([],[], ' ', label=f"Halofit model")
ax2.plot([],[], linestyle='solid', color=darker('gray'), label='  Takahashi, 2012')
# ax2.plot([],[], linestyle='solid', color=lighter('gray'), label='  HMcode-2020')
# ax2.plot([],[], ' ', label=f"linear theory")
ax2.plot([],[], linestyle=(0, (1, 1)), color='gray', label='linear theory')
# in snapshot-{i:03d}
# pdb.set_trace()


# ax2[1].plot(power_spec['lam'],power_spec['Pk'])
ax2.set_xlabel(r"$k$  ($h$ Mpc$^{-1}$)")
ax2.set_ylabel(r"$P(k)$  ($h^{-1}$Mpc)$^3$")
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(k_start*0.8,1.2*k_nyq)
ax2.set_ylim(top=1e5, bottom=2e-1)
ax2.grid(True)
ax2.set_title(f"Matter power spectrum at different redshifts")
ax2.legend()

plt.tight_layout()

fig1.savefig(os.path.join(plotsdir, f'{rundir_str:s}_{scheme}{theme:s}.pdf'))
# fig1.savefig(os.path.join(plotsdir, f'{rundir_str:s}_{scheme}{theme:s}.png'))
# fig1.savefig(os.path.join(plotsdir, f'single_snapshot_pk_{i:03d}{theme:s}.svg'))