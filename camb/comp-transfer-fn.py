#%%
import numpy as np
import pandas as pd
import os
import argparse

import camb
from camb import model, initialpower
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))


from munch import Munch


parser = argparse.ArgumentParser(description='Give cosmological parameters')
parser.add_argument('--params', type=str)
parser.add_argument('--cambfile', type=str)
args = parser.parse_args()

if args.params is None:
    params_list = [0.306, 0.694, 0.0484, 0.678, 0.9677, 0.815]
else:
    params_list = [float(x) for x in args.params.split(',')]

#%%
# params_list = [0.306, 0.694, 0.0484, 0.678, 0.9677, 0.815]
cos_pars = Munch()
cos_pars.Om0, cos_pars.Ode0, cos_pars.Ob0, cos_pars.h, cos_pars.ns, cos_pars.sig8 = params_list

cos_pars.Ombh2 = (cos_pars.Ob0)*cos_pars.h**2
cos_pars.Omch2 = (cos_pars.Om0-cos_pars.Ob0)*cos_pars.h**2

#Now get matter power spectra and sigma8 at redshift 0 and 0.8
pars = camb.CAMBparams()
pars.set_cosmology(H0=cos_pars.h*100, ombh2=cos_pars.Ombh2, omch2=cos_pars.Omch2)
pars.InitPower.set_params(ns=cos_pars.ns)
#Note non-linear corrections couples to smaller scales than you want
pars.set_matter_power(redshifts=[0],  kmax=100.0, k_per_logint=100)

#Non-Linear spectra (Halofit)
pars.NonLinear = model.NonLinear_both
pars.NonLinearModel.set_params(halofit_version='takahashi')
# pars.halofit_version='mead2020'
results = camb.get_transfer_functions(pars)

trans = results.get_matter_transfer_data()
trans_df = pd.DataFrame(data=trans.transfer_data[:,:,-1].T)

trans_df.to_csv(f'{args.cambfile:s}', sep='\t', index=False, header=False)