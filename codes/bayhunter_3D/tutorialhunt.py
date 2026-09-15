# #############################
#
# Copyright (C) 2018
# Jennifer Dreiling   (dreiling@gfz-potsdam.de)
#
#
# #############################

import os
# set os.environment variables to ensure that numerical computations
# do not do multiprocessing !! Essential !! Do not change !
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
import os.path as op
import matplotlib
matplotlib.use('PDF')

from BayHunter import PlotFromStorage
from BayHunter import Targets
from BayHunter import utils
from BayHunter import MCMC_Optimizer
from BayHunter import ModelMatrix
from BayHunter import SynthObs
import logging
from scipy.stats import binned_statistic
import subprocess
import os
import warnings
from glob import glob

import matplotlib.pyplot as plt

path = '/media/ivan/Seagate Hub/COPY_IVAN/NOISETOMO/Projects/Poas/observed'
max_depth = 7
#
# console printout formatting
#
formatter = ' %(processName)-12s: %(levelname)-8s |  %(message)s'
logging.basicConfig(format=formatter, level=logging.INFO)
logger = logging.getLogger()

indxs = np.arange(0, len(glob(path + '/disp_2Dmap_*.dat')), 3, dtype=int)

for ix in indxs:
    try:
        disp_file = sorted(glob(path + '/disp_2Dmap_*.dat'))[ix]
        if os.path.exists("results_{:s}".format(os.path.basename(disp_file))) == True:
           print(disp_file+ 'Already exists')
           continue
        else:
           print(disp_file)
        #
        # ------------------------------------------------------------  obs SYNTH DATA
        #
        # Load priors and initparams from config.ini or simply create dictionaries.
           initfile = 'config.ini'
           priors, initparams = utils.load_params(initfile)
           initparams["savepath"] += "{:s}".format(os.path.basename(disp_file))
           # Load observed data (synthetic test data)
           xsw, _ysw = np.loadtxt(disp_file).T
           #xrf, _yrf = np.loadtxt('observed/st3_prf.dat').T

           # add noise to create observed data
           # order of noise values (correlation, amplitude):
           # noise = [corr1, sigma1, corr2, sigma2] for 2 targets
           noise = [0.0, 0.012, 0.98, 0.005]
           ysw_err = SynthObs.compute_expnoise(_ysw, corr=noise[0], sigma=noise[1])
           ysw = _ysw + ysw_err
           #yrf_err = SynthObs.compute_gaussnoise(_yrf, corr=noise[2], sigma=noise[3])
           #yrf = _yrf + yrf_err


        #
        # -------------------------------------------  get reference model for BayWatch
        #
        # Create truemodel only if you wish to have reference values in plots
        # and BayWatch. You ONLY need to assign the values in truemodel that you
        # wish to have visible.
           dep, vs = np.loadtxt('observed/st3_mod.dat', usecols=[0, 2], skiprows=1).T
           pdep = np.concatenate((np.repeat(dep, 2)[1:], [150]))
           pvs = np.repeat(vs, 2)

           truenoise = np.concatenate(([noise[0]], [np.std(ysw_err)]))  # target 2

           explike = SynthObs.compute_explike(yobss=[ysw], ymods=[_ysw],
                                           noise=truenoise, gauss=[False, True],
                                           rcond=initparams['rcond'])
           truemodel = {'model': (pdep, pvs),
                     'nlays': 20,
                     'noise': truenoise,
                     'explike': explike,
                     }

           print(truenoise, explike)


        #
        #  -----------------------------------------------------------  DEFINE TARGETS
        #
        # Only pass x and y observed data to the Targets object which is matching
        # the data type. You can chose for SWD any combination of Rayleigh, Love, group
        # and phase velocity. Default is the fundamendal mode, but this can be updated.
        # For RF chose P or S. You can also use user defined targets or replace the
        # forward modeling plugin wih your own module.
           target1 = Targets.RayleighDispersionGroup(xsw, ysw, yerr=ysw_err)
        #target2 = Targets.PReceiverFunction(xrf, yrf)
        #target2.moddata.plugin.set_modelparams(gauss=1., water=0.01, p=6.4)

        # Join the targets. targets must be a list instance with all targets
        # you want to use for MCMC Bayesian inversion.
           targets = Targets.JointTarget(targets=[target1])


        #
        #  ---------------------------------------------------  Quick parameter update
        #
        # "priors" and "initparams" from config.ini are python dictionaries. You could
        # also simply define the dictionaries directly in the script, if you don't want
        # to use a config.ini file. Or update the dictionaries as follows, e.g. if you
        # have station specific values, etc.
        # See docs/bayhunter.pdf for explanation of parameters

           priors.update({#'mohoest': (38, 4),  # optional, moho estimate (mean, std)
                       'rfnoise_corr': 0.98,
                       'swdnoise_corr': 0.
                       # 'rfnoise_sigma': np.std(yrf_err),  # fixed to true value
                       # 'swdnoise_sigma': np.std(ysw_err),  # fixed to true value
                       })

           initparams.update({'nchains': 20,
                           'iter_burnin': (2048 * 100),
                           'iter_main': (2048 * 50),
                           'propdist': (0.015, 0.015, 0.015, 0.005, 0.005),
                           # maxmodels = per-chain storage array size for accepted
                           # p2 (main-phase) models. Must comfortably exceed the
                           # number of iterations, or chains can overflow it
                           # (IndexError in append_currentmodel) once iter_burnin/
                           # iter_main are increased. Sized here with margin.
                           'maxmodels': (2048 * 100) + (2048 * 50),
                           })


        #
        #  -------------------------------------------------------  MCMC BAY INVERSION
        #
        # Save configfile for baywatch. refmodel must not be defined.
           utils.save_baywatch_config(targets, path='.', priors=priors,
                                   initparams=initparams, refmodel=truemodel)
           optimizer = MCMC_Optimizer(targets, initparams=initparams, priors=priors,
                                   random_seed=None)
        # default for the number of threads is the amount of cpus == one chain per cpu.
        # if baywatch is True, inversion data is continuously send out (dtsend)
        # to be received by BayWatch (see below).
           optimizer.mp_inversion(nthreads=6, baywatch=True, dtsend=1)


        #
        # #  ---------------------------------------------- Model resaving and plotting
           savepath = initparams['savepath']
           cfile = '%s_config.pkl' % initparams['station']
           configfile = op.join(savepath, 'data', cfile)
           obj = PlotFromStorage(configfile)
        # The final distributions will be saved with save_final_distribution.
        # Beforehand, outlier chains will be detected and excluded.
        # Outlier chains are defined as chains with a likelihood deviation
        # of dev * 100 % from the median posterior likelihood of the best chain.
           obj.save_final_distribution(maxmodels=1000000, dev=0.1)
        # Save a selection of important plots
        # nchains here only controls how many individual-chain lines are drawn
        # in the diagnostic plots (iiter_*, currentmodels, currentdatafits).
        # It does NOT affect the final posterior (mean/median/mode), which
        # already combines all 30 chains via save_final_distribution above.
        # Capping it avoids loading/plotting all 30 chains' full iteration
        # histories at once, which was likely causing the OOM kill.
           plot_nchains = min(15, initparams['nchains'])
           obj.save_plots(depint=0.1, nchains=plot_nchains)
           obj.merge_pdfs()

        #
        # If you are only interested on the mean posterior velocity model, type:
           file = op.join(initparams['savepath'], 'data/c_models.npy')
           models, = obj._get_posterior_data(["models"], final=True)

           # malla inicial de profundidades para la interpolacion de los modelos
           # AJUSTA el paso (interp0) a la resolucion que quieras usar, en km
           interp0 = 0.25
           dep_int = np.arange(0, max_depth + interp0, interp0)

           singlemodels = ModelMatrix.get_singlemodels(models, dep_int=dep_int)

           for mod_type in ["mean", "median", "mode"]:  # "minmax", "stdminmax", "mode"]:
               vs, dep = singlemodels[mod_type]
               fname = os.path.join(savepath, f"model_{mod_type}.txt")
               np.savetxt(fname, np.column_stack((dep, vs)), delimiter=',', header='Depth(km),Vs(km/s)')

           # get histogram of interface depths
           models2 = ModelMatrix._replace_zvnoi_h(models)
           models2 = np.array([model[~np.isnan(model)] for model in models2], dtype='object')
           yinterf = np.array([np.cumsum(model[int(model.size / 4):-1])
                               for model in models2], dtype='object')
           yinterf = np.concatenate(yinterf)
           maxdepth = int(np.ceil(dep_int.max()))
           interp = dep_int[1] - dep_int[0]
           dep_int = np.arange(dep_int[0], dep_int[-1] + interp / 2., interp / 2.)
           depbins = np.arange(0, maxdepth + 2 * interp, interp)
           stats = binned_statistic(yinterf.astype(float), [], statistic='count', bins=depbins)
           fname = os.path.join(savepath, f"model_interfaces.txt")
           # np.save(finalmod_filename,finalmod)
           np.savetxt(fname, np.column_stack((dep, stats.statistic[:-1])), header="depth interface_count",
                      fmt="%5.1f %5.3f")

           vs_mean, dep = singlemodels['mean']
           vs_median, dep = singlemodels['median']
           vs_mode, depi = singlemodels['mode']
           vs_mode = np.interp(dep, depi, vs_mode)
           vs_std, dep = singlemodels['stdminmax']
           header = "depth,vs_mean,vs_median,vs_mode,vs_std_low,vs_std_high,interface_count"
           try:
               finalmod = np.column_stack(
                   (dep, vs_mean, vs_median, vs_mode, vs_std[0], vs_std[1], stats.statistic[:-1]))
           except:
               finalmod = np.column_stack((dep, vs_mean, vs_median, vs_mode, vs_std[0], vs_std[1], stats.statistic))
           fname = os.path.join(path, f"models_{os.path.basename(disp_file)}.txt")
           np.savetxt(fname, finalmod, delimiter=',', header=header)
           subprocess.run("rm -r ./*/data/c0*.npy", shell=True, check=True)

           #np.savetxt("./results_{:s}/model.txt".format(os.path.basename(disp_file)), np.column_stack((dep, vs)), delimiter=',', header='Depth(km),Vs(km/s)', comments='')
           print(f"Saved models to {fname}")

    except Exception as e:
        import traceback
        print(f"[ERROR] fallo procesando indice {ix}: {e}")
        traceback.print_exc()
        continue
