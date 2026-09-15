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
import subprocess
import os
import warnings
from glob import glob
import sys
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
# Slurm environment variables: used to create distinct savepaths
jobID = os.environ.get("SLURM_ARRAY_JOB_ID")
arrayID = os.environ.get("SLURM_ARRAY_TASK_ID")
nthreads = int(os.environ.get("SLURM_CPUS_PER_TASK"))

path = 'observed'
max_depth = 2.0
#

# console printout formatting
#
formatter = ' %(processName)-12s: %(levelname)-8s |  %(message)s'
logging.basicConfig(format=formatter, level=logging.INFO)
logger = logging.getLogger()

disp_file = sys.argv[1]
#grid_linear_index = os.path.basename(disp_file).split("_")[3] + "_" os.path.basename(disp_file).split("_")[4]

try:
    if os.path.exists("results_{:s}".format(os.path.basename(disp_file))) == True:
       print(disp_file+ 'Already exists')       
       sys.exit()
    else:
       print(disp_file)
    #
    # ------------------------------------------------------------  obs SYNTH DATA
    #
    # Load priors and initparams from config.ini or simply create dictionaries.
       initfile = 'config.ini'
       priors, initparams = utils.load_params(initfile)
       initparams["savepath"] += f"_{jobID}-{arrayID}"

       # Load observed data (synthetic test data)
       xsw, _ysw = np.loadtxt(disp_file).T
       xsw_ph, _ysw_ph = np.loadtxt(disp_file.replace('group', 'phase')).T
       #xrf, _yrf = np.loadtxt('observed/st3_prf.dat').T

       # add noise to create observed data
       # order of noise values (correlation, amplitude):
       # noise = [corr1, sigma1, corr2, sigma2] for 2 targets
       noise = [0.0, 0.012, 0.98, 0.005]
       ysw_err = SynthObs.compute_expnoise(_ysw, corr=noise[0], sigma=noise[1])
       ysw_err_ph = SynthObs.compute_expnoise(_ysw_ph, corr=noise[0], sigma=noise[1])
       ysw = _ysw
       ysw_ph = _ysw_ph
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
       target2 = Targets.RayleighDispersionPhase(xsw_ph, ysw_ph, yerr=ysw_err_ph)
    #target2 = Targets.PReceiverFunction(xrf, yrf)
    #target2.moddata.plugin.set_modelparams(gauss=1., water=0.01, p=6.4)

    # Join the targets. targets must be a list instance with all targets
    # you want to use for MCMC Bayesian inversion.
       targets = Targets.JointTarget(targets=[target1, target2])


    #
    #  ---------------------------------------------------  Quick parameter update
    #
    # "priors" and "initparams" from config.ini are python dictionaries. You could
    # also simply define the dictionaries directly in the script, if you don't want
    # to use a config.ini file. Or update the dictionaries as follows, e.g. if you
    # have station specific values, etc.
    # See docs/bayhunter.pdf for explanation of parameters
       """
       priors.update({#'mohoest': (38, 4),  # optional, moho estimate (mean, std)
                   'rfnoise_corr': 0.98,
                   'swdnoise_corr': 0.
                   # 'rfnoise_sigma': np.std(yrf_err),  # fixed to true value
                   # 'swdnoise_sigma': np.std(ysw_err),  # fixed to true value
                   })

       initparams.update({'nchains': 5,
                       'iter_burnin': (2048 * 32),
                       'iter_main': (2048 * 16),
                       'propdist': (0.025, 0.025, 0.015, 0.005, 0.005),
                       })
        """

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
       optimizer.mp_inversion(nthreads=nthreads, baywatch=False, dtsend=1)


    #
    # #  ---------------------------------------------- Model resaving and plotting
       path = initparams['savepath']
       cfile = '%s_config.pkl' % initparams['station']
       configfile = op.join(path, 'data', cfile)
       obj = PlotFromStorage(configfile)
    # The final distributions will be saved with save_final_distribution.
    # Beforehand, outlier chains will be detected and excluded.
    # Outlier chains are defined as chains with a likelihood deviation
    # of dev * 100 % from the median posterior likelihood of the best chain.
       obj.save_final_distribution(maxmodels=10000000, dev=5)
    # Save a selection of important plots
       obj.save_plots(depint=0.01, nchains=50)
       obj.merge_pdfs()

    #
    # If you are only interested on the mean posterior velocity model, type:
       file = op.join(initparams['savepath'], 'data/c_models.npy')
       models = np.load(file)
       singlemodels = ModelMatrix.get_singlemodels(models, dep_int=np.arange(0, 100, 0.01))
       vs, dep = singlemodels['median']
       vs = gaussian_filter(vs, sigma=1)
       argx = np.argmin(np.abs(dep - max_depth))
       vs = vs[:argx+1]
       dep = dep[:argx+1]
       
    #
    # #  ---------------------------------------------- WATCH YOUR INVERSION
    # if you want to use BayWatch, simply type "baywatch ." in the terminal in the
    # folder you saved your baywatch configfile or type the full path instead
    # of ".". Type "baywatch --help" for further options.

    # if you give your public address as option (default is local address of PC),
    # you can also use BayWatch via VPN from 'outside'.
    # address = '139.?.?.?'  # here your complete address !!!
       subprocess.run(f"mv results_{jobID}-{arrayID} " "results_{:s}".format(os.path.basename(disp_file)), shell=True, check=True)
       subprocess.run("rm results_{:s}/data/c0*_*.npy".format(os.path.basename(disp_file)), shell=True, check=True)
       np.savetxt("./results_{:s}/model.txt".format(os.path.basename(disp_file)), np.column_stack((dep, vs)), delimiter=',', header='Depth(km),Vs(km/s)', comments='')

except Exception as e:
    print("Failed to compute.", e)
    sys.exit()
    #warnings.simplefilter("always", e)
    
