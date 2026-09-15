
from Censor.covnet201902 import covariance as cncovariance
from Censor.covnet201902 import data as cndata
import matplotlib.pyplot as plt

from Censor.display import display_coherencestack
from datetime import datetime, timedelta


def check_window_coherence(starttime=datetime(2018, 1, 1, 12, 00), timedelta_hours=1,
                           ffilter_min=0.01, ffilter_max=10,
                           fstack_min=0.01, fstack_max=0.1,
                           subwindow_seconds=30, subwindow_number=10,
                           work_path='', stream=None, coherence_threshold=0.9,
                           fig_coherencestack=False):

    ############################################################################
    ## PARAMETERS

    # Time
    endtime = starttime + timedelta(hours=timedelta_hours)
    # Frequency
    fdatafilter = [ffilter_min, ffilter_max]
    fdisplay = [ffilter_min, ffilter_max]
    # Window / Subwindow
    subwindow_overlap = 0.5
    window_overlap = 0.5
    # Preprocessing
    preproc_spectral = True
    preproc_spectral_method = 'whiten'
    preproc_spectral_submethod = 'onebit'
    preproc_temporal = False
    preproc_temporal_method = 'demad'
    preproc_spectral_overlap = 0.5
    preproc_spectral_secs = subwindow_seconds * subwindow_number * preproc_spectral_overlap
    # Paths
    fig_path = work_path + 'fig/'
    # Figures
    fig_format = '.png'
    fig_coherencestack_name = fig_path + 'coherence_' +\
                '{:04d}-{:03d}-{:02d}:{:02d}_'.format(starttime.year,
                starttime.timetuple().tm_yday, starttime.hour, starttime.minute) +\
                '{:04d}-{:03d}-{:02d}:{:02d}'.format(endtime.year,
                endtime.timetuple().tm_yday, endtime.hour, endtime.minute) +\
                '_f{:.2f}-{:.2f}Hz'.format(fdatafilter[0], fdatafilter[1]) +\
                '_winov{:.2f}'.format(window_overlap) +\
                '_swin{:d}s'.format(subwindow_seconds) +\
                '-n{:d}'.format(subwindow_number) +\
                '-ov{:.2f}'.format(subwindow_overlap) +\
                '_pp{:.1f}'.format(preproc_spectral_overlap) +\
                '_stack' +\
                '{:s}'.format(fig_format)

    ############################################################################
    ## CALCULATION

    # Data
    stream = cndata.read(stream)
    # if traces_merge is True:
    #     stream.merge(method=1, interpolation_samples=0, fill_value=0)
    print(stream)
    stream.cut(starttime=starttime, endtime=endtime, pad=True, fill_value=0)
    # stream.decimate(fdecimate)
    stream.detrend(type='demean')
    stream.filter(type='bandpass', freqmin=fdatafilter[0], freqmax=fdatafilter[1],
        zerophase=True)
    # Pre-processing
    if preproc_spectral is True:
        if preproc_spectral_method == 'whiten':
            stream.whiten(preproc_spectral_secs, method=preproc_spectral_submethod)
    if preproc_temporal is True:
        if preproc_temporal_method == 'binarize':
            stream.binarize(epsilon=1e-10)
        elif preproc_temporal_method == 'stationarize':
            stream.stationarize(length=11, order=1, epsilon=1e-10)
        elif preproc_temporal_method == 'demad':
            stream.demad()
    # Fourier Transform
    times, freqs, spectra = stream.fft(subwindow_seconds, bandwidth=None,
        step=subwindow_overlap)
    del stream
    # Covariance
    times_cov, covariance = cncovariance.calculate(times, spectra,
        average=subwindow_number, overlap=window_overlap)
    del times, spectra
    # Coherence
    coherence = cncovariance.spectral_width()
    del covariance
    # Coherence stack
    gs = dict(width_ratios=[50, 1], height_ratios=[50, 10], wspace=0.1, hspace=0.2)
    fig_cohstack, ax_cohstack = plt.subplots(2, 2, figsize=(11, 9), gridspec_kw=gs)
    ax_cohstack = ax_cohstack.ravel()
    covariance.show_coherence(times_cov, freqs, coherence, cmap='RdYlBu',
        ax=ax_cohstack[0], cax=ax_cohstack[1], rasterized=True)
    fig_cohstack_argslist = list([fig_cohstack, ax_cohstack, starttime, endtime,
        fdisplay, times_cov, freqs, coherence, subwindow_number, fstack_min, fstack_max])
    if fig_coherencestack is True:
        fig_cohstack_argsdict = dict(coherence_threshold=coherence_threshold,
            figure_file_name=fig_coherencestack_name)
    else:
        fig_cohstack_argsdict = dict(coherence_threshold=coherence_threshold)
    test = display_coherencestack.show_coherencestack(*fig_cohstack_argslist,
        **fig_cohstack_argsdict)

    return test
