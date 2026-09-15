# coding: utf-8

import matplotlib.dates as dates
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def show_coherencestack(fig_coherencestack, ax_coherencestack, starttime,
                        endtime, fdisplay, times, frequencies,
                        coherence, Nsta, fstack_min, fstack_max,
                        coherence_threshold=0.9, figure_file_name=None,
                        **kwargs):
    # Parameters
    dpi = 300
    labelsize = 16
    labelsize_stack = 14
    tickslabelsize = 14
    detectionssize = 8
    yticks_stack = np.linspace(0.7, 1.0, num=3, endpoint=True)

    # Extract
    delta_coherencetime = coherence.shape[0] - len(times)
    if delta_coherencetime != 0:
        times = times[0:delta_coherencetime]
    fstack_min_id = np.argmin(np.abs(frequencies - fstack_min))
    fstack_max_id = np.argmin(np.abs(frequencies - fstack_max))
    coherence_stack = coherence[:, fstack_min_id : fstack_max_id + 1]
    coherence_stack_sum = np.sum(coherence_stack, axis=1)
    coherence_stack_norm = Nsta * coherence_stack.shape[1]
    coherence_stack = 1 - coherence_stack_sum / coherence_stack_norm

    ax_coherencestack[0].axhline(fstack_min, xmin=0, xmax=1, c='black', ls='--')
    ax_coherencestack[0].axhline(fstack_max, xmin=0, xmax=1, c='black', ls='--')
    ax_coherencestack[2].plot(times, coherence_stack, c='black', ls='-')
    ax_coherencestack[2].axhline(coherence_threshold, xmin=0, xmax=1, c='r')

    # Cosmetics
    ax_coherencestack[0].set_xlim([starttime, endtime])
    ax_coherencestack[2].set_xlim([starttime, endtime])
    ax_coherencestack[0].set_ylim(fdisplay)
    ax_coherencestack[2].set_ylim((yticks_stack[0], yticks_stack[-1]))
    ax_coherencestack[0].xaxis.set_major_locator(dates.AutoDateLocator())
    ax_coherencestack[2].xaxis.set_major_locator(dates.AutoDateLocator())
    ax_coherencestack[0].set_xticklabels([])
    ax_coherencestack[2].xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    ax_coherencestack[2].set_yticks(yticks_stack)
    ax_coherencestack[2].set_yticklabels(yticks_stack)
    ax_coherencestack[2].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
    ax_coherencestack[2].set_xlabel(starttime.strftime('%Y-%m-%d'), fontsize=labelsize)
    ax_coherencestack[0].set_ylabel('Frequency (Hz)', fontsize=labelsize)
    ax_coherencestack[2].set_ylabel('Stack \n {:.2f}-{:.2f}Hz'.format(fstack_min, fstack_max), fontsize=labelsize_stack)
    ax_coherencestack[1].set_ylabel('Spectral width', fontsize=labelsize)
    ax_coherencestack[0].tick_params(which='both', direction='out', labelsize=tickslabelsize)
    ax_coherencestack[1].tick_params(which='both', direction='out', labelsize=tickslabelsize)
    ax_coherencestack[2].tick_params(which='both', direction='out', labelsize=tickslabelsize)
    ax_coherencestack[3].axis('off')

    # Save
    if figure_file_name is not None:
        fig_stream_coh_stack.savefig(figure_file_name, dpi=dpi, bbox_inches='tight')
    #else:
    #    plt.show()
    plt.close('all')

    test = np.max(coherence_stack)
    return test
