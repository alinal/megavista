
import vista_utils as tsv # Get it at: https://github.com/arokem/vista_utils
from nitime.fmri.io import time_series_from_file as load_nii
import nitime.timeseries as ts
import nitime.viz as viz

from nitime.analysis import CorrelationAnalyzer, CoherenceAnalyzer, SpectralAnalyzer, FilterAnalyzer
#Import utility functions:
from nitime.utils import percent_change
from nitime.viz import drawmatrix_channels, drawgraph_channels
from nitime.viz import drawmatrix_channels, drawgraph_channels, plot_xcorr



def bp_data(array_data, TR, f_ub, f_lb):
    # Pass in time series data (voxels x TR)
    ts_data=ts.TimeSeries(array_data, sampling_interval=TR)
    ts_data_anal=FilterAnalyzer(ts_data, lb=f_lb, ub=f_ub)
    ts_box=ts_data_anal.filtered_boxcar

    return ts_box




