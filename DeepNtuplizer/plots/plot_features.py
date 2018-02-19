#######
# plot distributions  of different featues of deepntuples root tree
# add root files with ntuples to files_mc or files_data
# itinitialize feature to plot and bin range
#   for twodimensional features you have to choose 'all' if you want to plot all flattened or 'max' to plot the maximum value of each tuple
#
#######

import pdb
import ROOT
import os
from array import array

import matplotlib.pyplot as plt
import root_numpy as rn
import numpy as np

ROOT.gROOT.SetBatch()  # don't pop up canvases
ROOT.gROOT.SetStyle('Plain')  # white background
ROOT.gStyle.SetFillStyle(0)

def hist_bin_uncertainty(data, weights, bin_edges):
    """
    The statistical uncertainity per bin of the binned data.
    If there are weights then the uncertainity will be the root of the
    sum of the weights squared.
    If there are no weights (weights = 1) this reduces to the root of
    the number of events.
    Args:
        data: `array`, the data being histogrammed.
        weights: `array`, the associated weights of the `data`.
        bin_edges: `array`, the edges of the bins of the histogram.
    Returns:
        bin_uncertainties: `array`, the statistical uncertainity on the bins.
    """
    # Bound the data and weights to be within the bin edges
    in_range_index = [idx for idx in range(len(data))
                      if data[idx] > min(bin_edges) and data[idx] < max(bin_edges)]
    in_range_data = np.asarray([data[idx] for idx in in_range_index])

    if weights is None or np.array_equal(weights, np.ones(len(weights))):
        # Default to weights of 1 and thus uncertainty = sqrt(N)
        in_range_weights = np.ones(len(in_range_data))
    else:
        in_range_weights = np.asarray([weights[idx] for idx in in_range_index])

    # Bin the weights with the same binning as the data
    bin_index = np.digitize(in_range_data, bin_edges)
    # N.B.: range(1, bin_edges.size) is used instead of set(bin_index) as if
    # there is a gap in the data such that a bin is skipped no index would appear
    # for it in the set
    binned_weights = np.asarray(
        [in_range_weights[np.where(bin_index == idx)[0]] for idx in range(1, len(bin_edges))])
    bin_uncertainties = np.asarray(
        [np.sqrt(np.sum(np.square(w))) for w in binned_weights])
    return bin_uncertainties



class datacollection:
    """ DOC """


    def __init__(self, filenames_mc, filenames_data = None):

        self.features = []
        self.bins_list = []
        self.fillopt_list = []
        self.freqs_list = []
        self.scalefactor_data = 1.
        self.filenames_mc = filenames_mc
        self.filenames_data = filenames_data


        #pdb.set_trace()


    def fill(self):

        print 'load samples'
        self.collection_mc = rn.root2array(self.filenames_mc, "deepntuplizer/tree", list(set(self.features)))
        self.weight_mc = rn.root2array(self.filenames_mc, "deepntuplizer/tree", "event_weight")
        if self.filenames_data != None:
            self.collection_data = rn.root2array(self.filenames_data, "deepntuplizer/tree", list(set(self.features)))
            self.weight_data = rn.root2array(self.filenames_data, "deepntuplizer/tree", "event_weight")

        #pdb.set_trace()


    def add_feature(self,feature, bins, fillopt=''):
        self.features.append(feature)
        self.bins_list.append(bins)
        self.fillopt_list.append(fillopt)



    def scale_data(self,sf):
        self.scalefactor_data *= sf


    def plot_features(self):
        #pdb.set_trace()
        for i, var in enumerate(self.features):
            print 'make plot for feature '+var

            if 'all' in self.fillopt_list[i]:
                flat_col_mc = np.concatenate(self.collection_mc[var], axis=-1)
                freqlist = [len(j) for j in self.collection_mc[var]]
                flat_weight_mc = np.repeat(self.weight_mc, freqlist)

                flat_col_data = np.concatenate(self.collection_data[var], axis=-1)
                freqlist = [len(j) for j in self.collection_data[var]]
                flat_weight_data = np.repeat(self.weight_data, freqlist)

                entries_mc, binEdges = np.histogram(flat_col_mc, self.bins_list[i], weights=flat_weight_mc)
                entries_data, _ = np.histogram(flat_col_data, self.bins_list[i], weights=flat_weight_data)
                err_mc = hist_bin_uncertainty(flat_col_mc, flat_weight_mc, binEdges)
                err_data = hist_bin_uncertainty(flat_col_data, flat_weight_data, binEdges)

            if 'max' in self.fillopt_list[i]:
                max_col_mc = np.empty((self.collection_mc[var].shape))
                max_col_data = np.empty((self.collection_data[var].shape))
                for j,arr in enumerate(self.collection_mc[var]):
                    if len(arr) > 0:
                        max_col_mc[j] = np.amax(arr)

                for j, arr in enumerate(self.collection_data[var]):
                    if len(arr) > 0:
                        max_col_data[j] = np.amax(arr)
                entries_mc, binEdges = np.histogram(max_col_mc, self.bins_list[i], weights=self.weight_mc)
                entries_data, _ = np.histogram(max_col_data, self.bins_list[i], weights=self.weight_data)
                err_mc = hist_bin_uncertainty(max_col_mc, self.weight_mc, binEdges)
                err_data = hist_bin_uncertainty(max_col_data, self.weight_data, binEdges)

            if self.fillopt_list[i] == '':
                entries_mc, binEdges = np.histogram(self.collection_mc[var], self.bins_list[i], weights=self.weight_mc)
                entries_data, _ = np.histogram(self.collection_data[var], self.bins_list[i], weights=self.weight_data)
                err_mc = hist_bin_uncertainty(self.collection_mc[var],self.weight_mc,binEdges)
                err_data = hist_bin_uncertainty(self.collection_data[var], self.weight_data, binEdges)

            plt.cla()
            fig, ax = plt.subplots()

            bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
            binWidth = binEdges[1:] - binEdges[:-1]
            plt.bar(bincenters, entries_mc, yerr = err_mc, width = binWidth, align='center', color='w')
            plt.errorbar(bincenters,entries_data*self.scalefactor_data, yerr=err_data*self.scalefactor_data, fmt='k.',ecolor='k')
            ax.set_xlim(binEdges[0], binEdges[-1])
            ax.set_ylabel('weighted frequency')
            ax.set_xlabel(var)
            ax.set_title(self.fillopt_list[i]+' entries of '+var+' in GoodJets')
            ax.set_ylim(bottom=0)

            plt.savefig(var+"_"+self.fillopt_list[i]+".png")






files_mc = ["tt_1_0.root","tt_2_0.root","tt_3_0.root","tt_4_0.root","tt_5_0.root",
            "dy50_1_0.root","dy10to50_1_0.root","wantit_1_0.root","wt_1_0.root",
            "ww_1_0.root","wz_1_0.root","zz_1_0.root","wjets_1_0.root"]
files_data = "muonEG_H_0.root"

dc = datacollection(files_mc,files_data)
dc.add_feature("nCpfcand", np.arange(0,30,1))
dc.add_feature("nNpfcand", np.arange(0,30,1))
dc.add_feature("jet_pt", np.arange(0,320,20))
dc.add_feature("jet_eta", np.arange(-2.4,2.8,0.4))
dc.add_feature("npv", np.arange(0,50,1))
dc.add_feature("nsv", np.arange(0,10,1))

dc.add_feature("Cpfcan_pt", np.arange(0,105,5),fillopt='max')
dc.add_feature("Npfcan_pt", np.arange(0,105,5),fillopt='max')
dc.add_feature("sv_pt", np.arange(0,110,10),fillopt='max')

dc.add_feature("Cpfcan_pt", np.arange(0,55,2.5),fillopt='all')
dc.add_feature("Npfcan_pt", np.arange(0,55,2.5),fillopt='all')
dc.add_feature("sv_pt", np.arange(0,110,10),fillopt='all')
dc.add_feature("sv_deltaR", np.arange(0,110,10),fillopt='all')





dc.fill()


dc.scale_data(1./0.8701594758966494)#*35.9/8.651)


title='features'

directory = os.path.dirname('./plots_'+title+'/')
# make a canvas, draw, and save it
if not os.path.exists(directory):
    os.makedirs(directory)
os.chdir(directory)

dc.plot_features()