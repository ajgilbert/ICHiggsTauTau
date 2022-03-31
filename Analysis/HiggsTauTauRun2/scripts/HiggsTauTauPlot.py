import ROOT
import os
import glob
import json
from UserCode.ICHiggsTauTau.analysis import *
from UserCode.ICHiggsTauTau.uncertainties import ufloat
from optparse import OptionParser
import argparse
import ConfigParser
import UserCode.ICHiggsTauTau.plotting as plotting
from collections import OrderedDict
import copy

CHANNELS= ['et', 'mt', 'em','tt','zmm','zee','mj']
ANALYSIS= ['sm','mssm','Hhh']
METHODS= [8 ,9, 10, 11, 12 , 13, 14, 15, 16, 17, 18, 19]

conf_parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False
    )
conf_parser.add_argument("--cfg",
                    help="Specify config file", metavar="FILE")
options, remaining_argv = conf_parser.parse_known_args()

defaults = { 
    "channel":"mt" , "outputfolder":"output",
    "folder":"/vols/cms/dw515/Offline/output/MSSM/Jan11/" , "signal_folder":"", "embed_folder":"",
    "paramfile":"scripts/Params_2016_spring16.json", "cat":"inclusive", "year":"2016",
    "era":"mssmsummer16", "sel":"(1)", "set_alias":[], "analysis":"mssm", "var":"m_vis(7,0,140)",
    "method":8 , "do_ss":False, "sm_masses":"125", "ggh_masses":"", "bbh_masses":"",
    "bbh_nlo_masses":"", "nlo_qsh":False, "qcd_os_ss_ratio":-1, "add_sm_background":"",
    "syst_e_scale":"", "syst_mu_scale":"", "syst_tau_scale":"", "syst_tau_scale_0pi":"",
    "syst_e_res":"", "syst_mu_res":"", "syst_tau_res":"", 
    "syst_tau_scale_1pi":"", "syst_tau_scale_3prong":"", "syst_tau_scale_3prong1pi0":"", "syst_eff_t":"", "syst_tquark":"",
    "syst_zwt":"", "syst_w_fake_rate":"", "syst_scale_j":"", "syst_res_j":"", "syst_scale_j_rbal":"",
    "syst_scale_j_rsamp":"", "syst_scale_j_full":"", "syst_scale_j_cent":"", "syst_scale_j_hf":"", 
    "syst_scale_j_full_corr":"", "syst_scale_j_cent_corr":"", "syst_scale_j_hf_corr":"",
    "syst_scale_j_full_uncorr":"", "syst_scale_j_cent_uncorr":"", "syst_scale_j_hf_uncorr":"",
    "syst_scale_j_by_source":"","jes_sources":"1:27", "syst_eff_b":"", "syst_fake_b":"", "syst_eff_b_weights":"",
    "norm_bins":False, "blind":False, "x_blind_min":100, "x_blind_max":4000, "ratio":False,
    "y_title":"", "x_title":"", "custom_y_range":False, "y_axis_min":0.001,
    "y_axis_max":100,"custom_x_range":False, "x_axis_min":0.001, "x_axis_max":100, "log_x":False,
    "log_y":False, "extra_pad":0.0, "signal_scale":1, "draw_signal_mass":"", "draw_signal_tanb":10,
    "signal_scheme":"run2_mssm", "lumi":"12.9 fb^{-1} (13 TeV)", "no_plot":False,
    "ratio_range":"0.7,1.3", "datacard":"", "do_custom_uncerts":False, "uncert_title":"Background uncertainty", 
    "custom_uncerts_wt_up":"","custom_uncerts_wt_down":"", "add_flat_uncert":0,
    "add_stat_to_syst":False, "add_wt":"", "custom_uncerts_up_name":"", "custom_uncerts_down_name":"",
    "do_ff_systs":False, "syst_efake_0pi_scale":"", "syst_efake_1pi_scale":"",
    "syst_mufake_0pi_scale":"", "syst_mufake_1pi_scale":"", "scheme":"","scheme":"", "syst_zpt_es":"",
    "syst_zpt_tt":"", "syst_zpt_statpt0":"", "syst_zpt_statpt40":"", "syst_zpt_statpt80":"",
    "syst_jfake_m":"", "syst_jfake_e":"", "syst_z_mjj":"", "syst_qcd_scale":"", "syst_quarkmass":"", "syst_mssm_ggh":False,
    "syst_ps":"", "syst_ue":"", "doNLOScales":False, "gen_signal":False, "doPDF":False,
    "doMSSMReWeighting":False, "do_unrolling":1, "syst_tau_id_dm0":"", "syst_tau_id_dm1":"",
    "syst_tau_id_dm10":"", "syst_lfake_dm0":"","syst_lfake_dm1":"","syst_qcd_shape_wsf":"",
    "syst_scale_met_unclustered":"","syst_scale_met_clustered":"",
    "extra_name":"", "no_default":False, "no_systs":False, "embedding":False,"syst_embedding_tt":"", "syst_embed_pt":"",
    "vbf_background":False, "syst_em_qcd":"", "syst_prefire":"",
    "syst_scale_met":"", "syst_res_met":"", "split_sm_scheme": False,
    "ggh_scheme": "powheg", "symmetrise":False, "mergeXbins":False, 'em_qcd_weight':"",
    "syst_scale_j_corr":"","syst_scale_j_uncorr":"", "syst_qcd_bkg":"",
    "ff_ss_closure":False, "threePads":False,"auto_blind":False,
    "syst_tau_id_diff":"", "syst_tau_trg_diff":"","syst_lep_trg_diff":"",
    "syst_scale_j_regrouped":"", "syst_tau_scale_grouped":"","wp":"medium","singletau":False,"qcd_ff_closure":False,
    "w_ff_closure":False,"ggh_masses_powheg":"", "bbh_masses_powheg":"", "vlq_sig":"","ratio_log_y":False,"plot_signals":""

}

if options.cfg:
    config = ConfigParser.SafeConfigParser()
    config.read([options.cfg])
    defaults.update(dict(config.items("Defaults")))

parser = argparse.ArgumentParser(
    parents=[conf_parser]
    )
parser.set_defaults(**defaults)
parser.add_argument("--channel", dest="channel", type=str,
    help="Tau decay channel to process.  Supported channels: %(CHANNELS)s" % vars())
parser.add_argument("--outputfolder", dest="outputfolder", type=str,
    help="Name of output folder")
parser.add_argument("--folder", dest="folder", type=str,
    help="Name of input folder")
parser.add_argument("--signal_folder", dest="signal_folder", type=str,
    help="If specified will use as input folder for signal samples, else will use same directroy specified by \"folder\" option.")
parser.add_argument("--embed_folder", dest="embed_folder", type=str,
    help="If specified will use as input folder for embed samples, else will use same directroy specified by \"folder\" option.")
parser.add_argument("--paramfile", dest="paramfile", type=str,
    help="Name of parameter file")
parser.add_argument("--cat", dest="cat", type=str,
    help="Category")
parser.add_argument("--datacard", dest="datacard", type=str,
    help="Datacard name")
parser.add_argument("--year", dest="year", type=str,
    help="Year")
parser.add_argument("--era", dest="era", type=str,
    help="Era")
parser.add_argument("--sel", dest="sel", type=str,
    help="Selection")
parser.add_argument("--set_alias", action="append", dest="set_alias", type=str,
    help="Overwrite alias selection using this options. Specify with the form --set_alias=nameofaliastoreset:newselection")
parser.add_argument("--analysis", dest="analysis", type=str, 
    help="Analysis.  Supported options: %(CHANNELS)s" % vars())
parser.add_argument("--var", dest="var", type=str,
    help="Variable to plot")
parser.add_argument("--method", dest="method", type=int,
    help="Method.  Supported options: %(METHODS)s" % vars())
parser.add_argument("--do_ss", dest="do_ss", action='store_true',
    help="Do same-sign.")
parser.add_argument("--sm_masses", dest="sm_masses", type=str,
    help="Comma seperated list of SM signal masses.")
parser.add_argument("--ggh_masses", dest="ggh_masses", type=str,
    help="Comma seperated list of SUSY ggH signal masses.")
parser.add_argument("--bbh_nlo_masses", dest="bbh_nlo_masses", type=str,
    help="Comma seperated list of SUSY NLO bbH signal masses.")
parser.add_argument("--nlo_qsh", dest="nlo_qsh", action='store_true',
    help="Do the Up/Down Qsh variations for NLO samples.")
parser.add_argument("--doNLOScales", dest="doNLOScales", action='store_true',
    help="Do the Up/Down QCD scale variations for NLO samples and compute uncertainties.")
parser.add_argument("--doPDF", dest="doPDF", action='store_true',
    help="Do PDF and alphaS variations for NLO samples and compute uncertainties.")
parser.add_argument("--doMSSMReWeighting", dest="doMSSMReWeighting", action='store_true',
    help="Do mA-tanb dependent reweighting of MSSM ggH signal.")
parser.add_argument("--bbh_masses", dest="bbh_masses", type=str,
    help="Comma seperated list of SUSY bbH signal masses.")
parser.add_argument("--qcd_os_ss_ratio", dest="qcd_os_ss_ratio", type=float,
    help="QCD OS/SS ratio")
parser.add_argument("--add_sm_background", dest="add_sm_background", type=str,
    help="Add SM Higgs background for MSSM")
parser.add_argument("--syst_tau_scale", dest="syst_tau_scale", type=str,
    help="If this string is set then the systematic shift due to tau energy scale is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_e_scale", dest="syst_e_scale", type=str,
    help="If this string is set then the systematic shift due to electron energy scale is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_mu_scale", dest="syst_mu_scale", type=str,
    help="If this string is set then the systematic shift due to muon energy scale is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tau_scale_0pi", dest="syst_tau_scale_0pi", type=str,
    help="If this string is set then the systematic shift due to the 1 prong 0 pi tau energy scale is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tau_scale_1pi", dest="syst_tau_scale_1pi", type=str,
    help="If this string is set then the systematic shift due to the 1 prong 1 pi tau energy scale is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tau_scale_3prong", dest="syst_tau_scale_3prong", type=str,
    help="If this string is set then the systematic shift due to 3 prong tau energy scale is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tau_scale_3prong1pi0", dest="syst_tau_scale_3prong1pi0", type=str,
    help="If this string is set then the systematic shift due to 3 prong + 1 pi0 tau energy scale is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tau_scale_grouped", type=str,
    help="If this string is set then the systematic shift due to tau energy scale is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_e_res", dest="syst_e_res", type=str,
    help="If this string is set then the systematic shift due to electron energy resolution is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_mu_res", dest="syst_mu_res", type=str,
    help="If this string is set then the systematic shift due to muon energy resolution is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tau_res", dest="syst_tau_res", type=str,
    help="If this string is set then the systematic shift due to tau energy resolution is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_eff_t", dest="syst_eff_t", type=str, default='',
    help="If this string is set then the systematic shift due to tau ID is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tquark", dest="syst_tquark", type=str,
    help="If this string is set then the top-quark weight systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_zwt", dest="syst_zwt", type=str,
    help="If this string is set then the z-reweighting systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_w_fake_rate", dest="syst_w_fake_rate", type=str, default='',
    help="If this string is set then the W+jets fake-rate systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_scale_j", dest="syst_scale_j", type=str,
    help="If this string is set then the jet scale systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_res_j", dest="syst_res_j", type=str,
    help="If this string is set then the jER systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_scale_j_corr", dest="syst_scale_j_corr", type=str,
    help="If this string is set then the jet scale systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_scale_j_uncorr", dest="syst_scale_j_uncorr", type=str,
    help="If this string is set then the jet scale systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_scale_j_rbal", dest="syst_scale_j_rbal", type=str,
    help="If this string is set then the RelativeBal jet scale systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_scale_j_rsamp", dest="syst_scale_j_rsamp", type=str,
    help="If this string is set then the RelativeSample jet scale systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_scale_j_full", dest="syst_scale_j_full", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = full region (eta<5)")
parser.add_argument("--syst_scale_j_cent", dest="syst_scale_j_cent", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = central region (eta<3)")
parser.add_argument("--syst_scale_j_hf", dest="syst_scale_j_hf", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = full region (eta>3)")
parser.add_argument("--syst_scale_j_full_corr", dest="syst_scale_j_full_corr", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = full region (eta<5). 2016/2017 correlated part.")
parser.add_argument("--syst_scale_j_cent_corr", dest="syst_scale_j_cent_corr", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = central region (eta<3). 2016/2017 correlated part.")
parser.add_argument("--syst_scale_j_hf_corr", dest="syst_scale_j_hf_corr", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = full region (eta>3). 2016/2017 correlated part.")
parser.add_argument("--syst_scale_j_full_uncorr", dest="syst_scale_j_full_uncorr", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = full region (eta<5). 2016/2017 correlated part.")
parser.add_argument("--syst_scale_j_cent_uncorr", dest="syst_scale_j_cent_uncorr", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = central region (eta<3). 2016/2017 correlated part.")
parser.add_argument("--syst_scale_j_hf_uncorr", dest="syst_scale_j_hf_uncorr", type=str,
    help="If this string is set then the regional jet scale systematic is performed with the set string appended to the resulting histogram name. Region = full region (eta>3). 2016/2017 correlated part.")
parser.add_argument("--syst_scale_j_by_source", dest="syst_scale_j_by_source", type=str,
    help="If this string is set then the jet scale systematic is performed split by source with the set string appended to the resulting histogram name. The string should contrain the substring  \'SOUCE\' which will be replaced by the JES source name")
parser.add_argument("--syst_scale_j_regrouped", dest="syst_scale_j_regrouped", type=str,
    help="If this string is set then the set of regrouped jet scale systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--jes_sources", dest="jes_sources", type=str,
    help="JES sources to process specified by integers seperated by commas. Values seperated by x\':\'y will process all integers from x to y. e.g using --jes_sources=1:3,10 will process sources: 1,2,3,10")
parser.add_argument("--syst_eff_b", dest="syst_eff_b", type=str,
    help="If this string is set then the b-tag efficiency systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_eff_b_weights", dest="syst_eff_b_weights", type=str,
    help="If this string is set then the b-tag efficiency systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_fake_b", dest="syst_fake_b", type=str,
    help="If this string is set then the b-tag fake-rate systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--norm_bins", dest="norm_bins", action='store_true',
    help="Normalize bins by bin width.")
parser.add_argument("--blind", dest="blind", action='store_true',
    help="Blind histogram.")
parser.add_argument("--x_blind_min", dest="x_blind_min", type=float,
    help="Minimum x for blinding.")
parser.add_argument("--x_blind_max", dest="x_blind_max", type=float,
    help="Maximum x for blinding.")
parser.add_argument("--ratio", dest="ratio", action='store_true',
    help="Draw ratio.")
parser.add_argument("--y_title", dest="y_title", type=str,
    help="Y-axis title.")
parser.add_argument("--x_title", dest="x_title", type=str,
    help="X-axis title.")
parser.add_argument("--custom_y_range", dest="custom_y_range", action='store_true',
    help="Use custom y-axis range")
parser.add_argument("--y_axis_min", dest="y_axis_min", type=float,
    help="Minimum y-axis value.")
parser.add_argument("--y_axis_max", dest="y_axis_max", type=float,
    help="Maximum y-axis value.")
parser.add_argument("--custom_x_range", dest="custom_x_range", action='store_true',
    help="Use custom x-axis range")
parser.add_argument("--x_axis_min", dest="x_axis_min", type=float,
    help="Minimum x-axis value.")
parser.add_argument("--x_axis_max", dest="x_axis_max", type=float,
    help="Maximum x-axis value.")
parser.add_argument("--log_x", dest="log_x", action='store_true',
    help="Set log scale on x-axis.")
parser.add_argument("--log_y", dest="log_y", action='store_true',
    help="Set log scale on y-axis.")
parser.add_argument("--ratio_log_y", dest="ratio_log_y", action='store_true',
    help="Set log scale on ratio y-axis.")
parser.add_argument("--extra_pad", dest="extra_pad", type=float,
    help="Fraction of extra whitespace at top of plot.")
parser.add_argument("--signal_scale", dest="signal_scale", type=float,
    help="Signal scale.")
parser.add_argument("--draw_signal_mass", dest="draw_signal_mass", type=str,
    help="Signal mass.")
parser.add_argument("--draw_signal_tanb", dest="draw_signal_tanb", type=float,
    help="Signal tanb.")
parser.add_argument("--signal_scheme", dest="signal_scheme", type=str,
    help="Signal scale.")
parser.add_argument("--lumi", dest="lumi", type=str,
    help="Lumi.")
parser.add_argument("--no_plot", dest="no_plot", action='store_true',
    help="If option is set then no pdf or png plots will be created only the output root file will be produced.")
parser.add_argument("--ratio_range", dest="ratio_range", type=str,
    help="y-axis range for ratio plot in format MIN,MAX")
parser.add_argument("--do_custom_uncerts", dest="do_custom_uncerts", action='store_true',
    help="Do custom uncertainty band. Up and down weights for this uncertainty band should be set using \"custom_uncerts_wt_up\" and \"custom_uncerts_wt_down\" options")
parser.add_argument("--custom_uncerts_down_name", dest="custom_uncerts_down_name", type=str,
    help="Name of histogram to use for uncertainty down band")
parser.add_argument("--custom_uncerts_up_name", dest="custom_uncerts_up_name", type=str,
    help="Name of histogram to use for uncertainty up band")
parser.add_argument("--custom_uncerts_wt_up", dest="custom_uncerts_wt_up", type=str,
    help="Up weight for custom uncertainty band")
parser.add_argument("--custom_uncerts_wt_down", dest="custom_uncerts_wt_down", type=str,
    help="Down weight for custom uncertainty band")
parser.add_argument("--uncert_title", dest="uncert_title", type=str,
    help="Custom uncertainty band legend label")
parser.add_argument("--add_stat_to_syst", dest="add_stat_to_syst", action='store_true',
    help="Add custom uncertainty band to statistical uncertainty.")
parser.add_argument("--add_flat_uncert", dest="add_flat_uncert", type=float,
    help="If set to non-zero will add a flat uncertainty band in quadrature to the uncertainty.")
parser.add_argument("--add_wt", dest="add_wt", type=str,
    help="Name of additional weight to be applied to all templates.")
parser.add_argument("--do_ff_systs", dest="do_ff_systs", action='store_true',
    help="Do fake-factor systamatic shifts.")
parser.add_argument("--syst_efake_0pi_scale", dest="syst_efake_0pi_scale", type=str,
    help="If this string is set then the e->tau dm=0 fake-rate systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_efake_1pi_scale", dest="syst_efake_1pi_scale", type=str,
    help="If this string is set then the e->tau dm=1 fake-rate systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_mufake_0pi_scale", dest="syst_mufake_0pi_scale", type=str,
    help="If this string is set then the mu->tau dm=0 fake-rate systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_mufake_1pi_scale", dest="syst_mufake_1pi_scale", type=str,
    help="If this string is set then the mu->tau dm=1 fake-rate systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--scheme", dest="scheme", type=str,
    help="Set plotting scheme")
parser.add_argument("--syst_zpt_es", dest="syst_zpt_es", type=str,
    help="If this string is set then the zpT muon ES systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_zpt_tt", dest="syst_zpt_tt", type=str,
    help="If this string is set then the zpT tt X-section systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_zpt_statpt0", dest="syst_zpt_statpt0", type=str,
    help="If this string is set then the zpT statistical pt0 systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_zpt_statpt40", dest="syst_zpt_statpt40", type=str,
    help="If this string is set then the zpT statistical pt40 systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_zpt_statpt80", dest="syst_zpt_statpt80", type=str,
    help="If this string is set then the zpT statistical pt80 systematic is performed with the set string appended to the resulting histogram name")
parser.add_argument("--syst_jfake_e", dest="syst_jfake_e", type=str,
    help="If set, adds the e->jet fake rate uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_jfake_m", dest="syst_jfake_m", type=str,
    help="If set, adds the e->jet fake rate uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_z_mjj", dest="syst_z_mjj", type=str,
    help="If set then add the uncertainty on the Z mjj corrections")
parser.add_argument("--syst_qcd_scale", dest="syst_qcd_scale", type=str,
    help="If set then add the qcd scale uncertainty for ggH")
parser.add_argument("--syst_quarkmass", dest="syst_quarkmass", type=str,
    help="If set then add the finite quark mass uncertainty for ggH")
parser.add_argument("--syst_ps", dest="syst_ps", type=str,
    help="If set then add the PS uncertainty for ggH")
parser.add_argument("--syst_ue", dest="syst_ue", type=str,
    help="If set then add the UE uncertainty for ggH")
parser.add_argument("--syst_mssm_ggh", dest="syst_mssm_ggh", action='store_true',
    help="If set then add the ggH uncertainties for the MSSM (scale and Hdamp uncerts)")
parser.add_argument("--syst_tau_id_dm0", dest="syst_tau_id_dm0", type=str,
    help="If set, adds the tau dm = 0 id uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tau_id_dm1", dest="syst_tau_id_dm1", type=str,
    help="If set, adds the tau dm = 1 id uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_tau_id_dm10", dest="syst_tau_id_dm10", type=str,
    help="If set, adds the tau dm = 10 id uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_lfake_dm0", dest="syst_lfake_dm0", type=str,
    help="If set, adds the e/mu->tau dm = 0 fake-rate uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_lfake_dm1", dest="syst_lfake_dm1", type=str,
    help="If set, adds the e/mu->tau dm = 1 fake-rate uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_qcd_shape_wsf", dest="syst_qcd_shape_wsf", type=str,
    help="If set, adds QCD shape uncertainty relating to W subtraction from SS data with the set string appended to the resulting histogram name")
parser.add_argument("--syst_scale_met_unclustered", dest="syst_scale_met_unclustered", type=str,
    help="If set, adds the unclustered energy MET uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_scale_met_clustered", dest="syst_scale_met_clustered", type=str,
    help="If set, adds the clustered energy MET uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--gen_signal", dest="gen_signal", action='store_true',
    help="If set then use generator-level tree for signal")
parser.add_argument("--do_unrolling", dest="do_unrolling", type=int,
    help="If argument is set to true will unroll 2D histograms into 1D histogram.")
parser.add_argument("--no_default", dest="no_default", action='store_true',
    help="If option is speficied then don't do nominal histograms.")
parser.add_argument("--no_systs", dest="no_systs", action='store_true',
    help="If option is speficied then don't do systematics histograms.")
parser.add_argument("--extra_name", dest="extra_name", type=str,
    help="If set, adds an additional string to the output datacard name")
parser.add_argument("--embedding", dest="embedding", action='store_true',
    help="If option is speficied then use embedded samples for ZTT templates.")
parser.add_argument("--syst_embedding_tt", dest="syst_embedding_tt", type=str,
    help="If set, adds systematic templates for embedding corresponding to TTbar shift of +/-10\% ")
parser.add_argument("--syst_embed_pt", dest="syst_embed_pt", type=str,
    help="If set, adds systematic templates for embedding corresponding to Z->mumu non-closures ")
parser.add_argument("--vbf_background", dest="vbf_background", action='store_true',
    help="Add VBF to total background template")
parser.add_argument("--syst_em_qcd", dest="syst_em_qcd", type=str,
    help="If set, adds the QCD shape uncertainties for the em channel. You string should contain *BIN")
parser.add_argument("--syst_scale_met", dest="syst_scale_met", type=str,
    help="If set, adds the recoil corrected MET response uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--syst_res_met", dest="syst_res_met", type=str,
    help="If set, adds the recoil corrected MET resolution uncertainty with the set string appended to the resulting histogram name")
parser.add_argument("--split_sm_scheme", dest="split_sm_scheme", action='store_true',
    help="If set, splits the SM signal scheme into ggH, qqH and VH")
parser.add_argument("--ggh_scheme", dest="ggh_scheme", type=str,
    help="Decide which ggH scheme to plot with in split SM scheme mode (powheg or JHU)")
parser.add_argument("--symmetrise", dest="symmetrise", action='store_true',
    help="Use this option to symmetrise dijet Delta_phi bins in 2D histogram around centre bin")
parser.add_argument("--mergeXbins", dest="mergeXbins", action='store_true',
    help="Use this option to merge x bins to improve stats - to be used is the backgrounds are flat as a function of the X variable")
parser.add_argument("--em_qcd_weight", dest="em_qcd_weight", type=str,
    help="Define custom em QCD OSSS weight/function")
parser.add_argument("--syst_qcd_bkg", dest="syst_qcd_bkg", type=str,
    help="If set, adds systematic templates corresponding to shifting background subtraction in QCD method up/down by +/-10\% ")
parser.add_argument("--syst_prefire", dest="syst_prefire", type=str,
    help="If set, adds systematic templates corresponding to uncertainty on pre-firing correction.")
parser.add_argument("--syst_tau_id_diff", dest="syst_tau_id_diff", type=str,
    help="Do shape uncertainty corresponding to shifting the tau id SFs. The string you pass as the argument should contain either *DM, in which case the DM binned version will be used, or *PT in which base the pT binned version will be used.")
parser.add_argument("--syst_tau_trg_diff", dest="syst_tau_trg_diff", type=str,
    help="Do shape uncertainty corresponding to shifting the tau trigger SFs. The string you pass as the argument should contain either *DM which will be replaced with DMi for i=0,1,10,11")
parser.add_argument("--syst_lep_trg_diff", dest="syst_lep_trg_diff", type=str,
    help="Do shape uncertainty corresponding to shifting the lepton trigger SFs.")
parser.add_argument("--ff_ss_closure", dest="ff_ss_closure", action='store_true',
    help="If set then applies a non-closure correction to fake factor yields based on differences in SS data.")
parser.add_argument("--threePads", dest="threePads", action='store_true',
    help="If set then draws three pads (one ratio + additional).")
parser.add_argument("--wp", dest="wp", type=str,
    help="Tau ID working point to use (only works for mssmrun2).")
parser.add_argument("--singletau", dest="singletau", action='store_true',
    help="If set then use singletau baseline (only works for mssmrun2).")
parser.add_argument("--qcd_ff_closure", dest="qcd_ff_closure", action='store_true',
    help="Will run closure plot for QCD DR (mssmrun2)")
parser.add_argument("--w_ff_closure", dest="w_ff_closure", action='store_true',
    help="Will run closure plot for W DR (mssmrun2)")
parser.add_argument("--bbh_masses_powheg", dest="bbh_masses_powheg", type=str,
    help="SUSY bbh masses to run more powheg samples")
parser.add_argument("--ggh_masses_powheg", dest="ggh_masses_powheg", type=str,
    help="SUSY ggh masses to run more powheg samples")
parser.add_argument("--bkg_comp", dest="bkg_comp", action='store_true',
    help="Will plot the background composition on the 2nd pad and ratio on 3rd. Needs to be run with threePads.")
parser.add_argument("--ml_ff", dest="ml_ff", action='store_true',
    help="Use machine learning fake factors for mssmrun2.")
parser.add_argument("--vlq_sig", dest="vlq_sig", type=str,
    help="Comma separated list of signal parameter names i.e. vlq_betaRd33_minus1_mU4_gU1,vlq_betaRd33_minus1_mU4_gU2,vlq_betaRd33_minus1_mU4_gU3")
parser.add_argument("--plot_signals", dest="plot_signals", type=str,
    help="Comma separated list of what signals to plot")



options = parser.parse_args(remaining_argv)   

print ''
print '################### Options ###################'
print 'channel           = ' + options.channel
print 'outputfolder      = ' + options.outputfolder
print 'folder            = ' + options.folder
print 'paramfile         = ' + options.paramfile
print 'cat               = ' + options.cat
print 'datacard          = ' + options.datacard
print 'year              = ' + options.year
print 'era               = ' + options.era
print 'sel               = ' + options.sel
print 'analysis          = ' + options.analysis
print 'var               = ' + options.var
print 'method            ='  ,  options.method
print 'do_ss             ='  ,  options.do_ss
print 'sm_masses         = ' +  options.sm_masses
print 'ggh_masses        = ' +  options.ggh_masses
print 'bbh_masses        = ' +  options.bbh_masses
print 'qcd_os_ss_ratio   ='  ,  options.qcd_os_ss_ratio
print 'add_sm_background ='  ,  options.add_sm_background
print 'syst_tau_scale    ='  ,  options.syst_tau_scale
print 'syst_eff_t        ='  ,  options.syst_eff_t
print 'syst_tquark       ='  ,  options.syst_tquark
print 'syst_zwt          ='  ,  options.syst_zwt
print 'syst_w_fake_rate  ='  ,  options.syst_w_fake_rate
print 'syst_scale_j      ='  ,  options.syst_scale_j
print 'syst_eff_b        ='  ,  options.syst_eff_b
print 'syst_fake_b       ='  ,  options.syst_fake_b
print 'do_ff_systs       ='  ,  options.do_ff_systs
print 'singletau         ='  ,  options.singletau
print '###############################################'
print ''

# discrete x labels
discrete_x_axis = False
plot_var = options.var
discrete_x_labels=None
do_eq = False
do_geq = False
if  '[' in plot_var and ']' in plot_var:
  discrete_x_labels = plot_var.split('[')[1].split(']')[0].split(',')
  for i in discrete_x_labels:
    if ">=" in i:
      discrete_x_axis = True
      do_geq = True
    elif "==" in i:
      discrete_x_axis = True
      do_eq = True

if discrete_x_axis:
  if do_geq:
    run_bins = []
    for i in discrete_x_labels:
      if ">=" in i:
        run_bins.append(float(i.replace(">=","")))
      elif "==" in i:
        run_bins.append(float(i.replace("==","")))
      else:
        run_bins.append(float(i))
  elif do_eq:
    run_bins = range(int(discrete_x_labels[0].replace("==","")),int(discrete_x_labels[-1].replace("==",""))+2)  
  options.var = plot_var.split('[')[0]+str(run_bins)
else:
  discrete_x_labels = None 


# vbf_background = False
vbf_background = options.vbf_background

compare_w_shapes = False
compare_qcd_shapes = False
if options.scheme == "qcd_shape": compare_qcd_shapes = True
if options.scheme == "w_shape": compare_w_shapes = True
w_abs_shift=None # if not None then the QCD shape will be adjusted by shifting the W yield up and down by +/- w_abs_shift
if options.era in ["mssmsummer16","smsummer16","cpsummer16","cpdecay16","legacy16",'UL_16_preVFP','UL_16_postVFP',"tauid2016","mvadm2016"]: options.lumi = "35.9 fb^{-1} (13 TeV)"

# option to split the real tau events into rho, pi, a1, other
split_taus=False
if options.era == 'mvadm2016' and options.channel=='mt': split_taus=True

cats = {}
if options.analysis in ['sm','cpprod','cpdecay']:
    if options.channel == 'mt':
        cats['baseline'] = '(iso_1<0.15 && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto)'
        if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']: 
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_medium_2>0.5 && deepTauVsEle_vvloose_2>0.5 && deepTauVsMu_tight_2>0.5 && !leptonveto && ((trg_mutaucross&&pt_2>25&&pt_1<23)||(trg_singlemuon&&pt_1>23)))'
        if options.era in ['tauid2017']:
          cats['baseline'] = '(iso_1<0.15 && antiele_2 && antimu_2 && !leptonveto && pt_1>25 && trg_singlemuon &&pt_2>20)'
        if options.era in ['cpsummer17','UL_17','cp18','UL_18']:
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_medium_2>0.5 && deepTauVsEle_vvloose_2>0.5 && deepTauVsMu_tight_2>0.5 && !leptonveto && ((trg_mutaucross&&pt_2>32&&pt_1<25)||(trg_singlemuon&&pt_1>25)) && wt<2)'
        if options.era in ['tauid2016']: 
          cats['baseline'] = '(iso_1<0.15 && antiele_2 && antimu_2 && !leptonveto && trg_singlemuon && pt_1>23)'
          cats['baseline_loosemu'] = '(iso_1<0.15 && antiele_2 && antimu_loose_2 && !leptonveto && trg_singlemuon && pt_1>23)'
          cats['pass'] = 'mva_olddm_tight_2>0.5 && pzeta>-25'
          cats['fail'] = 'mva_olddm_tight_2<0.5 && pzeta>-25'
    elif options.channel == 'et': 
        cats['baseline'] = '(iso_1<0.1  && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto)'
        if options.era in ['smsummer16']: 
          cats['baseline'] = '(iso_1<0.1  && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto && trg_singleelectron)'
          cats['baseline_aisotau'] = '(iso_1<0.1 && mva_olddm_vloose_2>0.5 && mva_olddm_tight_2<0.5 && antiele_2 && antimu_2 && leptonveto==0 && trg_singleelectron && pt_2>30)'
        if options.era in ['tauid2016']: 
          cats['baseline'] = '(iso_1<0.1 && antiele_2 && antimu_2 && !leptonveto && trg_singleelectron)'
          cats['baseline_loosemu'] = '(iso_1<0.1 && antiele_2 && antimu_loose_2 && !leptonveto && trg_singleelectron)'
          cats['pass'] = 'mva_olddm_tight_2>0.5 && pzeta>-25'
          cats['fail'] = 'mva_olddm_tight_2<0.5 && pzeta>-25'
        if options.era in ['cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP']:
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_medium_2>0.5 && deepTauVsEle_tight_2>0.5 && deepTauVsMu_vloose_2>0.5 && leptonveto==0 && trg_singleelectron && fabs(wt<2))'
        if options.era in ['cpsummer17','UL_17']:
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_medium_2>0.5 && deepTauVsEle_tight_2>0.5 && deepTauVsMu_vloose_2>0.5 && leptonveto==0 && ((trg_etaucross&&pt_2>35&&pt_1<28&&fabs(eta_2)<2.1)||(trg_singleelectron&&pt_1>28)))'
        if options.era in ['cp18', 'UL_18']:
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_medium_2>0.5 && deepTauVsEle_tight_2>0.5 && deepTauVsMu_vloose_2>0.5 && !leptonveto && ((trg_etaucross&&pt_2>35&&pt_1<33&&fabs(eta_2)<2.1)||(trg_singleelectron&&pt_1>33)))'
          #cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_medium_2>0.5 && deepTauVsEle_tight_2>0.5 && deepTauVsMu_vloose_2>0.5 && ((trg_etaucross&&pt_2>35&&pt_1<33&&fabs(eta_2)<2.1 && fabs(eta_1-eta_2)>0.2)||(trg_singleelectron&&pt_1>33)))&&wt<2'
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_medium_2>0.5 && deepTauVsEle_tight_2>0.5 && deepTauVsMu_vloose_2>0.5 && ((trg_etaucross&&pt_2>35&&pt_1<33&&fabs(eta_2)<2.1)||(trg_singleelectron&&pt_1>33)))&&wt<2'

elif options.analysis in ['mssmrun2','vlq']:
    wp = options.wp

    if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']:
      m_lowpt=23
      e_lowpt=26
      t_highpt=120
      t_lowpt_mt=25
      t_lowpt_et=25
    if options.era in ['cpsummer17','UL_17']:
      m_lowpt=25
      e_lowpt=28
      t_highpt=180
      t_lowpt_mt=32
      t_lowpt_et=35
    if options.era in ['cp18','UL_18']:
      m_lowpt=25
      e_lowpt=33
      t_highpt=180
      t_lowpt_mt=32
      t_lowpt_et=35

    if options.channel == 'mt':
        if options.singletau:
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_%(wp)s_2>0.5 && deepTauVsEle_vvloose_2>0.5 && deepTauVsMu_tight_2>0.5 && !leptonveto && pt_2>30 && ((trg_mutaucross&&pt_2>%(t_lowpt_mt)s&&pt_2<%(t_highpt)s&&fabs(eta_2)<2.1&&pt_1<%(m_lowpt)s)||(trg_singlemuon&&pt_1>=%(m_lowpt)s)||(trg_singletau_2&&pt_2>=%(t_highpt)s&&fabs(eta_2)<2.1)))' % vars()
        else:
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_%(wp)s_2>0.5 && deepTauVsEle_vvloose_2>0.5 && deepTauVsMu_tight_2>0.5 && !leptonveto && pt_2>30 && ((trg_mutaucross&&pt_2>%(t_lowpt_mt)s&&fabs(eta_2)<2.1&&pt_1<%(m_lowpt)s)||(trg_singlemuon&&pt_1>=%(m_lowpt)s)))' % vars()

    if options.channel == 'et':
        if options.singletau:
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_%(wp)s_2>0.5 && deepTauVsEle_tight_2>0.5 && deepTauVsMu_vloose_2>0.5 && pt_2>30 && ((trg_etaucross&&pt_2>%(t_lowpt_et)s&&pt_2<%(t_highpt)s&&fabs(eta_2)<2.1&&pt_1<%(e_lowpt)s)||(trg_singleelectron&&pt_1>=%(e_lowpt)s)||(trg_singletau_2&&pt_2>=%(t_highpt)s&&fabs(eta_2)<2.1)))' % vars()
        else: 
          cats['baseline'] = '(iso_1<0.15 && deepTauVsJets_%(wp)s_2>0.5 && deepTauVsEle_tight_2>0.5 && deepTauVsMu_vloose_2>0.5 && pt_2>30 && ((trg_etaucross&&pt_2>%(t_lowpt_et)s&&fabs(eta_2)<2.1&&pt_1<%(e_lowpt)s)||(trg_singleelectron&&pt_1>=%(e_lowpt)s)))' % vars()
        
elif options.analysis == 'mssm':
    if options.channel == 'mt':        
        cats['baseline'] = '(iso_1<0.15 && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto)'
    elif options.channel == 'et':
        cats['baseline'] = '(iso_1<0.1  && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto)'
    if options.era == 'mssmsummer16':
        if options.channel == 'mt':        
            cats['baseline'] = '(iso_1<0.15 && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto)'
            cats['baseline_antiisotau'] = '(iso_1<0.15 && 1 && mva_olddm_tight_2<0.5 && antiele_2 && antimu_2 && !leptonveto && trg_singlemuon)'
            cats['ichep_baseline'] = '(iso_1<0.15 && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto && trg_singlemuon)'
        elif options.channel == 'et':
            cats['baseline'] = '(iso_1<0.1  && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto)'
            cats['baseline_antiisotau'] = '(iso_1<0.1 && mva_olddm_tight_2<0.5 && antiele_2 && antimu_2 && !leptonveto && trg_singleelectron)'
            cats['ichep_baseline'] = '(iso_1<0.1 && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto && trg_singleelectron)'
        elif options.channel == 'mj':        
            cats['baseline'] = '(iso_1<0.15 && !leptonveto)'
if options.channel == 'tt':
    cats['baseline'] = '(mva_olddm_tight_1>0.5 && mva_olddm_tight_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto)'
    if options.era == 'mssmsummer16': cats['baseline'] = '(mva_olddm_medium_1>0.5 && mva_olddm_medium_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto)'
    if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','tauid2016','mvadm2016']: 
        cats['baseline'] = '(mva_olddm_tight_1>0.5 && mva_olddm_tight_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto && trg_doubletau && pt_1>40)'
        cats['baseline_aisotau1'] = '(pt_1>50 && mva_olddm_vloose_1>0.5 && mva_olddm_tight_1<0.5 && mva_olddm_tight_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto && trg_doubletau)'
        cats['baseline_aisotau2'] = '(pt_1>50 && mva_olddm_vloose_2>0.5 && mva_olddm_tight_2<0.5 && mva_olddm_tight_1>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto && trg_doubletau)'
        cats['baseline_aisotau2_sb'] = '(mva_olddm_vloose_1>0.5 && mva_olddm_tight_1<0.5 && mva_olddm_tight_2<0.5 && mva_olddm_medium_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && leptonveto==0 && trg_doubletau)'
        cats['baseline_aisotau2_sb'] = '(mva_olddm_vloose_2>0.5 && mva_olddm_tight_2<0.5 && mva_olddm_tight_1<0.5 && mva_olddm_medium_1>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && leptonveto==0 && trg_doubletau)'
    if options.era in ['cpdecay16','legacy16','UL_16_preVFP','UL_16_postVFP']:
      cats['baseline'] = '(deepTauVsJets_medium_1>0.5 && deepTauVsJets_medium_2>0.5 && leptonveto==0 && trg_doubletau && deepTauVsEle_vvloose_1 && deepTauVsEle_vvloose_2 && deepTauVsMu_vloose_1 && deepTauVsMu_vloose_2)'
    if options.era in ['cpsummer17','UL_17','cp18','UL_18']:
      cats['baseline'] = '(deepTauVsJets_medium_1>0.5 && deepTauVsJets_medium_2>0.5 && leptonveto==0 && (trg_doubletau && pt_2>40) && deepTauVsEle_vvloose_1 && deepTauVsEle_vvloose_2 && deepTauVsMu_vloose_1 && deepTauVsMu_vloose_2)'

    if options.analysis in ['mssmrun2','vlq']:
      wp = options.wp
      if options.singletau:
        if options.era in['legacy16','UL_16_preVFP','UL_16_postVFP']:
          cats['baseline'] = '(deepTauVsJets_%(wp)s_1>0.5 && deepTauVsJets_%(wp)s_2>0.5 && leptonveto==0 && (trg_doubletau || (pt_1>120 && trg_singletau_1) || (pt_2>120 && trg_singletau_2)) && deepTauVsEle_vvloose_1 && deepTauVsEle_vvloose_2 && deepTauVsMu_vloose_1 && deepTauVsMu_vloose_2)' % vars()
        else:
          cats['baseline'] = '(deepTauVsJets_%(wp)s_1>0.5 && deepTauVsJets_%(wp)s_2>0.5 && leptonveto==0 && (trg_doubletau || (pt_1>180 && trg_singletau_1) || (pt_2>180 && trg_singletau_2)) && deepTauVsEle_vvloose_1 && deepTauVsEle_vvloose_2 && deepTauVsMu_vloose_1 && deepTauVsMu_vloose_2)' % vars()
      else:
        cats['baseline'] = '(deepTauVsJets_%(wp)s_1>0.5 && deepTauVsJets_medium_2>0.5 && leptonveto==0 && (trg_doubletau) && deepTauVsEle_vvloose_1 && deepTauVsEle_vvloose_2 && deepTauVsMu_vloose_1 && deepTauVsMu_vloose_2)' % vars()


elif options.channel == 'em':
    cats['baseline'] = '(iso_1<0.15 && iso_2<0.2 && !leptonveto)'
    if options.era == 'mssmsummer16':
      cats['loose_baseline'] = '(iso_1<0.5 && iso_2>0.2 && iso_2<0.5 && !leptonveto &&trg_muonelectron)'
    elif options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']:
      cats['loose_baseline'] = '(wt<2 && iso_1<0.15 && iso_2<0.5 && !leptonveto &&trg_muonelectron)'
    if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']: cats['baseline'] = '(iso_1<0.15 && iso_2<0.2 && !leptonveto && trg_muonelectron && pt_1>15 && pt_2>15)'
    if options.era in ['cpsummer17','UL_17','cp18','UL_18']:
        cats['baseline'] = '(iso_1<0.15 && iso_2<0.2 && !leptonveto && trg_muonelectron && pt_1>15 && pt_2>15)'
        cats['loose_baseline'] = '(wt<2 && iso_1<0.15 && iso_2<0.5 && !leptonveto && trg_muonelectron)'
elif options.channel == 'zmm':
    cats['baseline'] = '(iso_1<0.15 && iso_2<0.15)'
    if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']: cats['baseline'] = '(iso_1<0.15 && iso_2<0.15 && trg_singlemuon)'
    if options.era in ['cpsummer17','UL_17','cp18','UL_18']: cats['baseline'] = '(pt_1>25 && iso_1<0.15 && iso_2<0.15 && trg_singlemuon)'
elif options.channel == 'zee':
    cats['baseline'] = '(iso_1<0.1 && iso_2<0.1)'
    if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']: cats['baseline'] = '(iso_1<0.1 && iso_2<0.1 && trg_singleelectron && fabs(wt)<2)'

    if options.era in ['cpsummer17','UL_17']: cats['baseline'] = '(pt_1>28&&pt_2>13&&iso_1<0.15 && iso_2<0.15 && trg_singleelectron && fabs(wt)<2)'
    if options.era in ['cp18','UL_18']: 
        cats['baseline'] = '(pt_1>33 && pt_2>13 && iso_1<0.15 && iso_2<0.15 && trg_singleelectron && fabs(wt)<2)'

if options.analysis == 'cpdecay':
  if options.channel in ['mt','et']: cats['baseline'] += ' && mva_dm_2>=0 && (mva_dm_2>=1&&tau_decay_mode_2==0)==0 && m_vis>40'
  if options.channel in ['tt']: cats['baseline'] += ' && mva_dm_1>=0 && mva_dm_2>=0 && (mva_dm_1>=1&&tau_decay_mode_1==0)==0 && (mva_dm_2>=1&&tau_decay_mode_2==0)==0 && mva_dm_1<11 && mva_dm_2<11 && m_vis>40 && (mva_dm_1!=0 || ip_sig_1>1.5) && (mva_dm_2!=0 || ip_sig_2>1.5)'

  cats['tt_loose_baseline'] = '((deepTauVsJets_medium_1>0.5 && deepTauVsJets_vvvloose_2>0.5 && deepTauVsJets_medium_2<0.5 && leptonveto==0 && (trg_doubletau && pt_2>40) && deepTauVsEle_vvloose_1 && deepTauVsEle_vvloose_2 && deepTauVsMu_vloose_1 && deepTauVsMu_vloose_2) && mva_dm_1>=0 && mva_dm_2>=0 && (mva_dm_1>=1&&tau_decay_mode_1==0)==0 && (mva_dm_2>=1&&tau_decay_mode_2==0)==0 && mva_dm_1<11 && mva_dm_2<11 && m_vis>40 && (mva_dm_1!=0 || ip_sig_1>1.5) && (mva_dm_2!=0 || ip_sig_2>1.5))'
  cats['tt_loose_baseline_2'] = '((deepTauVsJets_medium_1>0.5 && deepTauVsJets_vvvloose_2>0.5 && deepTauVsJets_vvloose_2<0.5 && leptonveto==0 && (trg_doubletau && pt_2>40) && deepTauVsEle_vvloose_1 && deepTauVsEle_vvloose_2 && deepTauVsMu_vloose_1 && deepTauVsMu_vloose_2) && mva_dm_1>=0 && mva_dm_2>=0 && (mva_dm_1>=1&&tau_decay_mode_1==0)==0 && (mva_dm_2>=1&&tau_decay_mode_2==0)==0 && mva_dm_1<11 && mva_dm_2<11 && m_vis>40 && (mva_dm_1!=0 || ip_sig_1>1.5) && (mva_dm_2!=0 || ip_sig_2>1.5))'


if options.w_ff_closure:
  cats['baseline'] = '(' + cats['baseline'] + ' && mt_1>70 && n_deepbjets==0)'
elif options.qcd_ff_closure:
  if options.channel in ['et','mt']:
    cats['baseline'] = '(' + cats['baseline'] + ' && mt_1<50 && iso_1>0.05)'
    options.do_ss = True
  elif options.channel == 'tt':
    options.do_ss = True


if options.analysis == 'cpprod':
  if options.channel in ['mt','et']: cats['baseline'] += ' && pt_2>30'
  if options.channel in ['tt']: cats['baseline'] += ' && pt_1>50'
  cats['baseline'] = cats['baseline'].replace('deepTauVsEle_vvloose_2','deepTauVsEle_vvvloose_2').replace('deepTauVsEle_vvloose_1','deepTauVsEle_vvvloose_1')

cats['inclusive'] = '(1)' 
cats['w_os'] = 'os'
cats['w_sdb'] = 'mt_1>70.'
cats['pass'] = 'mva_olddm_tight_2>0.5 && pzeta>-25'
cats['fail'] = 'mva_olddm_tight_2<0.5 && pzeta>-25'
if options.era in ['smsummer16']: cats['w_sdb'] = 'mt_1>80.'
cats['w_sdb_os'] = 'os'
cats['tt_qcd_norm'] = '(mva_olddm_tight_1>0.5 && mva_olddm_medium_2>0.5 &&mva_olddm_tight_2<0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto)&&trg_doubletau'
if options.era == 'mssmsummer16': cats['tt_qcd_norm'] = '(mva_olddm_medium_1>0.5 && mva_olddm_loose_2>0.5 &&mva_olddm_medium_2<0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto)&&trg_doubletau'
if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']: cats['tt_qcd_norm'] = '(pt_1>40 && ((mva_olddm_loose_1>0.5 && mva_olddm_tight_1<0.5 && mva_olddm_medium_2>0.5) || (mva_olddm_loose_2>0.5 && mva_olddm_tight_2<0.5 && mva_olddm_medium_1>0.5))  && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto)&&trg_doubletau'
# if options.era in ['cpsummer17','cp18']: cats['tt_qcd_norm'] = '(mva_olddm_tight_1>0.5 && mva_olddm_tight_2<0.5 && mva_olddm_medium_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto && trg_doubletau)'
cats['qcd_loose_shape'] = '(iso_1>0.2 && iso_1<0.5 && mva_olddm_tight_2>0.5 && antiele_2 && antimu_2 && !leptonveto)'

cats['tt_qcd_norm'] = '(deepTauVsJets_medium_1<0.5 && deepTauVsJets_loose_1>0.5 && deepTauVsJets_medium_2>0.5 && leptonveto==0 && (trg_doubletau && pt_2>40) && deepTauVsEle_vvloose_1 && deepTauVsEle_vvloose_2 && deepTauVsMu_vloose_1 && deepTauVsMu_vloose_2)'


# CR categories
cats['ztt_control'] = '(m_sv>60&&m_sv<100)'
cats['ztt_control_dijet'] = '(m_sv>60&&m_sv<100 && n_jets>1)'
if options.channel == 'em':
  cats['ztt_control'] = '(m_sv>60&&m_sv<100 && n_bjets==0)'
  cats['ztt_control_dijet'] = '(m_sv>60&&m_sv<100 && n_jets>1 && n_bjets==0)'

# MSSM categories
#cats['btag'] = '(n_bjets>=1)'
#cats['nobtag'] = '(n_bjets==0)'
# loose/tight iso-MT categories
cats['atleast1bjet'] = '(n_bjets>0)'
cats['btag_wnobtag']='(n_lowpt_jets>=1)' # this is the one that is used for the b-tag method 16!
cats['0jet'] = '(n_jets==0)'
cats['1jet'] = '(n_jets==1)'
cats['ge2jet'] = '(n_jets>=2)'
cats['btag_tight_wnobtag']='(n_lowpt_jets>=1)'
cats['w_shape']=''
cats['qcd_shape']=''
cats['w_shape_comp']=''
cats['qcd_shape_comp']=''

# MSSM run 2 categories
if options.channel == "tt":
  cats['NbtagGt1'] = '(n_deepbjets>0)'
  cats['Nbtag0'] = '(n_deepbjets==0)'
elif options.channel in ["et","mt"]:
  cats['NbtagGt1'] = '(n_deepbjets>0 && mt_1<70)'
  cats['Nbtag0'] = '(n_deepbjets==0 && mt_1<70)'
elif options.channel == "em":
  cats['NbtagGt1'] = '(n_deepbjets>0 && pzeta>-35)'
  cats['Nbtag0'] = '(n_deepbjets==0 && pzeta>-35)'




cats['Nbtag0_MTLt40'] = '(n_deepbjets==0 && mt_1<40)'
cats['Nbtag0_MT40To70'] = '(n_deepbjets==0 && mt_1>40 && mt_1<70)'
cats['NbtagGt1_MTLt40'] = '(n_deepbjets>0 && mt_1<40)'
cats['NbtagGt1_MT40To70'] = '(n_deepbjets>0 && mt_1>40 && mt_1<70)'

cats['Nbtag0_MHGt250'] = '(n_deepbjets==0 && svfit_mass>250)'
cats['Nbtag0_MHGt200'] = '(n_deepbjets==0 && svfit_mass>200)'

cats['Nbtag0_MTLt40_MHGt250'] = '(n_deepbjets==0 && mt_1<40 && svfit_mass>250)'
cats['Nbtag0_MT40To70_MHGt250'] = '(n_deepbjets==0 && mt_1>40 && mt_1<70 && svfit_mass>250)'

cats['Nbtag0_MTLt40_MHGt200'] = '(n_deepbjets==0 && mt_1<40 && svfit_mass>200)'
cats['Nbtag0_MT40To70_MHGt200'] = '(n_deepbjets==0 && mt_1>40 && mt_1<70 && svfit_mass>200)'

cats['Nbtag1'] = '(n_deepbjets==1)'
cats['Nbtag1_MTLt40'] = '(n_deepbjets==1 && mt_1<40)'
cats['Nbtag1_MT40To70'] = '(n_deepbjets==1 && mt_1>40 && mt_1<70)'

cats['Nbtag0_MTLt70'] = '(n_deepbjets==0 && mt_1<70)'
cats['NbtagGt1_MTLt70'] = '(n_deepbjets>0 && mt_1<70)'
cats['MTLt70'] = '(mt_1<70)'

cats['tightmt'] = '(mt_1<40)'
cats['loosemt'] = '(mt_1>40 && mt_1<70)'

cats['wjets_control'] = '(mt_1>70 && n_deepbjets==0)'
if options.channel == "tt":
  cats['qcd_control'] = '(1)'
elif options.channel in ["et","mt"]:
  cats['qcd_control'] = 'mt_1<50 && iso_1>0.05'

if options.cat == 'qcd_control':
  options.do_ss = True

cats['Nbtag0_DZetaGt30'] = '(n_deepbjets==0 && pzeta>30)'
cats['Nbtag0_DZetam10To30'] = '(n_deepbjets==0 && pzeta<=30 && pzeta>-10)'
cats['Nbtag0_DZetam35Tom10'] = '(n_deepbjets==0 && pzeta<=-10 && pzeta>-35)'
cats['NbtagGt1_DZetaGt30'] ='(n_deepbjets>0 && pzeta>30)'
cats['NbtagGt1_DZetam10To30'] ='(n_deepbjets>0 && pzeta<=30 && pzeta>-10)'
cats['NbtagGt1_DZetam35Tom10'] ='(n_deepbjets>0 && pzeta<=-10 && pzeta>-35)'
cats['NbtagGt1_DZetaLtm35'] = '(n_deepbjets>0 && pzeta<=-35)'

cats['Nbtag0_DZetaGt30_MHGt250'] = '(n_deepbjets==0 && pzeta>30 && svfit_mass>250)'
cats['Nbtag0_DZetam10To30_MHGt250'] = '(n_deepbjets==0 && pzeta<=30 && pzeta>-10 && svfit_mass>250)'
cats['Nbtag0_DZetam35Tom10_MHGt250'] = '(n_deepbjets==0 && pzeta<=-10 && pzeta>-35 && svfit_mass>250)'

cats['DZetaGtm35'] = '(pzeta>-35)'


if options.channel == 'et': cats['baseline_loose'] = '(iso_1<0.3 && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto &&  trg_singleelectron)'
if options.channel == 'mt': cats['baseline_loose'] = '(iso_1<0.3 && mva_olddm_medium_2>0.5 && antiele_2 && antimu_2 && !leptonveto && (trg_singlemuon*(pt_1>23) || trg_mutaucross*(pt_1<23)))'

# MVA DM categories for SF measurments
cats['mvadm_pi'] = 'mvadm_pi_2>mvadm_rho_2&&mvadm_pi_2>mvadm_a1_2'
cats['mvadm_rho'] = 'mvadm_rho_2>mvadm_pi_2&&mvadm_rho_2>mvadm_a1_2'
cats['mvadm_a1'] = 'mvadm_a1_2>mvadm_pi_2&&mvadm_a1_2>mvadm_rho_2'
cats['mvadm_rho'] = 'mvadm_rho_2>mvadm_pi_2&&mvadm_rho_2>mvadm_a1_2'
cats['mvadm_notrho'] = '!(mvadm_rho_2>mvadm_pi_2&&mvadm_rho_2>mvadm_a1_2)'

# SM categories
cats['0jet'] = '(n_jets==0)'
if options.channel == 'em': cats['0jet'] = '(n_jets==0 && n_bjets==0)'
cats['vbf'] = '(0)'
cats['twojets'] = '(n_jets>=2)'
if options.channel == 'em': cats['vbf'] = '(n_jets==2 && mjj>300 && n_bjets==0)'
if options.channel == 'et': cats['vbf'] = '(n_jets>=2 && mjj>300 && pt_tt>50)'
if options.channel == 'mt': cats['vbf'] = '(n_jets>=2 && mjj>300 && pt_tt>50 && pt_2>40)'
if options.channel == 'tt': cats['vbf'] = '(n_jets>=2 && pt_tt>100 && jdeta>2.5)'
cats['boosted'] = '(!(%s) && !(%s))' % (cats['0jet'], cats['vbf'])
if options.channel == 'em': cats['boosted'] = '(!(%s) && !(%s) && n_bjets==0)' % (cats['0jet'], cats['vbf'])

if options.era in ['cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','cpsummer17','UL_17','cp18','UL_18','mvadm2016'] or options.analysis in ['cpprod']:
  cats['ttbar'] = 'n_jets>0'  
  if options.channel in ['mt','et']: 
      cats['0jet'] = '(n_jets==0 && n_bjets==0 &&n_loose_bjets<2)'
      cats['dijet']='n_jets>=2 && mjj>300 && n_bjets==0 &&n_loose_bjets<2'
      cats['dijet_boosted']='%s && pt_tt>200' % cats['dijet']
      cats['dijet_lowboost']='%s && pt_tt<200' % cats['dijet']
      cats['dijet_loosemjj_boosted']='%s && mjj<500 && pt_tt>150' % cats['dijet']
      cats['dijet_loosemjj_lowboost']='%s && mjj<500 && pt_tt<150' % cats['dijet']
      cats['dijet_tightmjj_boosted']='%s && mjj>500 && pt_tt>150' % cats['dijet']
      cats['dijet_tightmjj_lowboost']='%s && mjj>500 && pt_tt<150' % cats['dijet']
      cats['boosted'] = '(!(%s) && !(%s) && n_bjets==0 &&n_loose_bjets<2)' % (cats['0jet'], cats['dijet'])
      #cats['boosted'] = '(n_jets==1&& n_bjets==0)'
  elif options.channel in ['em']:
      cats['0jet'] = '(n_jets==0 && n_loose_bjets<2 && n_bjets==0)'
      cats['dijet']='n_jets>=2 && mjj>300 && n_loose_bjets<2 && n_bjets==0 '
      cats['dijet_boosted']='%s && pt_tt>200' % cats['dijet']
      cats['dijet_lowboost']='%s && pt_tt<200' % cats['dijet']
      cats['dijet_loosemjj_boosted']='%s && mjj<500 && pt_tt>150' % cats['dijet']
      cats['dijet_loosemjj_lowboost']='%s && mjj<500 && pt_tt<150' % cats['dijet']
      cats['dijet_tightmjj_boosted']='%s && mjj>500 && pt_tt>150' % cats['dijet']
      cats['dijet_tightmjj_lowboost']='%s && mjj>500 && pt_tt<150' % cats['dijet']
      cats['boosted'] = '(!(%s) && !(%s) &&n_loose_bjets<2 && n_bjets==0)' % (cats['0jet'], cats['dijet'])
      #cats['boosted'] = '(n_jets==1&& n_bjets==0)'
  else:    
    cats['0jet'] = '(n_jets==0)'
    cats['dijet']='n_jets>=2 && mjj>300'
    cats['dijet_boosted']='%s && pt_tt>200' % cats['dijet']
    cats['dijet_lowboost']='%s && pt_tt<200' % cats['dijet']
    cats['dijet_loosemjj_boosted']='%s && mjj<500 && pt_tt>150' % cats['dijet']
    cats['dijet_loosemjj_lowboost']='%s && mjj<500 && pt_tt<150' % cats['dijet']
    cats['dijet_tightmjj_boosted']='%s && mjj>500 && pt_tt>150' % cats['dijet']
    cats['dijet_tightmjj_lowboost']='%s && mjj>500 && pt_tt<150' % cats['dijet']

    cats['boosted'] = '(!(%s) && !(%s))' % (cats['0jet'], cats['dijet'])
    #cats['boosted'] = '(n_jets==1)'

#cats['0jet'] = '(n_jets==0 && n_bjets==0)'
cats['1jet'] = '(n_jets==1 && n_bjets==0)'
cats['2jet'] = '(n_jets>=2 && n_bjets==0)'

# SM ML categories

if options.channel == 'tt':

    mva_highMjj_ggh =      '(IC_highMjj_Oct05_max_index==0)'
    mva_highMjj_jetFakes = '(IC_highMjj_Oct05_max_index==1)'
    mva_highMjj_zttEmbed = '(IC_highMjj_Oct05_max_index==2)'

    mva_lowMjj_ggh =      '(IC_lowMjj_Oct05_max_index==0)'
    mva_lowMjj_jetFakes = '(IC_lowMjj_Oct05_max_index==1)'
    mva_lowMjj_misc =     '(IC_lowMjj_Oct05_max_index==2)'
    mva_lowMjj_qqh =      '(IC_lowMjj_Oct05_max_index==2)'
    mva_lowMjj_zttEmbed = '(IC_lowMjj_Oct05_max_index==3)'

    # with simple mjj and pt_tt cuts like in Cutbased analysis
    # cats['ggh_lowboost_highMjj'] =  '({} && n_jets>=2 && mjj>300 && mjj<500 && pt_tt<150)'.format(mva_highMjj_ggh)
    # cats['ggh_boosted_highMjj'] =  '({} && n_jets>=2 && mjj>300 && mjj<500 && pt_tt>150)'.format(mva_highMjj_ggh)
    # cats['qqh_lowboost_highMjj'] =  '({} && n_jets>=2 && mjj>500 && pt_tt<150)'.format(mva_highMjj_ggh)
    # cats['qqh_boosted_highMjj'] =  '({} && n_jets>=2 && mjj>500 && pt_tt>150)'.format(mva_highMjj_ggh)
    
    # with binary ggH vs qqH mva
    # cats['ggh_highMjj'] = '({} && n_jets>=2 && mjj>300 && IC_binary_Oct11_score>0.7)'.format(mva_highMjj_ggh)
    # cats['ggh_lowboost_highMjj'] = '({} && n_jets>=2 && pt_tt<150 && mjj>300)'.format(mva_highMjj_ggh)
    # cats['ggh_boosted_highMjj'] = '({} && n_jets>=2 && pt_tt>150 && mjj>300)'.format(mva_highMjj_ggh)
    # cats['ggh_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt<150 && IC_binary_Oct11_score>0.5)'.format(mva_highMjj_ggh)
    # cats['ggh_boosted_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt>150 && IC_binary_Oct11_score>0.5)'.format(mva_highMjj_ggh)

    # simplest mjj based qqH vs ggH
    cats['ggh_highMjj'] =  '({} && n_jets>=2 && mjj>300 && mjj<500)'.format(mva_highMjj_ggh)
    cats['qqh_highMjj'] =  '({} && n_jets>=2 && mjj>500)'.format(mva_highMjj_ggh)
    cats['zttEmbed_highMjj'] =  '({} && n_jets>=2 && mjj>300)'.format(mva_highMjj_zttEmbed)
    cats['jetFakes_highMjj'] =  '({} && n_jets>=2 && mjj>300)'.format(mva_highMjj_jetFakes)

    cats['ggh_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_ggh)
    cats['misc_lowMjj'] = '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_misc)
    cats['qqh_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_qqh)
    cats['zttEmbed_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_zttEmbed)
    cats['jetFakes_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_jetFakes)

if options.channel in ['et','mt']:

    mva_highMjj_ggh =      '(IC_highMjj_Oct05_max_index==0)'
    mva_highMjj_jetFakes = '(IC_highMjj_Oct05_max_index==1)'
    mva_highMjj_tt =       '(IC_highMjj_Oct05_max_index==2)'
    if options.channel == 'et':
        mva_highMjj_zll      = '(IC_highMjj_Oct05_max_index==3)'
        mva_highMjj_zttEmbed = '(IC_highMjj_Oct05_max_index==4)'
    else:
        mva_highMjj_zttEmbed = '(IC_highMjj_Oct05_max_index==3)'

    mva_lowMjj_ggh =      '(IC_lowMjj_Oct05_max_index==0)'
    mva_lowMjj_jetFakes = '(IC_lowMjj_Oct05_max_index==1)'
    mva_lowMjj_qqh =      '(IC_lowMjj_Oct05_max_index==2)'
    mva_lowMjj_tt =       '(IC_lowMjj_Oct05_max_index==3)'
    mva_lowMjj_zll =      '(IC_lowMjj_Oct05_max_index==4)'
    mva_lowMjj_zttEmbed = '(IC_lowMjj_Oct05_max_index==5)'

    # if options.channel == 'et':
        # latest one with binary MVA
        # cats['ggh_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt<150 && IC_binary_Oct11_score>0.5)'.format(mva_highMjj_ggh)
        # cats['ggh_boosted_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt>150 && IC_binary_Oct11_score>0.5)'.format(mva_highMjj_ggh)
        # cats['qqh_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt<150 && IC_binary_Oct11_score<0.5)'.format(mva_highMjj_ggh)
        # cats['qqh_boosted_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt>150 && IC_binary_Oct11_score<0.5)'.format(mva_highMjj_ggh)

    # else:
        # latest ones with binary MVA
        # cats['ggh_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt<150 &&IC_binary_Oct11_score>0.45)'.format(mva_highMjj_ggh)
        # cats['ggh_boosted_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt>150 && IC_binary_Oct11_score>0.45)'.format(mva_highMjj_ggh)
        # cats['qqh_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt<150 && IC_binary_Oct11_score<0.45)'.format(mva_highMjj_ggh)
        # cats['qqh_boosted_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt>150 && IC_binary_Oct11_score<0.45)'.format(mva_highMjj_ggh)

    cats['ggh_highMjj'] = '({} && n_jets>=2 && mjj>300 && mjj<500)'.format(mva_highMjj_ggh)
    cats['qqh_highMjj'] = '({} && n_jets>=2 && mjj>500)'.format(mva_highMjj_ggh)
    cats['tt_highMjj'] =   '({} && n_jets>=2 && mjj>300)'.format(mva_highMjj_tt)
    if options.channel == 'et':
        cats['zll_highMjj'] =  '({} && n_jets>=2 && mjj>300)'.format(mva_highMjj_zll)
    cats['zttEmbed_highMjj'] =  '({} && n_jets>=2 && mjj>300)'.format(mva_highMjj_zttEmbed)
    cats['jetFakes_highMjj'] =  '({} && n_jets>=2 && mjj>300)'.format(mva_highMjj_jetFakes)

    cats['jetFakes_lowMjj'] = '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_jetFakes)
    cats['ggh_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_ggh)
    cats['qqh_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_qqh)
    cats['tt_lowMjj'] =   '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_tt)
    cats['zll_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_zll)
    cats['zttEmbed_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_zttEmbed)

if options.channel == 'em':

    mva_highMjj_ggh =      '(IC_highMjj_Oct05_max_index==0)'
    mva_highMjj_tt =       '(IC_highMjj_Oct05_max_index==1)'
    mva_highMjj_zttEmbed = '(IC_highMjj_Oct05_max_index==2)'

    mva_lowMjj_ggh =      '(IC_lowMjj_Oct05_max_index==0)'
    mva_lowMjj_qcd =      '(IC_lowMjj_Oct05_max_index==1)'
    mva_lowMjj_qqh =      '(IC_lowMjj_Oct05_max_index==2)'
    mva_lowMjj_tt =       '(IC_lowMjj_Oct05_max_index==3)'
    mva_lowMjj_zttEmbed = '(IC_lowMjj_Oct05_max_index==4)'

    # latest one with binary MVA
    # cats['ggh_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt<150 && IC_binary_Oct11_score>0.55)'.format(mva_highMjj_ggh)
    # cats['ggh_boosted_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt>150 && IC_binary_Oct11_score>0.55)'.format(mva_highMjj_ggh)
    # cats['qqh_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt<150 && IC_binary_Oct11_score<0.55)'.format(mva_highMjj_ggh)
    # cats['qqh_boosted_highMjj'] = '({} && n_jets>=2 && mjj>300 && pt_tt>150 && IC_binary_Oct11_score<0.55)'.format(mva_highMjj_ggh)

    cats['ggh_highMjj'] = '({} && n_jets>=2 && mjj>300 && mjj<500)'.format(mva_highMjj_ggh)
    cats['qqh_highMjj'] = '({} && n_jets>=2 && mjj>500)'.format(mva_highMjj_ggh)
    cats['tt_highMjj'] =   '({} && n_jets>=2 && mjj>300)'.format(mva_highMjj_tt)
    cats['zttEmbed_highMjj'] =  '({} && n_jets>=2 && mjj>300)'.format(mva_highMjj_zttEmbed)

    cats['ggh_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_ggh)
    cats['qcd_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_qcd)
    cats['qqh_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_qqh)
    cats['tt_lowMjj'] =   '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_tt)
    cats['zttEmbed_lowMjj'] =  '({} && !(n_jets>=2 && mjj>300))'.format(mva_lowMjj_zttEmbed)

# CP in decays categories
if options.channel in ["mt","et"]:
    cats["inclusive_mixed"]         = "(tau_decay_mode_2==1)"
    cats["dijet_mixed"]             = "(n_jets>=2 && mjj>300)"
    cats["0jet_mixed"]              = "({} && {})".format(cats["0jet"], cats["inclusive_mixed"])
    cats["dijet_boosted_mixed"]     = "({} && pt_tt>100 && {})".format(cats["dijet_mixed"] ,cats["inclusive_mixed"])
    cats["dijet_lowboost_mixed"]    = "({} && pt_tt<100 && {})".format(cats["dijet_mixed"] ,cats["inclusive_mixed"])
    cats["boosted_mixed"]               = "(!{} && !({}) && {})".format(cats["0jet"], cats["dijet_mixed"], cats["inclusive_mixed"])

    # MVA multiclass categories
    cut_string = ""
    if options.channel == "et":
        cut_string = "IC_Apr02"
    elif options.channel == "mt":
        cut_string = "IC_Nov25_tauspinner"
    mva_ggh      = '({}_max_index==0)'.format(cut_string)
    mva_jetFakes = '({}_max_index==1)'.format(cut_string)
    mva_zttEmbed = '({}_max_index==2)'.format(cut_string)

    # apply b jet veto as done in training (to remove ttbar)
    cats['higgs']     = '({} && n_bjets==0)'.format(mva_ggh)
    cats['zttEmbed']  = '({} && n_bjets==0)'.format(mva_zttEmbed)
    cats['jetFakes']  = '({} && n_bjets==0)'.format(mva_jetFakes)

    cats["mva_mupi"]  = "(tau_decay_mode_2==0 && mva_dm_2==0)" # aco_angle_6 (ip)
    cats["mva_murho"] = "(tau_decay_mode_2==1 && mva_dm_2==1)" # aco_angle_5 (planes) 6 (ip)
    cats["mva_mua1"]  = "(tau_decay_mode_2>=10 && mva_dm_2==10)" # aco_angle_5 (planes)

    cats["higgs_mvaMuPi"]  = "{} && {} && n_bjets==0".format(cats["higgs"],cats["mva_mupi"])
    cats["higgs_mvaMuRho"] = "{} && {} && n_bjets==0".format(cats["higgs"],cats["mva_murho"])
    cats["higgs_mvaMuA1"]  = "{} && {} && n_bjets==0".format(cats["higgs"],cats["mva_mua1"])

    cats["inclusive_mupi"]  = "(tau_decay_mode_2==0)"
    cats["inclusive_murho"] = "(tau_decay_mode_2==1)"

if options.channel == 'tt':

# start fresh here
    cats["inclusive_rhorho"]       = "(tau_decay_mode_1==1 && tau_decay_mode_2==1)"
    cats["inclusive_pirho"]       = "((tau_decay_mode_1==1 && tau_decay_mode_2==0 && ip_sig_2>=1.5) || (tau_decay_mode_1==0 && ip_sig_1>=1.5 && tau_decay_mode_2==1))"
    cats["inclusive_a1rho"]     = "((tau_decay_mode_1==10 && tau_decay_mode_2==1) || (tau_decay_mode_1==1 && tau_decay_mode_2==10))"
    cats["inclusive_a1pi"]     = "((tau_decay_mode_1==10 && tau_decay_mode_2==0 && ip_sig_2>=1.5) || (tau_decay_mode_1==0 && ip_sig_1>=1.5 && tau_decay_mode_2==10))"
    cats["inclusive_a1a1"]     = "(tau_decay_mode_1==10 && tau_decay_mode_2==10)"
    cats["inclusive_pipi"]     = "(tau_decay_mode_1==0 && ip_sig_1>=1.5 && tau_decay_mode_2==0 && ip_sig_2>=1.5)"


    cats["inclusive_mvapirho_up"]       = "((tau_decay_mode_1==1 && mva_dm_1==1 && ip_sig_2_up>=1.5 && mva_dm_2==0) || (ip_sig_1_up>=1.5 && mva_dm_1==0 && tau_decay_mode_2==1 && mva_dm_2==1))"
    cats["inclusive_mvaa1pi_up"]     = "((mva_dm_1==10 && ip_sig_2_up>=1.5 && mva_dm_2==0) || (ip_sig_1_up>=1.5 && mva_dm_1==0 && mva_dm_2==10))"
    cats["inclusive_mvapipi_up"]     = "(mva_dm_1==0 && ip_sig_1_up>=1.5 && ip_sig_2_up>=1.5 && mva_dm_2==0)"
    cats["inclusive_mvapi0a1_up"]     = "((mva_dm_1==0 && ip_sig_1_up>=1.5 && tau_decay_mode_2==1 && mva_dm_2==2) || (tau_decay_mode_1==1 && mva_dm_1==2 && mva_dm_2==0 && ip_sig_2_up>=1.5))"

    cats["inclusive_mvapirho_down"]       = "((tau_decay_mode_1==1 && mva_dm_1==1 && ip_sig_2_down>=1.5 && mva_dm_2==0) || (ip_sig_1_down>=1.5 && mva_dm_1==0 && tau_decay_mode_2==1 && mva_dm_2==1))"
    cats["inclusive_mvaa1pi_down"]     = "((mva_dm_1==10 && ip_sig_2_down>=1.5 && mva_dm_2==0) || (ip_sig_1_down>=1.5 && mva_dm_1==0 && mva_dm_2==10))"
    cats["inclusive_mvapipi_down"]     = "(mva_dm_1==0 && ip_sig_1_down>=1.5 && ip_sig_2_down>=1.5 && mva_dm_2==0)"
    cats["inclusive_mvapi0a1_down"]     = "((mva_dm_1==0 && ip_sig_1_down>=1.5 && tau_decay_mode_2==1 && mva_dm_2==2) || (tau_decay_mode_1==1 && mva_dm_1==2 && mva_dm_2==0 && ip_sig_2_down>=1.5))"

    cats["inclusive_mvarhorho"]       = "(tau_decay_mode_1==1 && mva_dm_1==1 && tau_decay_mode_2==1 && mva_dm_2==1)"
    cats["inclusive_mvapirho"]       = "((tau_decay_mode_1==1 && mva_dm_1==1 && ip_sig_2>=1.5 && mva_dm_2==0) || (ip_sig_1>=1.5 && mva_dm_1==0 && tau_decay_mode_2==1 && mva_dm_2==1))"
    cats["inclusive_mvaa1rho"]     = "((mva_dm_1==10 && tau_decay_mode_2==1 && mva_dm_2==1) || (tau_decay_mode_1==1 && mva_dm_1==1 && mva_dm_2==10))"
    cats["inclusive_mvaa1pi"]     = "((mva_dm_1==10 && ip_sig_2>=1.5 && mva_dm_2==0) || (ip_sig_1>=1.5 && mva_dm_1==0 && mva_dm_2==10))"
    cats["inclusive_mvaa1a1"]     = "(mva_dm_1==10 && mva_dm_2==10)"
    cats["inclusive_mvapipi"]     = "(mva_dm_1==0 && ip_sig_1>=1.5 && ip_sig_2>=1.5 && mva_dm_2==0)"
    cats["inclusive_mvapi0a1"]     = "((mva_dm_1==0 && ip_sig_1>=1.5 && tau_decay_mode_2==1 && mva_dm_2==2) || (tau_decay_mode_1==1 && mva_dm_1==2 && mva_dm_2==0 && ip_sig_2>=1.5))"
    cats["inclusive_mvarho0a1"]     = "(tau_decay_mode_1==1 && tau_decay_mode_2==1 && ((mva_dm_1==1&&mva_dm_2==2) || (mva_dm_1==2&&mva_dm_2==1) || (mva_dm_1==2&&mva_dm_2==2)))"
    cats["inclusive_mvaa10a1"]     = "((mva_dm_1==10 && tau_decay_mode_2==1 && mva_dm_2==2) || (tau_decay_mode_1==1 && mva_dm_1==2 && mva_dm_2==10))"

    #mva_ggh                = '(svfit_mass>100&&svfit_mass<150)'
    #mva_jetFakes           = '(svfit_mass>150)'
    #mva_zttEmbed           = '(svfit_mass<100)'

    cut_string = "IC_01Jun2020"
    mva_ggh      = '({}_max_index==0)'.format(cut_string)
    mva_jetFakes = '({}_max_index==1)'.format(cut_string)
    mva_zttEmbed = '({}_max_index==2)'.format(cut_string)


    cats['higgs_mvarhorho']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvarhorho"])
    cats['zttEmbed_mvarhorho']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvarhorho"])
    cats['jetFakes_mvarhorho']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvarhorho"])
    cats['higgs_mvapirho']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapirho"])
    cats['zttEmbed_mvapirho']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvapirho"])
    cats['jetFakes_mvapirho']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvapirho"])
    cats['higgs_mvaa1rho']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvaa1rho"])
    cats['zttEmbed_mvaa1rho']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvaa1rho"])
    cats['jetFakes_mvaa1rho']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvaa1rho"])
    cats['higgs_mvaa1pi']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvaa1pi"])
    cats['zttEmbed_mvaa1pi']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvaa1pi"])
    cats['jetFakes_mvaa1pi']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvaa1pi"])
    cats['higgs_mvaa1a1']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvaa1a1"])
    cats['higgs_mvaa1a1_sync']      = '({} && {} && tau_decay_mode_1==10 && tau_decay_mode_2==10 && a1_flag)'.format(mva_ggh, cats["inclusive_mvaa1a1"])
    cats['zttEmbed_mvaa1a1']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvaa1a1"])
    cats['jetFakes_mvaa1a1']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvaa1a1"])
    cats['higgs_mvapipi']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapipi"])
    cats['zttEmbed_mvapipi']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvapipi"])
    cats['jetFakes_mvapipi']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvapipi"])
    cats['higgs_mvapi0a1']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapi0a1"])
    cats['zttEmbed_mvapi0a1']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvapi0a1"])
    cats['jetFakes_mvapi0a1']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvapi0a1"])
    cats['higgs_mvarho0a1']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvarho0a1"])
    cats['zttEmbed_mvarho0a1']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvarho0a1"])
    cats['jetFakes_mvarho0a1']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvarho0a1"])


    cats['higgs_mvaa10a1']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvaa10a1"])

    cats['zttEmbed_mvaa10a1']      = '({} && {})'.format(mva_zttEmbed, cats["inclusive_mvaa10a1"])
    cats['jetFakes_mvaa10a1']      = '({} && {})'.format(mva_jetFakes, cats["inclusive_mvaa10a1"])

    cats['higgs']      = '({} && (mva_dm_1==0&&ip_sig_1<1.5)==0 && (mva_dm_2==0&&ip_sig_2<1.5)==0)'.format(mva_ggh)
    cats['zttEmbed']   = '({} && (mva_dm_1==0&&ip_sig_1<1.5)==0 && (mva_dm_2==0&&ip_sig_2<1.5)==0)'.format(mva_zttEmbed)
    cats['jetFakes']   = '({} && (mva_dm_1==0&&ip_sig_1<1.5)==0 && (mva_dm_2==0&&ip_sig_2<1.5)==0)'.format(mva_jetFakes)

    cats['higgs_up']      = '({} && (mva_dm_1==0&&ip_sig_1_up<1.5)==0 && (mva_dm_2==0&&ip_sig_2_up<1.5)==0)'.format(mva_ggh)
    cats['zttEmbed_up']   = '({} && (mva_dm_1==0&&ip_sig_1_up<1.5)==0 && (mva_dm_2==0&&ip_sig_2_up<1.5)==0)'.format(mva_zttEmbed)
    cats['jetFakes_up']   = '({} && (mva_dm_1==0&&ip_sig_1_up<1.5)==0 && (mva_dm_2==0&&ip_sig_2_up<1.5)==0)'.format(mva_jetFakes)

    cats['higgs_down']      = '({} && (mva_dm_1==0&&ip_sig_1_down<1.5)==0 && (mva_dm_2==0&&ip_sig_2_down<1.5)==0)'.format(mva_ggh)
    cats['zttEmbed_down']   = '({} && (mva_dm_1==0&&ip_sig_1_down<1.5)==0 && (mva_dm_2==0&&ip_sig_2_down<1.5)==0)'.format(mva_zttEmbed)
    cats['jetFakes_down']   = '({} && (mva_dm_1==0&&ip_sig_1_down<1.5)==0 && (mva_dm_2==0&&ip_sig_2_down<1.5)==0)'.format(mva_jetFakes)

    cats['higgs_mvaother']    = '({} && !({}||{}||{}||{}||{}||{}||{}||{}||{}))'\
            .format(mva_ggh, cats["higgs_mvarhorho"], cats["higgs_mvapirho"], cats["higgs_mvaa1rho"], cats["higgs_mvaa1pi"], cats["higgs_mvaa1a1"], cats["higgs_mvapipi"], cats["higgs_mvapi0a1"], cats["higgs_mvarho0a1"], cats["higgs_mvaa10a1"])

    cats['higgs_mvapirho_up']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapirho_up"])
    cats['higgs_mvaa1pi_up']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvaa1pi_up"])
    cats['higgs_mvapipi_up']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapipi_up"])
    cats['higgs_mvapi0a1_up']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapi0a1_up"])

    cats['higgs_mvapirho_down']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapirho_down"])
    cats['higgs_mvaa1pi_down']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvaa1pi_down"])
    cats['higgs_mvapipi_down']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapipi_down"])
    cats['higgs_mvapi0a1_down']      = '({} && {})'.format(mva_ggh, cats["inclusive_mvapi0a1_down"])

########################## old below ##########################
#
#    cats["inclusive_rho"]       = "(tau_decay_mode_1==1 && tau_decay_mode_2==1)"
#    cats["inclusive_a1rho"]     = "((tau_decay_mode_1>=10 && tau_decay_mode_2==1) || (tau_decay_mode_1==1 && tau_decay_mode_2>=10))"
#    cats["mva_a1rho"]           = "((tau_decay_mode_1>=10 && mva_dm_1==10 && tau_decay_mode_2==1 && mva_dm_2==1) || (tau_decay_mode_1==1 && mva_dm_1==1 && tau_decay_mode_2>=10 && mva_dm_2==10))"
#    cats["mva_a1pi"]            = "((tau_decay_mode_1>=10 && mva_dm_1==10 && (tau_decay_mode_2==1 || tau_decay_mode_2==0) && mva_dm_2==0) || (tau_decay_mode_2>=10 && mva_dm_2==10 && (tau_decay_mode_1==1 || tau_decay_mode_1==0) && mva_dm_1==0) )"
#    cats["mva_0a1pi"]           = "(tau_decay_mode_1==1 && mva_dm_1==2 && (tau_decay_mode_2==1 || tau_decay_mode_2==0) && mva_dm_2==0 ) || (tau_decay_mode_2==1 && mva_dm_2==2 && (tau_decay_mode_1==1 || tau_decay_mode_1==0) && mva_dm_1==0 )"
#    cats["inclusive_a1"]        = "(tau_decay_mode_1>=10 && tau_decay_mode_2>=10)"
#    # lift conditions on DM (use DM=0||1 for 1pr, DM10||11 for 3pr)
#    # use in combination with MVA DM scores to determine genuine rhos/a1s
#    cats["DM0or1"]              = "((tau_decay_mode_1==1 || tau_decay_mode_1==0) && (tau_decay_mode_2==1 || tau_decay_mode_2==0))"
#    cats["mva_rhopi"]           = "( ((tau_decay_mode_1==1 || tau_decay_mode_1==0) && mva_dm_1==0 && tau_decay_mode_2==1 && mva_dm_2==1)  ||  ((tau_decay_mode_2==1 || tau_decay_mode_2==0) && mva_dm_2==0 && tau_decay_mode_1==1 && mva_dm_1==1) )"
#    cats["DM10or11and0or1"]     = "(((tau_decay_mode_1==10 || tau_decay_mode_1==11) && (tau_decay_mode_2==1 || tau_decay_mode_2==0)) || ((tau_decay_mode_1==1 || tau_decay_mode_1==0) && (tau_decay_mode_2==10 || tau_decay_mode_2==11)))"
#    # cats["dijet_rho"]             = "(n_jets>=2 && mjj>300)"
#    cats["idg0p5"]              = "(rho_id_1>0.5 && rho_id_2>0.5)"
#    cats["idgl0p5"]             = "!(rho_id_1>0.5 && rho_id_2>0.5)" #"((rho_id_1>0.5 && rho_id_2<0.5) || (rho_id_1<0.5 && rho_id_2>0.5))"
#    cats["idl0p5"]              = "(rho_id_1<0.5 && rho_id_2<0.5)"
#    cats["rho_idg0p5"]          = "{} && {}".format(cats["inclusive_rho"], cats["idg0p5"])
#    cats["rho_idl0p5"]          = "{} && !({})".format(cats["inclusive_rho"], cats["idg0p5"])
#
#    cats["mvadm_leadrho"]       = "(mvadm_rho_1>mvadm_a1_1 && mvadm_rho_1>mvadm_pi_1 && mvadm_rho_1>mvadm_other_1)"
#    cats["mvadm_leada1_1pr"]    = "(mvadm_a1_1>mvadm_rho_1 && mvadm_a1_1>mvadm_pi_1 && mvadm_a1_1>mvadm_other_1)"
#    cats["mvadm_subleadrho"]    = "(mvadm_rho_2>mvadm_a1_2 && mvadm_rho_2>mvadm_pi_2 && mvadm_rho_2>mvadm_other_2)"
#    cats["mvadm_subleada1_1pr"] = "(mvadm_a1_2>mvadm_rho_2 && mvadm_a1_2>mvadm_pi_2 && mvadm_a1_2>mvadm_other_2)"
#
#    cats["mvadm_leada1_3pr"]    = "(mvadm_a1_1>mvadm_3pipi0_1 && mvadm_a1_1>mvadm_other_1)"
#    cats["mvadm_subleada1_3pr"] = "(mvadm_a1_2>mvadm_3pipi0_2 && mvadm_a1_2>mvadm_other_2)"
#
#    cats["inclusive_mvarho"]    = "({} && {})"\
#            .format(cats["mvadm_leadrho"], cats["mvadm_subleadrho"])
#    cats["inclusive_mvaa1rho"]  = "(({} && {}) || ({} && {}))"\
#            .format(cats["mvadm_leada1_3pr"], cats["mvadm_subleadrho"], 
#                    cats["mvadm_subleada1_3pr"], cats["mvadm_leadrho"])
#    cats["inclusive_mvaa1"]     = "({} && {})"\
#            .format(cats["mvadm_leada1_1pr"], cats["mvadm_subleada1_1pr"])
#
#    # cut based categories
#    cats["0jet_rho"]              = "({} && {})".format(cats["0jet"], cats["inclusive_rho"])
#    cats["1jet_rho"]              = "(n_jets==1 && {})".format(cats["inclusive_rho"])
#    cats["2jet_rho"]              = "(n_jets==2 && {})".format(cats["inclusive_rho"])
#    cats["vbf_rho"]                   = "({} && mjj>500 && n_jets>=2)".format(cats["inclusive_rho"])
#    cats["dijet_rho"]             = "({} && {})".format(cats["dijet"] ,cats["inclusive_rho"])
#    cats["dijet_boosted_rho"]     = "({} && pt_tt>100 && {})".format(cats["dijet"] ,cats["inclusive_rho"])
#    cats["dijet_lowboost_rho"]    = "({} && pt_tt<100 && {})".format(cats["dijet"] ,cats["inclusive_rho"])
#    cats["boosted_rho"]           = "(!{} && !({}) && {})".format(cats["0jet"], cats["dijet_rho"], cats["inclusive_rho"])
#
#    cats["0jet_other"]              = "({} && !{})".format(cats["0jet"], cats["inclusive_rho"])
#    cats["dijet_boosted_other"]     = "({} && pt_tt>100 && !{})".format(cats["dijet"] ,cats["inclusive_rho"])
#    cats["dijet_lowboost_other"]    = "({} && pt_tt<100 && !{})".format(cats["dijet"] ,cats["inclusive_rho"])
#    cats["boosted_other"]           = "(!{} && !({}) && !{})".format(cats["0jet"], cats["dijet_rho"], cats["inclusive_rho"])
#
#    cats["0jet_rho_idg0p5"]           = "({} && {})".format(cats["0jet_rho"], cats["idg0p5"])
#    cats["1jet_rho_idg0p5"]           = "({} && {})".format(cats["1jet_rho"], cats["idg0p5"])
#    cats["2jet_rho_idg0p5"]           = "({} && {})".format(cats["2jet_rho"], cats["idg0p5"])
#    cats["dijet_boosted_rho_idg0p5"]  = "({} && {})".format(cats["dijet_boosted_rho"], cats["idg0p5"])
#    cats["dijet_lowboost_rho_idg0p5"] = "({} && {})".format(cats["dijet_lowboost_rho"], cats["idg0p5"])
#    cats["dijet_rho_idg0p5"]          = "({} && {})".format(cats["dijet_rho"], cats["idg0p5"])
#    cats["boosted_rho_idg0p5"]        = "({} && {})".format(cats["boosted_rho"], cats["idg0p5"])
#
#    cats["0jet_rho_idl0p5"]           = "({} && !({}))".format(cats["0jet_rho"], cats["idg0p5"])
#    cats["1jet_rho_idl0p5"]           = "({} && !({}))".format(cats["1jet_rho"], cats["idg0p5"])
#    cats["2jet_rho_idl0p5"]           = "({} && !({}))".format(cats["2jet_rho"], cats["idg0p5"])
#    cats["dijet_boosted_rho_idl0p5"]  = "({} && !({}))".format(cats["dijet_boosted_rho"], cats["idg0p5"])
#    cats["dijet_lowboost_rho_idl0p5"] = "({} && !({}))".format(cats["dijet_lowboost_rho"], cats["idg0p5"])
#    cats["dijet_rho_idl0p5"]          = "({} && !({}))".format(cats["dijet_rho"], cats["idg0p5"])
#    cats["boosted_rho_idl0p5"]        = "({} && !({}))".format(cats["boosted_rho"], cats["idg0p5"])

    # BDT multiclass categories 
    # old is Feb13_fix1
    # mva_ggh                = '(IC_Jun13_dR_tauspinner_split_max_index==0)'
    # mva_jetFakes           = '(IC_Jun13_dR_tauspinner_split_max_index==1)'
    # mva_qqh                = '(IC_Jun13_dR_tauspinner_split_max_index==2)'
    # mva_zttEmbed           = '(IC_Jun13_dR_tauspinner_split_max_index==3)'

    # mva_ggh                = '(IC_Feb13_fix1_max_index==0)'
    # mva_jetFakes           = '(IC_Feb13_fix1_max_index==1)'
    # mva_zttEmbed           = '(IC_Feb13_fix1_max_index==2)'

    # testing SM vs PS training on tauspinner samples
    # mva_ggh                = '(IC_Oct07_tauspinnerSM_max_index==0)'
    # mva_jetFakes           = '(IC_Oct07_tauspinnerSM_max_index==1)'
    # mva_zttEmbed           = '(IC_Oct07_tauspinnerSM_max_index==2)'

    # mva_ggh                = '(IC_Oct07_tauspinnerPS_max_index==0)'
    # mva_jetFakes           = '(IC_Oct07_tauspinnerPS_max_index==1)'
    # mva_zttEmbed           = '(IC_Oct07_tauspinnerPS_max_index==2)'

    # test tauspinner SM with split Higgs
    # mva_ggh                = '(IC_Oct07_tauspinnerSM_split_max_index==0)'
    # mva_jetFakes           = '(IC_Oct07_tauspinnerSM_split_max_index==1)'
    # mva_qqh                = '(IC_Oct07_tauspinnerSM_split_max_index==2)'
    # mva_zttEmbed           = '(IC_Oct07_tauspinnerSM_split_max_index==3)'

    # tauspinner with individual weights for signals, still single higgs
    # mva_ggh                = '(IC_Oct07_tauspinnerSM_individualSigWts_max_index==0)'
    # mva_jetFakes           = '(IC_Oct07_tauspinnerSM_individualSigWts_max_index==1)'
    # mva_zttEmbed           = '(IC_Oct07_tauspinnerSM_individualSigWts_max_index==2)'

   # # tauspinner with slimmed vars and classic bkg methods
   # # note the names are the same but actually jetFakes is only QCD
   # mva_ggh                = '(IC_Nov13_tauspinner_max_index==0)'
   # mva_jetFakes           = '(IC_Nov13_tauspinner_max_index==1)'
   # mva_zttEmbed           = '(IC_Nov13_tauspinner_max_index==2)'

   # cats['higgs']      = '({})'.format(mva_ggh)
   # cats['zttEmbed']   = '({})'.format(mva_zttEmbed)
   # cats['jetFakes']   = '({})'.format(mva_jetFakes)

   # cats['higgs_rho']      = '({} && {})'.format(mva_ggh, cats["inclusive_rho"])
   # # cats['ggh_rho']        = '({} && {})'.format(mva_ggh, cats["inclusive_rho"])
   # # cats['qqh_rho']        = '({} && {})'.format(mva_qqh, cats["inclusive_rho"])
   # cats['zttEmbed_rho']   = '({} && {})'.format(mva_zttEmbed, cats["inclusive_rho"])
   # cats['jetFakes_rho']   = '({} && {})'.format(mva_jetFakes, cats["inclusive_rho"])

   # cats['higgs_a1rho']    = '({} && {})'.format(mva_ggh, cats["inclusive_a1rho"])
   # # cats['ggh_a1rho']      = '({} && {})'.format(mva_ggh, cats["inclusive_a1rho"])
   # # cats['qqh_a1rho']      = '({} && {})'.format(mva_qqh, cats["inclusive_a1rho"])
   # cats['zttEmbed_a1rho'] = '({} && {})'.format(mva_zttEmbed, cats["inclusive_a1rho"])
   # cats['jetFakes_a1rho'] = '({} && {})'.format(mva_jetFakes, cats["inclusive_a1rho"])

   # cats['higgs_other']    = '({} && !({}||{}))'.format(mva_ggh, cats["inclusive_rho"], cats["inclusive_a1rho"])
   # # cats['ggh_other']      = '({} && !({}||{}))'.format(mva_ggh, cats["inclusive_rho"], cats["inclusive_a1rho"])
   # # cats['qqh_other']      = '({} && !({}||{}))'.format(mva_qqh, cats["inclusive_rho"], cats["inclusive_a1rho"])
   # cats['zttEmbed_other'] = '({} && !({}||{}))'.format(mva_zttEmbed, cats["inclusive_rho"], cats["inclusive_a1rho"])
   # cats['jetFakes_other'] = '({} && !({}||{}))'.format(mva_jetFakes, cats["inclusive_rho"], cats["inclusive_a1rho"])

   # # mvadm for rho and a1
   # cats['higgs_mvarho']      = '({} && {} && mva_dm_1==1 && mva_dm_2==1)'.format(mva_ggh, cats["inclusive_rho"])
   # cats['zttEmbed_mvarho']   = '({} && {} && mva_dm_1==1 && mva_dm_2==1 )'.format(mva_zttEmbed, cats["inclusive_rho"])
   # cats['jetFakes_mvarho']   = '({} && {} && mva_dm_1==1 && mva_dm_2==1)'.format(mva_jetFakes, cats["inclusive_rho"])

   # cats['higgs_mvaa1rho']      = '({} && {} )'.format(mva_ggh, cats["mva_a1rho"])
   # # cats['ggh_mvaa1rho']        = '({} && {} && {})'.format(mva_ggh, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
   # # cats['qqh_mvaa1rho']        = '({} && {} && {})'.format(mva_qqh, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
   # cats['zttEmbed_mvaa1rho']   = '({} && {} && {})'.format(mva_zttEmbed, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
   # cats['jetFakes_mvaa1rho']   = '({} && {} && {})'.format(mva_jetFakes, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
   # 
   # #rho-a1(1prong) and  a1(1prong)-a1(1prong)
   # cats['higgs_mvarho0A1'] = '{} && {} && ( (mva_dm_1==2 && (mva_dm_2==2 || mva_dm_2==1)) || (mva_dm_2==2 && (mva_dm_1==2 || mva_dm_1==1)) )'.format(mva_ggh, cats["inclusive_rho"])

   # #a1(3prong)-a1(3prong)
   # cats['higgs_mvaA1'] = '{} && {} && mva_dm_1==10 && mva_dm_2==10'.format(mva_ggh, cats["inclusive_a1"])

   # #pi-pi
   # cats["higgs_mvapi"]= '{} && {} && mva_dm_1==0 && mva_dm_2==0 && ip_sig_1>1.5 && ip_sig_2>1.5'.format(mva_ggh, cats["DM0or1"])

   # #rho-pi
   # cats["higgs_mvarhopi"]='{} && {}'.format(mva_ggh, cats["mva_rhopi"])

   # #pi-a1(3prong)
   # cats["higgs_mvaA1pi"]= '{} && {}'.format(mva_ggh, cats["mva_a1pi"])

   # #pi-a1(1prong)
   # cats["higgs_mva0A1pi"]='{} && {}'.format(mva_ggh, cats["mva_0a1pi"])

   # cats['higgs_mvaother']    = '({} && !({}||{}||{}||{}||{}||{}||{}||{}))'\
   #         .format(mva_ggh, cats["higgs_mvarho"], cats["higgs_mvaa1rho"], cats["higgs_mvarho0A1"], 
   #                 cats["higgs_mvaA1"], cats["higgs_mvapi"], cats["higgs_mvarhopi"],
   #                 cats["higgs_mvaA1pi"], cats["higgs_mva0A1pi"])
   # # cats['ggh_mvaother']    = '({} && !({}||{}||{}||{}))'\
   # #         .format(mva_ggh, cats["inclusive_rho"], cats["inclusive_a1rho"], 
   # #                 cats["inclusive_mvarho"], cats["inclusive_mvaa1rho"])
   # # cats['qqh_mvaother']    = '({} && !({}||{}||{}||{}))'\
   # #         .format(mva_qqh, cats["inclusive_rho"], cats["inclusive_a1rho"], 
   # #                 cats["inclusive_mvarho"], cats["inclusive_mvaa1rho"])
   # cats['zttEmbed_mvaother'] = '({} && !({}||{}||{}||{}))'\
   #         .format(mva_zttEmbed, cats["inclusive_rho"], cats["inclusive_a1rho"], 
   #                 cats["inclusive_mvarho"], cats["inclusive_mvaa1rho"])
   # cats['jetFakes_mvaother'] = '({} && !({}||{}||{}||{}))'\
   #         .format(mva_jetFakes, cats["inclusive_rho"], cats["inclusive_a1rho"], 
   #                 cats["inclusive_mvarho"], cats["inclusive_mvaa1rho"])

    # # mvadm for rho and a1 (lift reco DM conditions) -- NOT VERY GOOD ATM
    # # rho-rho
    # cats['higgs_mvarho']      = '({} && {} && {})'.format(mva_ggh, cats["DM0or1"], cats["inclusive_mvarho"])
    # cats['zttEmbed_mvarho']   = '({} && {} && {})'.format(mva_zttEmbed, cats["DM0or1"], cats["inclusive_mvarho"])
    # cats['jetFakes_mvarho']   = '({} && {} && {})'.format(mva_jetFakes, cats["DM0or1"], cats["inclusive_mvarho"])

    # # a1-a1 (1pr)
    # # cats['higgs_mvaa1']      = '({} && {} && {})'.format(mva_ggh, cats["DM0or1"], cats["inclusive_mvaa1"])
    # # cats['zttEmbed_mvaa1']   = '({} && {} && {})'.format(mva_zttEmbed, cats["DM0or1"], cats["inclusive_mvaa1"])
    # # cats['jetFakes_mvaa1']   = '({} && {} && {})'.format(mva_jetFakes, cats["DM0or1"], cats["inclusive_mvaa1"])

    # # a1(3pr)-rho
    # cats['higgs_mvaa1rho']      = '({} && {} && {})'.format(mva_ggh, cats["DM10or11and0or1"], cats["inclusive_mvaa1rho"])
    # cats['zttEmbed_mvaa1rho']   = '({} && {} && {})'.format(mva_zttEmbed, cats["DM10or11and0or1"], cats["inclusive_mvaa1rho"])
    # cats['jetFakes_mvaa1rho']   = '({} && {} && {})'.format(mva_jetFakes, cats["DM10or11and0or1"], cats["inclusive_mvaa1rho"])

    # cats['higgs_mvaother']    = '({} && !({}||{}||{}||{}))'\
    #         .format(mva_ggh, cats["DM0or1"], cats["DM10or11and0or1"], 
    #                 cats["inclusive_mvarho"], cats["inclusive_mvaa1rho"])
    # cats['zttEmbed_mvaother'] = '({} && !({}||{}||{}||{}))'\
    #         .format(mva_zttEmbed, cats["DM0or1"], cats["DM10or11and0or1"], 
    #                 cats["inclusive_mvarho"], cats["inclusive_mvaa1rho"])
    # cats['jetFakes_mvaother'] = '({} && !({}||{}||{}||{}))'\
    #         .format(mva_jetFakes, cats["DM0or1"], cats["DM10or11and0or1"], 
    #                 cats["inclusive_mvarho"], cats["inclusive_mvaa1rho"])

    # # both taus pass ID>0.5
    # cats['higgs_idg0p5']     = '({} && {})'.format(cats["higgs"], cats["idg0p5"])
    # cats['zttEmbed_idg0p5']  = '({} && {})'.format(cats["zttEmbed"], cats["idg0p5"])
    # cats['jetFakes_idg0p5']  = '({} && {})'.format(cats["jetFakes"], cats["idg0p5"])

    # # one tau passes ID>0.5
    # cats['higgs_idgl0p5']    = '({} && {})'.format(cats["higgs"], cats["idgl0p5"])
    # cats['zttEmbed_idgl0p5'] = '({} && {})'.format(cats["zttEmbed"], cats["idgl0p5"])
    # cats['jetFakes_idgl0p5'] = '({} && {})'.format(cats["jetFakes"], cats["idgl0p5"])

    # # both taus fail ID>0.5
    # cats['higgs_idl0p5']     = '({} && {})'.format(cats["higgs"], cats["idl0p5"])
    # cats['zttEmbed_idl0p5']  = '({} && {})'.format(cats["zttEmbed"], cats["idl0p5"])
    # cats['jetFakes_idl0p5']  = '({} && {})'.format(cats["jetFakes"], cats["idl0p5"])

    # # !(both taus pass ID>0.5)
    # cats['higgs_NOTidg0p5']     = '({} && !({}))'.format(cats["higgs"], cats["idg0p5"])
    # cats['zttEmbed_NOTidg0p5']  = '({} && !({}))'.format(cats["zttEmbed"], cats["idg0p5"])
    # cats['jetFakes_NOTidg0p5']  = '({} && !({}))'.format(cats["jetFakes"], cats["idg0p5"])

    
#    # MVA multiclass SM keras categories
#    nn_ggh              = '(IC_keras_sm7_max_index==0)'
#    nn_jetFakes         = '(IC_keras_sm7_max_index==1)'
#    nn_qqh              = '(IC_keras_sm7_max_index==2)'
#    nn_zttEmbed         = '(IC_keras_sm7_max_index==3)'
#
#    # nn_misc             = '(IC_keras_sm1_max_index==2)'
#    # nn_zttEmbed         = '(IC_keras_sm1_max_index==4)'
#
#    cats['NN_higgs']    = '({} && {})'.format(nn_ggh, cats["inclusive_rho"])
#    cats['NN_ggh']      = '({} && {})'.format(nn_ggh, cats["inclusive_rho"])
#    cats['NN_qqh']      = '({} && {})'.format(nn_qqh, cats["inclusive_rho"])
#    cats['NN_zttEmbed'] = '({} && {})'.format(nn_zttEmbed, cats["inclusive_rho"])
#    # cats['NN_misc']     = '({} && {})'.format(nn_misc, cats["inclusive_rho"])
#    cats['NN_jetFakes'] = '({} && {})'.format(nn_jetFakes, cats["inclusive_rho"])
#
#    cats['NN_higgs_other']    = '({} && !{})'.format(nn_ggh, cats["inclusive_rho"])
#    cats['NN_zttEmbed_other'] = '({} && !{})'.format(nn_zttEmbed, cats["inclusive_rho"])
#    cats['NN_jetFakes_other'] = '({} && !{})'.format(nn_jetFakes, cats["inclusive_rho"])
#
#    cats['NN_ggh_inclusive']      = '({})'.format(nn_ggh)
#    cats['NN_qqh_inclusive']      = '({})'.format(nn_qqh)
#    cats['NN_zttEmbed_inclusive'] = '({})'.format(nn_zttEmbed)
#    cats['NN_jetFakes_inclusive'] = '({})'.format(nn_jetFakes)
#
#    # MVA multiclass Vienna (25-01-2019 training)
#    nn_sm_ggh      = '(IC_Vienna_fix_max_index==0)'
#    nn_sm_qqh      = '(IC_Vienna_fix_max_index==1)'
#    nn_sm_zttEmbed = '(IC_Vienna_fix_max_index==2)'
#    nn_sm_jetFakes = '(IC_Vienna_fix_max_index==3)'
#    nn_sm_misc     = '(IC_Vienna_fix_max_index==4)'
#
#    cats['NN_sm_higgs_rho']    = '(({}||{}) && {})'.format(nn_sm_ggh, nn_sm_qqh, cats["inclusive_rho"])
#    cats['NN_sm_ggh_rho']      = '({} && {})'.format(nn_sm_ggh, cats["inclusive_rho"])
#    cats['NN_sm_qqh_rho']      = '({} && {})'.format(nn_sm_qqh, cats["inclusive_rho"])
#    cats['NN_sm_zttEmbed_rho'] = '({} && {})'.format(nn_sm_zttEmbed, cats["inclusive_rho"])
#    cats['NN_sm_jetFakes_rho'] = '({} && {})'.format(nn_sm_jetFakes, cats["inclusive_rho"])
#    cats['NN_sm_misc_rho']     = '({} && {})'.format(nn_sm_misc, cats["inclusive_rho"])
#
#    cats['NN_sm_higgs_a1rho']    = '(({}||{}) && {})'.format(nn_sm_ggh, nn_sm_qqh, cats["inclusive_a1rho"])
#    cats['NN_sm_ggh_a1rho']      = '({} && {})'.format(nn_sm_ggh, cats["inclusive_a1rho"])
#    cats['NN_sm_qqh_a1rho']      = '({} && {})'.format(nn_sm_qqh, cats["inclusive_a1rho"])
#    cats['NN_sm_zttEmbed_a1rho'] = '({} && {})'.format(nn_sm_zttEmbed, cats["inclusive_a1rho"])
#    cats['NN_sm_jetFakes_a1rho'] = '({} && {})'.format(nn_sm_jetFakes, cats["inclusive_a1rho"])
#    cats['NN_sm_misc_a1rho']     = '({} && {})'.format(nn_sm_misc, cats["inclusive_a1rho"])
#
#    cats['NN_sm_higgs_other']    = '(({}||{}) && !({}||{}))'.format(nn_sm_ggh, nn_sm_qqh, cats["inclusive_rho"], cats["inclusive_a1rho"])
#    cats['NN_sm_ggh_other']      = '({} && !({}||{}))'.format(nn_sm_ggh, cats["inclusive_rho"], cats["inclusive_a1rho"])
#    cats['NN_sm_qqh_other']      = '({} && !({}||{}))'.format(nn_sm_qqh, cats["inclusive_rho"], cats["inclusive_a1rho"])
#    cats['NN_sm_zttEmbed_other'] = '({} && !({}||{}))'.format(nn_sm_zttEmbed, cats["inclusive_rho"], cats["inclusive_a1rho"])
#    cats['NN_sm_jetFakes_other'] = '({} && !({}||{}))'.format(nn_sm_jetFakes, cats["inclusive_rho"], cats["inclusive_a1rho"])
#    cats['NN_sm_misc_other']     = '({} && !({}||{}))'.format(nn_sm_misc, cats["inclusive_rho"], cats["inclusive_a1rho"])
#
#
#    # with MVA DM for rho and a1
#    cats['NN_sm_higgs_mvarho']      = '(({}||{}) && {} && {})'.format(nn_sm_ggh, nn_sm_qqh, cats["inclusive_rho"], cats["inclusive_mvarho"])
#    cats['NN_sm_ggh_mvarho']        = '({} && {} && {})'.format(nn_sm_ggh, cats["inclusive_rho"], cats["inclusive_mvarho"])
#    cats['NN_sm_qqh_mvarho']        = '({} && {} && {})'.format(nn_sm_qqh, cats["inclusive_rho"], cats["inclusive_mvarho"])
#    cats['NN_sm_zttEmbed_mvarho']   = '({} && {} && {})'.format(nn_sm_zttEmbed,cats["inclusive_rho"], cats["inclusive_mvarho"])
#    cats['NN_sm_jetFakes_mvarho']   = '({} && {} && {})'.format(nn_sm_jetFakes,cats["inclusive_rho"], cats["inclusive_mvarho"])
#    cats['NN_sm_misc_mvarho']       = '({} && {} && {})'.format(nn_sm_misc, cats["inclusive_rho"], cats["inclusive_mvarho"])
#
#    cats['NN_sm_higgs_mvaa1rho']    = '(({}||{}) && {} && {})'.format(nn_sm_ggh, nn_sm_qqh, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_ggh_mvaa1rho']      = '({} && {} && {})'.format(nn_sm_ggh, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_qqh_mvaa1rho']      = '({} && {} && {})'.format(nn_sm_qqh, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_zttEmbed_mvaa1rho'] = '({} && {} && {})'.format(nn_sm_zttEmbed, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_jetFakes_mvaa1rho'] = '({} && {} && {})'.format(nn_sm_jetFakes, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_misc_mvaa1rho']     = '({} && {} && {})'.format(nn_sm_misc, cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#
#    cats['NN_sm_higgs_mvaother']    = '(({}||{}) && !({}||{}||{}))'.format(nn_sm_ggh, nn_sm_qqh, cats["inclusive_rho"], cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_ggh_mvaother']      = '({} && !({}||{}||{}))'.format(nn_sm_ggh, cats["inclusive_rho"], cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_qqh_mvaother']      = '({} && !({}||{}||{}))'.format(nn_sm_qqh, cats["inclusive_rho"], cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_zttEmbed_mvaother'] = '({} && !({}||{}||{}))'.format(nn_sm_zttEmbed, cats["inclusive_rho"], cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_jetFakes_mvaother'] = '({} && !({}||{}||{}))'.format(nn_sm_jetFakes, cats["inclusive_rho"], cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#    cats['NN_sm_misc_mvaother']     = '({} && !({}||{}||{}))'.format(nn_sm_misc, cats["inclusive_rho"], cats["inclusive_a1rho"], cats["inclusive_mvaa1rho"])
#
#    # cats['NN_sm_higgs_other']    = '(({}||{}) && !({}))'.format(nn_sm_ggh, nn_sm_qqh, cats["inclusive_rho"])
#    # cats['NN_sm_ggh_other']      = '({} && !({}))'.format(nn_sm_ggh, cats["inclusive_rho"])
#    # cats['NN_sm_qqh_other']      = '({} && !({}))'.format(nn_sm_qqh, cats["inclusive_rho"])
#    # cats['NN_sm_zttEmbed_other'] = '({} && !({}))'.format(nn_sm_zttEmbed, cats["inclusive_rho"])
#    # cats['NN_sm_jetFakes_other'] = '({} && !({}))'.format(nn_sm_jetFakes, cats["inclusive_rho"])
#    # cats['NN_sm_misc_other']     = '({} && !({}))'.format(nn_sm_misc, cats["inclusive_rho"])
#
#    cats['NN_sm_ggh_inclusive']      = '({})'.format(nn_sm_ggh)
#    cats['NN_sm_qqh_inclusive']      = '({})'.format(nn_sm_qqh)
#    cats['NN_sm_zttEmbed_inclusive'] = '({})'.format(nn_sm_zttEmbed)
#    cats['NN_sm_jetFakes_inclusive'] = '({})'.format(nn_sm_jetFakes)
#    cats['NN_sm_misc_inclusive']     = '({})'.format(nn_sm_misc)


# 2016 sm analysis uses relaxed shape selections for W + QCD processes in et and mt channel, these are set here
if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']: # Remove the False when finished!!!!!
    if options.channel in ['et','mt'] and options.cat in ['boosted','vbf','dijet','dijet_lowboost','dijet_boosted',
            # 'ggh_lowMjj','qqh_lowMjj','misc_lowMjj','qcd_lowMjj','qqh_lowMjj','tt_lowMjj','zll_lowMjj','ztt_lowMjj','fake_lowMjj','jetFakes_lowMjj','zttEmbed_lowMjj',
            # 'ggh_highMjj','ggh_boosted_highMjj','qqh_boosted_highMjj','ggh_tight_highMjj','ggh_loose_highMjj','ggh_loose_boost_highMjj','ggh_tight_boost_highMjj','qqh_highMjj','misc_highMjj','qcd_highMjj','qqh_highMjj','tt_highMjj','zll_highMjj','ztt_highMjj','fake_highMjj','jetFakes_highMjj','zttEmbed_highMjj','jetFakes_tight_highMjj','jetFakes_loose_highMjj'
            ]:
        cats['qcd_shape'] = '('+cats['baseline_loose']+')*('+cats[options.cat]+')'
    if options.cat in ['boosted','vbf','dijet','dijet_lowboost','dijet_boosted',
            # 'ggh_lowMjj','qqh_lowMjj','misc_lowMjj','qcd_lowMjj','qqh_lowMjj','tt_lowMjj','zll_lowMjj','ztt_lowMjj','fake_lowMjj','jetFakes_lowMjj','zttEmbed_lowMjj',
            # 'ggh_highMjj','ggh_boosted_highMjj','qqh_boosted_highMjj','ggh_tight_highMjj','ggh_loose_highMjj','ggh_loose_boost_highMjj','ggh_tight_boost_highMjj','qqh_highMjj','misc_highMjj','qcd_highMjj','qqh_highMjj','tt_highMjj','zll_highMjj','ztt_highMjj','fake_highMjj','jetFakes_highMjj','zttEmbed_highMjj','jetFakes_tight_highMjj','jetFakes_loose_highMjj'
            ]:
        cats['w_shape'] = cats['qcd_shape']


# Overwrite selection depending on whether tight or loose-mt categories is chosen - this can still be overwritten from command line using the --set_alias=sel:(...) option
#if options.cat == 'nobtag_tightmt' or options.cat == 'btag_tightmt':
#    if options.channel == 'mt' or options.channel == 'et': options.sel = '(mt_1<50)'
#if options.cat == 'nobtag_loosemt' or options.cat == 'btag_loosemt':
#    if options.channel == 'mt' or options.channel == 'et': options.sel = '(mt_1<70 && mt_1>50)'

if options.era == "mssmsummer16":
    if options.channel == "em": cats['baseline']+=" && trg_muonelectron"
    if options.channel == "et" or options.channel == 'zee': cats['baseline']+=" && trg_singleelectron"
    if options.channel in ['mt','zmm','mj']: cats['baseline']+=" && trg_singlemuon"
    if options.channel == "tt": cats['baseline']+=" && trg_doubletau"

# Overwrite any category selections if the --set_alias option is used
for i in options.set_alias:
    cat_to_overwrite = i.split(':')[0]
    cat_to_overwrite=cat_to_overwrite.replace("\"","")
    overwrite_with = i.split(':')[1]
    overwrite_with=overwrite_with.replace("\"","")
    start_index=overwrite_with.find("{")
    end_index=overwrite_with.find("}")
    while start_index >0:
        replace_with=overwrite_with[start_index:end_index+1]
        replace_with=replace_with.replace("{","")
        replace_with=replace_with.replace("}","")
        replace_string = cats[replace_with]
        overwrite_with=overwrite_with[0:start_index] + replace_string  + overwrite_with[end_index+1:]
        start_index=overwrite_with.find("{")
        end_index=overwrite_with.find("}")

    print 'Overwriting alias: \"'+cat_to_overwrite+'\" with selection: \"'+overwrite_with+'\"'
    if cat_to_overwrite == 'sel':
        options.sel = overwrite_with
    else:
        cats[cat_to_overwrite] = overwrite_with

# Additional selections to seperate MC samples by gen flags

z_sels = {}
if options.channel == 'et':
    z_sels['ztt_sel'] = '(gen_match_2==5)'
    z_sels['zl_sel'] = '(gen_match_2<5)'
    z_sels['zj_sel'] = '(gen_match_2==6)'
elif options.channel in ['mt','mj']:
    z_sels['ztt_sel'] = '(gen_match_2==5)'
    z_sels['zl_sel'] = '(gen_match_2<5)'
    z_sels['zj_sel'] = '(gen_match_2==6)'
elif options.channel in ['mj']:
    z_sels['ztt_sel'] = '(0)'
    z_sels['zl_sel'] = '(0)'
    z_sels['zj_sel'] = '(1)'
elif options.channel == 'tt':
    z_sels['ztt_sel'] = '(gen_match_1==5&&gen_match_2==5)'
    z_sels['zl_sel'] = '(gen_match_2<6&&gen_match_1<6&&!(gen_match_1==5&&gen_match_2==5))'
    z_sels['zj_sel'] = '(gen_match_2==6||gen_match_1==6)'
elif options.channel == 'em':
    z_sels['ztt_sel'] = '(gen_match_1>2 && gen_match_2>3)'
    z_sels['zll_sel'] = '(gen_match_1<3 || gen_match_2<4)'
elif options.channel == 'zee' or  options.channel == 'zmm':
    z_sels['ztt_sel'] = '(gen_match_1>2&&gen_match_1<6 && gen_match_2>2&&gen_match_2<6)'
    if options.channel == 'zmm': z_sels['zl_sel'] = '(gen_match_1==2&&gen_match_2==2)'
    else: z_sels['zl_sel'] = '(gen_match_1==1&&gen_match_2==1)'
    z_sels['zj_sel'] = '(!('+z_sels['zl_sel']+') && !('+z_sels['ztt_sel']+'))'

top_sels = {}
vv_sels = {}
top_sels['ttt_sel'] = z_sels['ztt_sel']
top_sels['ttj_sel'] = '!('+z_sels['ztt_sel']+')'
vv_sels['vvt_sel'] = z_sels['ztt_sel']
vv_sels['vvj_sel'] = '!('+z_sels['ztt_sel']+')'

if options.channel == 'zee' or  options.channel == 'zmm':
  top_sels['ttt_sel'] = z_sels['zl_sel']
  top_sels['ttj_sel'] = '!('+z_sels['zl_sel']+')'
  vv_sels['vvt_sel'] = z_sels['zl_sel']
  vv_sels['vvj_sel'] = '!('+z_sels['zl_sel']+')'


top_sels_embed = {}
if options.channel == 'mt': top_sels_embed['ttt_sel'] = '((gen_match_1 == 4) && (gen_match_2 == 5))'
if options.channel == 'et': top_sels_embed['ttt_sel'] = '((gen_match_1 == 3) && (gen_match_2 == 5))'
if options.channel == 'tt': top_sels_embed['ttt_sel'] = '((gen_match_1 == 5) && (gen_match_2 == 5))'
if options.channel == 'em': top_sels_embed['ttt_sel'] = '((gen_match_1 == 3) && (gen_match_2 == 4))'
if options.channel in ['zmm']: top_sels_embed['ttt_sel'] = '((gen_match_1 == 2) && (gen_match_2 == 2))'
if options.channel in ['zee']: top_sels_embed['ttt_sel'] = '((gen_match_1 == 1) && (gen_match_2 == 1))'

top_sels_embed['vvt_sel'] = top_sels_embed['ttt_sel']

if options.channel in ['et','mt','mj']:
  vv_sels['vvt_sel'] = '(gen_match_2<6)'
  vv_sels['vvj_sel'] = '(gen_match_2==6)'
  top_sels['ttt_sel'] = '(gen_match_2<6)' 
  top_sels['ttj_sel'] = '(gen_match_2==6)'
elif options.channel == 'tt':
  vv_sels['vvt_sel'] = '(gen_match_1<6 && gen_match_2<6)'
  vv_sels['vvj_sel'] = '(!(gen_match_1<6 && gen_match_2<6))'
  top_sels['ttt_sel'] = '(gen_match_1<6 && gen_match_2<6)' 
  top_sels['ttj_sel'] = '(!(gen_match_1<6 && gen_match_2<6))'
if options.channel in ['mj']:
  vv_sels['vvt_sel'] = '(0)'
  vv_sels['vvj_sel'] = '(1)'
  top_sels['ttt_sel'] = '(0)' 
  top_sels['ttj_sel'] = '(1)'
  
if options.era in ["smsummer16"]:
  if options.channel in ['mt','et']:
    z_sels['ztt_sel'] = '(gen_match_2==5)'
    z_sels['zl_sel'] = '(gen_match_2!=6&&gen_match_2!=5)'
    z_sels['zj_sel'] = '(gen_match_2==6)'  
    vv_sels['vvt_sel'] = '(gen_match_2==5)'
    vv_sels['vvj_sel'] = '(gen_match_2!=5)'
    top_sels['ttt_sel'] = '(gen_match_2==5)' 
    top_sels['ttj_sel'] = '(gen_match_2!=5)'    
  elif options.channel == 'tt':
    z_sels['ztt_sel'] = '(gen_match_1==5&&gen_match_2==5)'
    z_sels['zl_sel'] = '(!(gen_match_1==6 || gen_match_2==6) && !(gen_match_1==5&&gen_match_2==5))'
    z_sels['zj_sel'] = '(gen_match_1==6 || gen_match_2==6)'  
    vv_sels['vvt_sel'] = '(gen_match_1==5 && gen_match_2==5)'
    vv_sels['vvj_sel'] = '(!(gen_match_1==5 && gen_match_2==5))'
    top_sels['ttt_sel'] = '(gen_match_1==5 && gen_match_2==5)' 
    top_sels['ttj_sel'] = '(!(gen_match_1==5 && gen_match_2==5))'

if options.era in ['cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','tauid2016','cpsummer17','UL_17','tauid2017','cp18','UL_18','mvadm2016']:
  # define these selections to make them more compatible with the fake-factor method and embedded samples
  if options.channel =='em':
    z_sels['zll_sel'] = '(!(gen_match_2==4 && gen_match_1==3))'
    z_sels['ztt_sel'] = '(gen_match_2==4 && gen_match_1==3)'
    top_sels['ttt_sel'] = '(gen_match_2==4 && gen_match_1==3)'
    top_sels['ttj_sel'] = '!(gen_match_2==4 && gen_match_1==3)'
    vv_sels['vvt_sel'] = '(gen_match_2==4 && gen_match_1==3)'
    vv_sels['vvj_sel'] = '!(gen_match_2==4 && gen_match_1==3)'
  if options.channel =='mt':
    z_sels['ztt_sel'] = '(gen_match_2==5 && gen_match_1==4)'
  if options.channel =='et':
    z_sels['ztt_sel'] = '(gen_match_2==5 && gen_match_1==3)'
  if options.channel in ['mt','et']:
    z_sels['zj_sel'] = '(gen_match_2==6)'
    z_sels['zl_sel'] = '(gen_match_2!=6&&!(%s))' % z_sels['ztt_sel']
    vv_sels['vvt_sel'] = '(gen_match_2<6)'
    vv_sels['vvj_sel'] = '(gen_match_2==6)'
    top_sels['ttt_sel'] = '(gen_match_2<6)'
    top_sels['ttj_sel'] = '(gen_match_2==6)'
  elif options.channel == 'tt':
    z_sels['ztt_sel'] = '(gen_match_1==5&&gen_match_2==5)'
    z_sels['zl_sel'] = '(!(gen_match_1==6 || gen_match_2==6) && !(gen_match_1==5&&gen_match_2==5))'
    z_sels['zj_sel'] = '(gen_match_1==6 || gen_match_2==6)'
    vv_sels['vvt_sel'] = '(!(gen_match_1==6 || gen_match_2==6))'
    vv_sels['vvj_sel'] = '(gen_match_1==6 || gen_match_2==6)'
    top_sels['ttt_sel'] = '(!(gen_match_1==6 || gen_match_2==6))'
    top_sels['ttj_sel'] = '(gen_match_1==6 || gen_match_2==6)'
    if options.method==0 and True in ['baseline_aisotau1' in x for x in options.set_alias]:
      top_sels['ttj_sel'] = '(gen_match_1==6)'
      top_sels['ttt_sel'] = '(gen_match_1!=6)'
      vv_sels['vvj_sel'] = '(gen_match_1==6)'
      vv_sels['vvt_sel'] = '(gen_match_1!=6)'
      z_sels['zj_sel'] = '(gen_match_1==6)' 
      z_sels['zl_sel'] = '(!(gen_match_1==6) && !(gen_match_1==5&&gen_match_2==5))'
    if options.method==0 and True in ['baseline_aisotau2' in x for x in options.set_alias]:
      top_sels['ttj_sel'] = '(gen_match_2==6)'
      top_sels['ttt_sel'] = '(gen_match_2!=6)'
      vv_sels['vvj_sel'] = '(gen_match_2==6)'
      vv_sels['vvt_sel'] = '(gen_match_2!=6)'
      z_sels['zj_sel'] = '(gen_match_2==6)' 
      z_sels['zl_sel'] = '(!(gen_match_2==6) && !(gen_match_1==5&&gen_match_2==5))'
       

extra_top_sel = '1' 
if options.embedding:
    if options.channel == 'mt': extra_top_sel = '!((gen_match_1 == 4) && (gen_match_2 == 5))'
    if options.channel == 'et': extra_top_sel = '!((gen_match_1 == 3) && (gen_match_2 == 5))'
    if options.channel == 'tt': extra_top_sel = '!((gen_match_1 == 5) && (gen_match_2 == 5))'
    if options.channel == 'em': extra_top_sel = '!((gen_match_1 == 3) && (gen_match_2 == 4))'
    if options.channel =='zee': extra_top_sel = '!((gen_match_1 == 1) && (gen_match_2 == 1))'
    top_sels["ttt_ztt_sel"] = top_sels["ttt_sel"]+'&&!'+extra_top_sel
    vv_sels["vvt_ztt_sel"] = vv_sels["vvt_sel"]+'&&!'+extra_top_sel
    for sel in top_sels: 
      if "ztt_sel" not in sel: top_sels[sel]+='&&'+extra_top_sel
    for sel in vv_sels: 
      if "ztt_sel" not in sel: vv_sels[sel]+='&&'+extra_top_sel
    top_sels_embed["ttt_ztt_sel"] = top_sels["ttt_ztt_sel"]
    top_sels_embed["vvt_ztt_sel"] = vv_sels["vvt_ztt_sel"]
 
# Add data sample names
if options.channel == 'mt': 
    data_samples = ['SingleMuonB','SingleMuonC','SingleMuonD']
if options.channel == 'em': 
    data_samples = ['MuonEGB','MuonEGC','MuonEGD']
if options.channel == 'et': 
    data_samples = ['SingleElectronB','SingleElectronC','SingleElectronD']
if options.channel == 'tt': 
    data_samples = ['TauB','TauC','TauD']
    
# Add MC sample names   
ztt_samples = ['DYJetsToLL-LO-ext','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO']
vv_samples = ['T-tW', 'Tbar-tW','Tbar-t','T-t','WWTo1L1Nu2Q','WZJToLLLNu','VVTo2L2Nu','ZZTo2L2Q','ZZTo4L','WZTo2L2Q','WZTo1L3Nu','WZTo1L1Nu2Q']
wgam_samples = ['WGToLNuG','WGstarToLNuEE','WGstarToLNuMuMu']
top_samples = ['TT']
ztt_shape_samples = ['DYJetsToLL-LO-ext','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO']
wjets_samples = ['WJetsToLNu-LO','W1JetsToLNu-LO','W2JetsToLNu-LO','W3JetsToLNu-LO','W4JetsToLNu-LO']
ewkz_samples = []
gghww_samples = []
qqhww_samples = []

if options.era in ["mssmsummer16","smsummer16",'cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','tauid2016','mvadm2016']:
    # Add data sample names
    if options.channel in ['mt','zmm','mj']:
        data_samples = ['SingleMuonB','SingleMuonC','SingleMuonD','SingleMuonE','SingleMuonF','SingleMuonG','SingleMuonHv2','SingleMuonHv3']
	if options.era in ["legacy16"]: 
          data_samples = ['SingleMuonB','SingleMuonC','SingleMuonD','SingleMuonE','SingleMuonF','SingleMuonG','SingleMuonH']
          if options.analysis in ['mssmrun2','vlq'] and options.channel == 'mt': data_samples += ['TauB','TauC','TauD','TauE','TauF','TauG','TauH'] 
	if options.era in ['UL_16_preVFP']:
	  data_samples = ['SingleMuonB','SingleMuonC','SingleMuonD','SingleMuonE','SingleMuonF']
	  if options.analysis in ['mssmrun2','vlq'] and options.channel == 'mt': data_samples += ['TauB','TauC','TauD','TauE','TauF'] 
	if options.era in ['UL_16_postVFP']:
          data_samples = ['SingleMuonG','SingleMuonH']
	  if options.analysis in ['mssmrun2','vlq'] and options.channel == 'mt': data_samples += ['TauG','TauH']  

    if options.channel == 'em':
        data_samples = ['MuonEGB','MuonEGC','MuonEGD','MuonEGE','MuonEGF','MuonEGG','MuonEGHv2','MuonEGHv3']
        if options.era in ["legacy16"]: data_samples = ['MuonEGB','MuonEGC','MuonEGD','MuonEGE','MuonEGF','MuonEGG','MuonEGH']
	if options.era in ["UL_16_preVFP"]: data_samples = ['MuonEGB','MuonEGC','MuonEGD','MuonEGE','MuonEGF']
	if options.era in ["UL_16_postVFP"]: data_samples = ['MuonEGG','MuonEGH']
   
    if options.channel == 'et' or options.channel == 'zee':
        data_samples = ['SingleElectronB','SingleElectronC','SingleElectronD','SingleElectronE','SingleElectronF','SingleElectronG','SingleElectronHv2','SingleElectronHv3']
        if options.era in ["legacy16"]: 
          data_samples = ['SingleElectronB','SingleElectronC','SingleElectronD','SingleElectronE','SingleElectronF','SingleElectronG','SingleElectronH']
          if options.analysis in ['mssmrun2','vlq'] and options.channel == 'et': data_samples += ['TauB','TauC','TauD','TauE','TauF','TauG','TauH']
        if options.era in ["UL_16_preVFP"]: 
          data_samples = ['SingleElectronB','SingleElectronC','SingleElectronD','SingleElectronE','SingleElectronF']
          if options.analysis in ['mssmrun2','vlq'] and options.channel == 'et': data_samples += ['TauB','TauC','TauD','TauE','TauF']
        if options.era in ["UL_16_postVFP"]: 
          data_samples = ['SingleElectronG','SingleElectronH']
          if options.analysis in ['mssmrun2','vlq'] and options.channel == 'et': data_samples += ['TauG','TauH']


    if options.channel == 'tt':
        data_samples = ['TauB','TauC','TauD','TauE','TauF','TauG','TauHv2','TauHv3']
        if options.era in ["legacy16"]: data_samples = ['TauB','TauC','TauD','TauE','TauF','TauG','TauH']
	if options.era in ['UL_16_preVFP']: data_samples = ['TauB','TauC','TauD','TauE','TauF']
	if options.era in ['UL_16_postVFP']: data_samples = ['TauG','TauH']


    # Add MC sample names   
    if options.era in ["legacy16"]: 
	ztt_samples = ['DYJetsToLL-LO-ext1','DYJetsToLL-LO-ext2','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO']
        vv_samples = ['T-tW', 'Tbar-tW','Tbar-t','WWTo1L1Nu2Q','WZJToLLLNu','VVTo2L2Nu','VVTo2L2Nu-ext1','ZZTo2L2Q','ZZTo4L-amcat','WZTo2L2Q','WZTo1L3Nu','WZTo1L1Nu2Q']
	wgam_samples = ['WGToLNuG-ext1','WGToLNuG-ext2','WGToLNuG-ext3']

    if options.era in ["UL_16_preVFP"]:
	ztt_samples = ['DYJetsToLL','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10to50-LO','DYJetsToLL_M-10to50-LO','DY2JetsToLL_M-10to50-LO','DY3JetsToLL_M-10to50-LO','DY4JetsToLL_M-10to50-LO']
	vv_samples = ['WZTo1L1Nu2Q','WZTo3LNu','WZTo3LNu','WWTo2L2Nu','ZZTo2L2Nu','ZZTo4L','T-tW', 'Tbar-tW','Tbar-t','T-t']
	wgam_samples = ['WGToLNuG']

    if options.era in ["UL_16_postVFP"]:
        ztt_samples = ['DYJetsToLL','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10to50-LO','DYJetsToLL_M-10to50-LO','DY2JetsToLL_M-10to50-LO','DY3JetsToLL_M-10to50-LO','DY4JetsToLL_M-10to50-LO']
        vv_samples = ['WZTo3LNu','WZTo3LNu','WWTo2L2Nu','ZZTo2L2Nu','ZZTo4L','T-tW', 'Tbar-tW','Tbar-t','T-t']
        wgam_samples = ['WGToLNuG']

    if options.analysis in ['mssmrun2','vlq']: top_samples = ['TT']
    else: top_samples = ['TTTo2L2Nu', 'TTToHadronic', 'TTToSemiLeptonic']
    ztt_shape_samples = ['DYJetsToLL-LO-ext2','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO']
    wjets_samples = ['WJetsToLNu-LO', 'WJetsToLNu-LO-ext','W1JetsToLNu-LO','W2JetsToLNu-LO','W2JetsToLNu-LO-ext','W3JetsToLNu-LO','W3JetsToLNu-LO-ext','W4JetsToLNu-LO','W4JetsToLNu-LO-ext1','W4JetsToLNu-LO-ext2']
    if options.channel == 'mt': embed_samples = ['EmbeddingMuTauB','EmbeddingMuTauC','EmbeddingMuTauD','EmbeddingMuTauE','EmbeddingMuTauF','EmbeddingMuTauG','EmbeddingMuTauH']
    if options.channel == 'et': embed_samples = ['EmbeddingElTauB','EmbeddingElTauC','EmbeddingElTauD','EmbeddingElTauE','EmbeddingElTauF','EmbeddingElTauG','EmbeddingElTauH']
    if options.channel == 'em': embed_samples = ['EmbeddingElMuB','EmbeddingElMuC','EmbeddingElMuD','EmbeddingElMuE','EmbeddingElMuF','EmbeddingElMuG','EmbeddingElMuH']
    if options.channel == 'tt': embed_samples = ['EmbeddingTauTauB','EmbeddingTauTauC','EmbeddingTauTauD','EmbeddingTauTauE','EmbeddingTauTauF','EmbeddingTauTauG','EmbeddingTauTauH']
    if options.channel == 'zmm': embed_samples = ['EmbeddingMuMuB','EmbeddingMuMuC','EmbeddingMuMuD','EmbeddingMuMuE','EmbeddingMuMuF','EmbeddingMuMuG','EmbeddingMuMuH']
    if options.channel == 'zee': embed_samples = ['EmbeddingElElB','EmbeddingElElC','EmbeddingElElD','EmbeddingElElE','EmbeddingElElF','EmbeddingElElG','EmbeddingElElH']
    
if options.era in ["smsummer16",'cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','tauid2016','mvadm2016']:
    vv_samples = ['T-tW', 'Tbar-tW','Tbar-t','T-t','WWTo1L1Nu2Q','WZJToLLLNu','VVTo2L2Nu','VVTo2L2Nu-ext1','ZZTo2L2Q','ZZTo4L-amcat','WZTo2L2Q','WZTo1L3Nu','WZTo1L1Nu2Q']
    if options.era in ["legacy16",'UL_16_preVFP','UL_16_postVFP']:
        vv_samples = ['T-tW', 'Tbar-tW','Tbar-t','T-t','WWTo1L1Nu2Q','WZJToLLLNu','VVTo2L2Nu','VVTo2L2Nu-ext1','ZZTo2L2Q','ZZTo4L-amcat','WZTo2L2Q','WZTo1L3Nu','WZTo1L1Nu2Q']
        if options.analysis in ['mssmrun2','vlq']:
          vv_samples = ['T-tW', 'Tbar-tW','Tbar-t','T-t','WZJToLLLNu','VVTo2L2Nu','VVTo2L2Nu-ext1','ZZTo2L2Q','ZZTo4L-amcat','WZTo2L2Q']

    wjets_samples = ['WJetsToLNu-LO', 'WJetsToLNu-LO-ext','W1JetsToLNu-LO','W2JetsToLNu-LO','W2JetsToLNu-LO-ext','W3JetsToLNu-LO','W3JetsToLNu-LO-ext','W4JetsToLNu-LO','W4JetsToLNu-LO-ext1','W4JetsToLNu-LO-ext2', 'EWKWMinus2Jets_WToLNu','EWKWMinus2Jets_WToLNu-ext1','EWKWMinus2Jets_WToLNu-ext2','EWKWPlus2Jets_WToLNu','EWKWPlus2Jets_WToLNu-ext1','EWKWPlus2Jets_WToLNu-ext2']
    ewkz_samples = ['EWKZ2Jets_ZToLL','EWKZ2Jets_ZToLL-ext']
    if options.era in ["legacy16",'UL_16_preVFP','UL_16_postVFP']:
        ewkz_samples = ['EWKZ2Jets_ZToLL','EWKZ2Jets_ZToLL-ext1','EWKZ2Jets_ZToLL-ext2']

        wgam_samples = ['WGToLNuG-ext1','WGToLNuG-ext1','WGToLNuG-ext2','WGToLNuG-ext3','WGstarToLNuEE','WGstarToLNuMuMu']
    gghww_samples = ['GluGluHToWWTo2L2Nu_M-125']
    qqhww_samples = ['VBFHToWWTo2L2Nu_M-125']


# if options.era in ['cpdecay16']:
#     # Add data sample names
#     if options.channel in ['mt','zmm','mj']:
#         data_samples = ['SingleMuonB','SingleMuonC','SingleMuonD','SingleMuonE','SingleMuonF','SingleMuonG','SingleMuonH']
#     if options.channel == 'em':
#         data_samples = ['MuonEGB','MuonEGC','MuonEGD','MuonEGE','MuonEGF','MuonEGG','MuonEGH']
#     if options.channel == 'et' or options.channel == 'zee':
#         data_samples = ['SingleElectronB','SingleElectronC','SingleElectronD','SingleElectronE','SingleElectronF','SingleElectronG','SingleElectronH']
#     if options.channel == 'tt':
#         data_samples = ['TauB','TauC','TauD','TauE','TauF','TauG','TauH']

#     # Add MC sample names   
#     ztt_samples = ['DYJetsToLL-LO-ext1','DYJetsToLL-LO-ext2','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO']
#     #ztt_samples = ['DYJetsToLL'] # NL0 filelists 
#     wgam_samples = ['WGToLNuG','WGToLNuG-ext1','WGToLNuG-ext2','WGToLNuG-ext3','WGstarToLNuEE','WGstarToLNuMuMu']
#     top_samples = ['TT']
#     ztt_shape_samples = ['DYJetsToLL-LO-ext2','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO']
#     if options.channel == 'mt': embed_samples = ['EmbeddingMuTauB','EmbeddingMuTauC','EmbeddingMuTauD','EmbeddingMuTauE','EmbeddingMuTauF','EmbeddingMuTauG','EmbeddingMuTauH']
#     if options.channel == 'et': embed_samples = ['EmbeddingElTauB','EmbeddingElTauC','EmbeddingElTauD','EmbeddingElTauE','EmbeddingElTauF','EmbeddingElTauG','EmbeddingElTauH']
#     if options.channel == 'em': embed_samples = ['EmbeddingElMuB','EmbeddingElMuC','EmbeddingElMuD','EmbeddingElMuE','EmbeddingElMuF','EmbeddingElMuG','EmbeddingElMuH']
#     if options.channel == 'tt': embed_samples = ['EmbeddingTauTauB','EmbeddingTauTauC','EmbeddingTauTauD','EmbeddingTauTauE','EmbeddingTauTauF','EmbeddingTauTauG','EmbeddingTauTauH']
#     if options.channel == 'zmm': embed_samples = ['EmbeddingMuMuB','EmbeddingMuMuC','EmbeddingMuMuD','EmbeddingMuMuE','EmbeddingMuMuF','EmbeddingMuMuG','EmbeddingMuMuH']
#     if options.channel == 'zee': embed_samples = ['EmbeddingElElB','EmbeddingElElC','EmbeddingElElD','EmbeddingElElE','EmbeddingElElF','EmbeddingElElG','EmbeddingElElH']

#     vv_samples = ['T-tW', 'Tbar-tW','Tbar-t','T-t','WWTo1L1Nu2Q','WZJToLLLNu','VVTo2L2Nu','VVTo2L2Nu-ext1','ZZTo2L2Q','ZZTo4L-amcat','WZTo2L2Q','WZTo1L3Nu','WZTo1L1Nu2Q']
#     wjets_samples = ['WJetsToLNu-LO', 'WJetsToLNu-LO-ext','W1JetsToLNu-LO','W2JetsToLNu-LO','W2JetsToLNu-LO-ext','W3JetsToLNu-LO','W3JetsToLNu-LO-ext','W4JetsToLNu-LO','W4JetsToLNu-LO-ext1','W4JetsToLNu-LO-ext2', 'EWKWMinus2Jets_WToLNu','EWKWMinus2Jets_WToLNu-ext1','EWKWMinus2Jets_WToLNu-ext2','EWKWPlus2Jets_WToLNu','EWKWPlus2Jets_WToLNu-ext1','EWKWPlus2Jets_WToLNu-ext2']
#     ewkz_samples = ['EWKZ2Jets_ZToLL','EWKZ2Jets_ZToLL-ext1','EWKZ2Jets_ZToLL-ext2']
#     gghww_samples = ['GluGluHToWWTo2L2Nu_M-125']
#     qqhww_samples = ['VBFHToWWTo2L2Nu_M-125']

    
if options.era in ['cpsummer17','UL_17','tauid2017']:

    ztt_samples = ['DYJetsToLL-LO','DYJetsToLL-LO-ext1','DY1JetsToLL-LO','DY1JetsToLL-LO-ext','DY2JetsToLL-LO','DY2JetsToLL-LO-ext','DY3JetsToLL-LO','DY3JetsToLL-LO-ext','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO','DYJetsToLL_M-10-50-LO-ext1']
    if options.channel == "tt": # remove 'DYJetsToLL_M-10-50-LO (zero entries)
        ztt_samples = ['DYJetsToLL-LO','DYJetsToLL-LO-ext1','DY1JetsToLL-LO','DY1JetsToLL-LO-ext','DY2JetsToLL-LO','DY2JetsToLL-LO-ext','DY3JetsToLL-LO','DY3JetsToLL-LO-ext','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO-ext1']
    # ztt_samples = ['DYJetsToLL','DYJetsToLL-ext'] # NL0 filelists
    top_samples = ['TTTo2L2Nu', 'TTToHadronic', 'TTToSemiLeptonic']
    vv_samples = ['T-tW', 'Tbar-tW','Tbar-t','T-t','WWToLNuQQ','WZTo2L2Q','WZTo1L1Nu2Q','WZTo1L3Nu','WZTo3LNu', 'WWTo2L2Nu', 'ZZTo2L2Nu', 'ZZTo2L2Q','ZZTo4L-ext','ZZTo4L']
    if options.analysis in ['mssmrun2','vlq']:
      vv_samples = ['T-tW', 'Tbar-tW','Tbar-t','T-t','WZTo2L2Q','WZTo3LNu', 'ZZTo2L2Q','ZZTo4L-ext','ZZTo4L','VVTo2L2Nu']

    wjets_samples = ['WJetsToLNu-LO','WJetsToLNu-LO-ext','W1JetsToLNu-LO','W2JetsToLNu-LO','W3JetsToLNu-LO','W4JetsToLNu-LO','EWKWMinus2Jets','EWKWPlus2Jets']
    #wjets_samples = ['WJetsToLNu-LO','WJetsToLNu-LO-ext','EWKWMinus2Jets','EWKWPlus2Jets']
    wgam_samples = ['WGToLNuG']
    ewkz_samples = ['EWKZ2Jets']
    gghww_samples = []
    qqhww_samples = []
    gghww_samples = ['GluGluHToWWTo2L2Nu_M-125']
    qqhww_samples = ['VBFHToWWTo2L2Nu_M-125']   
 
    if options.channel in ['mt','zmm','mj']: 
        data_samples = ['SingleMuonB','SingleMuonC','SingleMuonD','SingleMuonE','SingleMuonF']
        if options.analysis in ['mssmrun2','vlq'] and options.channel == 'mt': data_samples += ['TauB','TauC','TauD','TauE','TauF']

    if options.channel == 'em': 
        data_samples = ['MuonEGB','MuonEGC','MuonEGD','MuonEGE','MuonEGF']
    if options.channel == 'et' or options.channel == 'zee': 
        data_samples = ['SingleElectronB','SingleElectronC','SingleElectronD','SingleElectronE','SingleElectronF']
        if options.analysis in ['mssmrun2','vlq'] and options.channel == 'et': data_samples += ['TauB','TauC','TauD','TauE','TauF']
    if options.channel == 'tt': 
        data_samples = ['TauB','TauC','TauD','TauE','TauF']

    if options.channel == 'mt': embed_samples = ['EmbeddingMuTauB','EmbeddingMuTauC','EmbeddingMuTauD','EmbeddingMuTauE','EmbeddingMuTauF']
    if options.channel == 'et': embed_samples = ['EmbeddingElTauB','EmbeddingElTauC','EmbeddingElTauD','EmbeddingElTauE','EmbeddingElTauF']
    if options.channel == 'em': embed_samples = ['EmbeddingElMuB','EmbeddingElMuC','EmbeddingElMuD','EmbeddingElMuE','EmbeddingElMuF']
    if options.channel == 'tt': embed_samples = ['EmbeddingTauTauB','EmbeddingTauTauC','EmbeddingTauTauD','EmbeddingTauTauE','EmbeddingTauTauF']
    if options.channel == 'zmm': embed_samples = ['EmbeddingMuMuB','EmbeddingMuMuC','EmbeddingMuMuD','EmbeddingMuMuE','EmbeddingMuMuF']
    if options.channel == 'zee': embed_samples = ['EmbeddingElElB','EmbeddingElElC','EmbeddingElElD','EmbeddingElElE','EmbeddingElElF']

if options.era in ['cp18','UL_18']:

    ztt_samples = ['DYJetsToLL-LO','DY1JetsToLL-LO','DY2JetsToLL-LO','DY3JetsToLL-LO','DY4JetsToLL-LO','DYJetsToLL_M-10-50-LO']
    # ztt_samples = ['DYJetsToLL-LO','DYJetsToLL_M-10-50-LO']
    # ztt_samples = ['DYJetsToLL'] # NL0 filelists
    top_samples = ['TTTo2L2Nu', 'TTToHadronic', 'TTToSemiLeptonic']
    vv_samples = [
            'T-tW-ext1', 'Tbar-tW-ext1','Tbar-t','WWTo2L2Nu','T-t',
            'WWToLNuQQ','WZTo1L3Nu','WZTo3LNu','WZTo3LNu-ext1','WZTo2L2Q',
            'ZZTo2L2Nu-ext1','ZZTo2L2Nu-ext2','ZZTo2L2Q','ZZTo4L-ext','ZZTo4L'
            ]
    if options.analysis in ['mssmrun2','vlq']:
      vv_samples = [
              'T-tW-ext1', 'Tbar-tW-ext1','Tbar-t','T-t',
              'WZTo3LNu','WZTo3LNu-ext1','WZTo2L2Q',
              'ZZTo2L2Q','ZZTo4L-ext','ZZTo4L', 'VVTo2L2Nu',
              ]

    wjets_samples = ['WJetsToLNu-LO','W1JetsToLNu-LO','W2JetsToLNu-LO','W3JetsToLNu-LO','W4JetsToLNu-LO','EWKWMinus2Jets','EWKWPlus2Jets']
    #wjets_samples = ['WJetsToLNu-LO','EWKWMinus2Jets','EWKWPlus2Jets']
    # wjets_samples = ['WJetsToLNu-LO',]
    wgam_samples = ['WGToLNuG']
    ewkz_samples = ['EWKZ2Jets']
    gghww_samples = []
    qqhww_samples = []
    gghww_samples = ['GluGluHToWWTo2L2Nu_M-125']
    qqhww_samples = ['VBFHToWWTo2L2Nu_M-125']   

 
    if options.channel in ['mt','zmm','mj']:
        data_samples = ['SingleMuonA','SingleMuonB','SingleMuonC','SingleMuonD']
        if options.analysis in ['mssmrun2','vlq']  and options.channel == 'mt': data_samples += ['TauA','TauB','TauC','TauD']
    if options.channel == 'em':
        data_samples = ['MuonEGA','MuonEGB','MuonEGC','MuonEGD']
    if options.channel == 'et' or options.channel == 'zee':
        data_samples = ['EGammaA','EGammaB','EGammaC','EGammaD']
        if options.analysis in ['mssmrun2','vlq'] and options.channel == 'et': data_samples += ['TauA','TauB','TauC','TauD']
    if options.channel == 'tt':
        data_samples = ['TauA','TauB','TauC','TauD']

    if options.channel == 'mt':  embed_samples = ['EmbeddingMuTauA','EmbeddingMuTauB','EmbeddingMuTauC','EmbeddingMuTauD']
    if options.channel == 'et':  embed_samples = ['EmbeddingElTauA','EmbeddingElTauB','EmbeddingElTauC','EmbeddingElTauD']
    if options.channel == 'em':  embed_samples = ['EmbeddingElMuA','EmbeddingElMuB','EmbeddingElMuC','EmbeddingElMuD']
    if options.channel == 'tt':  embed_samples = ['EmbeddingTauTauA','EmbeddingTauTauB','EmbeddingTauTauC','EmbeddingTauTauD']
    if options.channel == 'zmm': embed_samples = ['EmbeddingMuMuA','EmbeddingMuMuB','EmbeddingMuMuC','EmbeddingMuMuD']
    if options.channel == 'zee': embed_samples = ['EmbeddingElElA','EmbeddingElElB','EmbeddingElElC','EmbeddingElElD']

if options.method==0: ztt_samples+=ewkz_samples

sm_samples = { 'ggH' : 'GluGluHToTauTau_M-*', 'qqH' : 'VBFHToTauTau_M-*', 'WplusH' : 'WplusHToTauTau_M-*', 'WminusH' : 'WminusHToTauTau_M-*', 'ZH' : 'ZHToTauTau_M-*', 'TTH' : 'TTHToTauTau_M-*' }
if options.era in ["smsummer16"]: sm_samples = { 'ggH_htt*' : 'GluGluToHToTauTau_M-*', 'qqH_htt*' : 'VBFHToTauTau_M-*', 'WplusH_htt*' : 'WplusHToTauTau_M-*', 'WminusH_htt*' : 'WminusHToTauTau_M-*', 'ZH_htt*' : 'ZHToTauTau_M-*'}


if options.analysis in ['cpprod']: 
# 2016
  if options.era in ['legacy16','UL_16_preVFP','UL_16_postVFP']:
    sm_samples = { 

         'ggh*_powheg' : 'GluGluToHToTauTau_M-125', 
   #      'vbf*_powheg' : 'VBFHToTauTau_M-125',
    #     'wplush*_powheg': 'WplusHToTauTau_M-125',
     #    'wminush*_powheg': 'WminusHToTauTau_M-125',
      #   'zh*_powheg': 'ZHToTauTau_M-125',
#         'reweighted_ggH_htt_0PM*' : ['JJH0PMToTauTauPlusZeroJets_Filtered','JJH0PMToTauTauPlusOneJets_Filtered','JJH0PMToTauTauPlusTwoJets_Filtered'],
 #        'reweighted_ggH_htt_0Mf05ph0*' : ['JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered','JJH0Mf05ph0ToTauTauPlusOneJets_Filtered','JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered'],
  #       'reweighted_ggH_htt_0M*' : ['JJH0MToTauTauPlusZeroJets_Filtered','JJH0MToTauTauPlusOneJets_Filtered','JJH0MToTauTauPlusTwoJets_Filtered'],
         'ggH_sm_htt*' : ['JJH0PMToTauTauPlusZeroJets_Filtered','JJH0PMToTauTauPlusOneJets_Filtered','JJH0PMToTauTauPlusTwoJets_Filtered'],
         'ggH_mm_htt*' : ['JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered','JJH0Mf05ph0ToTauTauPlusOneJets_Filtered','JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered'],
         'ggH_ps_htt*' : ['JJH0MToTauTauPlusZeroJets_Filtered','JJH0MToTauTauPlusOneJets_Filtered','JJH0MToTauTauPlusTwoJets_Filtered'],
         'qqH_htt*' : 'VBFHToTauTau_M-125',
         'WH_htt*': ['WplusHToTauTau_M-125','WminusHToTauTau_M-125'],
         'ZH_htt*': 'ZHToTauTau_M-125',
       #  "ggH_cpdecay_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",

         'qqH_mm_htt*' : 'VBFHiggs0Mf05ph0_HToTauTau',
         'qqH_ps_htt*' : 'VBFHiggs0M_HToTauTau',
         'qqH_sm_htt*' : 'VBFHiggs0PM_HToTauTau',
         'WH_mm_htt*' : 'WHiggs0MfWH05ph0',
         'WH_ps_htt*' : 'WHiggs0M',
         'WH_sm_htt*' : 'WHiggs0PM',
         'ZH_mm_htt*' : 'ZHiggs0MfZH05ph0',
         'ZH_ps_htt*' : 'ZHiggs0M',
         'ZH_sm_htt*' : 'ZHiggs0PM',
    }

# 2017

  if (options.era == 'cpsummer17' or options.era == 'UL_17'):
    sm_samples = {
         'ggh*_powheg' : ['GluGluHToTauTau_M-125','GluGluHToTauTau_M-125-ext'],
     #    'vbf*_powheg' : 'VBFHToTauTau_M-125',
      #   'wplush*_powheg': 'WplusHToTauTau_M-125',
      #   'wminush*_powheg': 'WminusHToTauTau_M-125',
       #  'zh*_powheg': 'ZHToTauTau_M-125',
         #'ZHps_htt*' : 'ZHiggs0MToTauTau',
         #'ZHsm_htt*' : 'ZHiggs0PMToTauTau',
         #'ZHmm_htt*' : 'ZHiggs0Mf05ph0ToTauTau',
         #'qqHsm_htt*' : 'VBFHiggs0PMToTauTau',
         #'qqHps_htt*' :'VBFHiggs0MToTauTau',
         #'qqHmm_htt*' : 'VBFHiggs0Mf05ph0ToTauTau',
         #'WHps_htt*' :'WHiggs0MToTauTau',
         #'WHmm_htt*' : 'WHiggs0Mf05ph0ToTauTau',
         #'WHsm_htt*' : 'WHiggs0PMToTauTau',
        # 'reweighted_ggH_htt_0PM*' : ['JJH0PMToTauTauPlusZeroJets_Filtered','JJH0PMToTauTauPlusOneJets_Filtered','JJH0PMToTauTauPlusTwoJets_Filtered'],
        # 'reweighted_ggH_htt_0Mf05ph0*' : ['JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered','JJH0Mf05ph0ToTauTauPlusOneJets_Filtered','JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered'],
        # 'reweighted_ggH_htt_0M*' : ['JJH0MToTauTauPlusZeroJets_Filtered','JJH0MToTauTauPlusOneJets_Filtered','JJH0MToTauTauPlusTwoJets_Filtered'],
         'ggH_sm_htt*' : ['JJH0PMToTauTauPlusZeroJets_Filtered','JJH0PMToTauTauPlusOneJets_Filtered','JJH0PMToTauTauPlusTwoJets_Filtered'],
         'ggH_mm_htt*' : ['JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered','JJH0Mf05ph0ToTauTauPlusOneJets_Filtered','JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered'],
         'ggH_ps_htt*' : ['JJH0MToTauTauPlusZeroJets_Filtered','JJH0MToTauTauPlusOneJets_Filtered','JJH0MToTauTauPlusTwoJets_Filtered'],
         'qqH_htt*' : 'VBFHToTauTau_M-125',
         'WH_htt*': ['WplusHToTauTau_M-125','WminusHToTauTau_M-125'],
         'ZH_htt*': 'ZHToTauTau_M-125',
     #    "ggH_cpdecay_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
         'qqH_mm_htt*' : 'VBFHiggs0Mf05ph0ToTauTau',
         'qqH_ps_htt*' : 'VBFHiggs0MToTauTau',
         'qqH_sm_htt*' : 'VBFHiggs0PMToTauTau',
         'WH_mm_htt*' : 'WHiggs0Mf05ph0ToTauTau',
         'WH_ps_htt*' : 'WHiggs0MToTauTau',
         'WH_sm_htt*' : 'WHiggs0PMToTauTau',
         'ZH_mm_htt*' : 'ZHiggs0Mf05ph0ToTauTau',
         'ZH_ps_htt*' : 'ZHiggs0MToTauTau',
         'ZH_sm_htt*' : 'ZHiggs0PMToTauTau',
    }

# 2018

  if (options.era == 'cp18' or options.era == 'UL_18'):
    sm_samples = {
         #'vbf_new' : 'VBFHToTauTau_M125_withDipoleRecoil',
         'ggh*_powheg' : 'GluGluHToTauTau_M-125',
     #    'vbf*_powheg' : 'VBFHToTauTau_M-125-ext1',
      #   'wplush*_powheg': 'WplusHToTauTau_M-125',
       #  'wminush*_powheg': 'WminusHToTauTau_M-125',
        # 'zh*_powheg': 'ZHToTauTau_M-125',
         #'reweighted_ggH_htt_0PM*' : ['JJH0PMToTauTauPlusZeroJets_Filtered','JJH0PMToTauTauPlusOneJets_Filtered','JJH0PMToTauTauPlusTwoJets_Filtered'],
         #'reweighted_ggH_htt_0Mf05ph0*' : ['JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered','JJH0Mf05ph0ToTauTauPlusOneJets_Filtered','JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered'],
         #'reweighted_ggH_htt_0M*' : ['JJH0MToTauTauPlusZeroJets_Filtered','JJH0MToTauTauPlusOneJets_Filtered','JJH0MToTauTauPlusTwoJets_Filtered'],
         'ggH_sm_htt*' : ['JJH0PMToTauTauPlusZeroJets_Filtered','JJH0PMToTauTauPlusOneJets_Filtered','JJH0PMToTauTauPlusTwoJets_Filtered'],
         'ggH_mm_htt*' : ['JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered','JJH0Mf05ph0ToTauTauPlusOneJets_Filtered','JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered'],
         'ggH_ps_htt*' : ['JJH0MToTauTauPlusZeroJets_Filtered','JJH0MToTauTauPlusOneJets_Filtered','JJH0MToTauTauPlusTwoJets_Filtered'],
         'qqH_htt*' : 'VBFHToTauTau_M-125-ext1',
         'WH_htt*': ['WplusHToTauTau_M-125','WminusHToTauTau_M-125'],
         'ZH_htt*': 'ZHToTauTau_M-125',
        # "ggH_cpdecay_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
         'qqH_mm_htt*' : 'VBFHiggs0Mf05ph0ToTauTau',
         'qqH_ps_htt*' : 'VBFHiggs0MToTauTau',
         'qqH_sm_htt*' : 'VBFHiggs0PMToTauTau',
         'WH_mm_htt*' : 'WHiggs0Mf05ph0ToTauTau',
         'WH_ps_htt*' : 'WHiggs0MToTauTau',
         'WH_sm_htt*' : 'WHiggs0PMToTauTau',
         'ZH_mm_htt*' : 'ZHiggs0Mf05ph0ToTauTau',
         'ZH_ps_htt*' : 'ZHiggs0MToTauTau',
         'ZH_sm_htt*' : 'ZHiggs0PMToTauTau',
    }


if options.analysis in ['cpdecay']:

  if options.era in ['legacy16','UL_16_preVFP','UL_16_postVFP']:
    sm_samples = {
       #'ggH_ph_htt*' : 'GluGluToHToTauTau_M-125',
       #'qqH_ph_htt*' : 'VBFHToTauTau_M-125',
       "qqH_sm_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
       "qqH_ps_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
       "qqH_mm_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
       "ggH_sm_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
       "ggH_ps_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
       "ggH_mm_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
       "WH_sm_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
       "WH_ps_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
       "WH_mm_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
       "ZH_sm_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
       "ZH_ps_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
       "ZH_mm_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
       #"WH_ph_htt*": ['WplusHToTauTau_M-125','WminusHToTauTau_M-125'], 
       #"ZH_ph_htt*": 'ZHToTauTau_M-125', 
       #"ggH_flat_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
       #"qqH_flat_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
    }

# for 2017
  if (options.era == 'cpsummer17' or options.era == 'UL_17'):
    sm_samples = {
        'ggH_ph_htt*' : ['GluGluHToTauTau_M-125','GluGluHToTauTau_M-125-ext'],
        'qqH_ph_htt*' : 'VBFHToTauTau_M-125',
        "qqH_sm_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
        "qqH_ps_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
        "qqH_mm_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
        "ggH_sm_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
        "ggH_ps_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
        "ggH_mm_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered", 
        "WH_sm_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
        "WH_ps_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
        "WH_mm_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
        "ZH_sm_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
        "ZH_ps_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
        "ZH_mm_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
        "WH_ph_htt*": ['WplusHToTauTau_M-125','WminusHToTauTau_M-125'],
        "ZH_ph_htt*": 'ZHToTauTau_M-125',
        "ggH_flat_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
        "qqH_flat_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered", 
    }

  if (options.era == 'cp18' or options.era == 'UL_18'):
    sm_samples = {
        # test CP in decay samples
        'ggH_ph_htt*' : 'GluGluHToTauTau_M-125',
        'qqH_ph_htt*' : 'VBFHToTauTau_M-125-ext1',
        "qqH_sm_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
        "qqH_ps_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
        "qqH_mm_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
        "ggH_sm_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
        "ggH_ps_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
        "ggH_mm_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
        "WH_sm_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
        "WH_ps_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
        "WH_mm_htt*" : ['WplusHToTauTauUncorrelatedDecay_Filtered','WminusHToTauTauUncorrelatedDecay_Filtered'],
        "ZH_sm_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
        "ZH_ps_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
        "ZH_mm_htt*": 'ZHToTauTauUncorrelatedDecay_Filtered',
        "WH_ph_htt*": ['WplusHToTauTau_M-125','WminusHToTauTau_M-125'],
        "ZH_ph_htt*": 'ZHToTauTau_M-125',  
        "ggH_flat_htt*": "GluGluHToTauTauUncorrelatedDecay_Filtered",
        "qqH_flat_htt*": "VBFHToTauTauUncorrelatedDecay_Filtered",
    }

if options.analysis in ['mssmrun2','vlq']:
  if (options.era == 'cp18' or options.era == 'UL_18'):
    sm_samples = { 'ggH125_SM' : 'GluGluHToTauTau_M-125',
                   'qqH125' : ['VBFHToTauTau_M-125-ext1','ZHToTauTau_M-125','WplusHToTauTau_M-125','WminusHToTauTau_M-125'],
                   'qqH95' : 'VBFHToTauTau_M-95',
                   #'qqH125' : 'VBFHToTauTau_M-125-ext1',
                   #'ZH125' : 'ZHToTauTau_M-125',
                   #'WplusH125' : 'WplusHToTauTau_M-125',
                   #'WminusH125' : 'WminusHToTauTau_M-125',
                   #'ttH125' : 'ttHToTauTau',
                   #'ggHWW125' : 'GluGluHToWWTo2L2Nu_M-125',
                   #'qqHWW125' : 'VBFHToWWTo2L2Nu_M-125',
                   #'WminusHWW125' : 'HWminusJ_HToWW',
                   #'WplusHWW125' : 'HWplusJ_HToWW',
                 }
  elif (options.era == 'cpsummer17' or options.era == 'UL_17'):
    sm_samples = { 'ggH125_SM' : ['GluGluHToTauTau_M-125','GluGluHToTauTau_M-125-ext'],
                   'qqH125' : ['VBFHToTauTau_M-125','ZHToTauTau_M-125','WplusHToTauTau_M-125','WminusHToTauTau_M-125'],
                   'qqH95' : 'VBFHToTauTau_M-95',
                   #'qqH125' : 'VBFHToTauTau_M-125',
                   #'ZH125' : 'ZHToTauTau_M-125',
                   #'WplusH125' : 'WplusHToTauTau_M-125',
                   #'WminusH125' : 'WminusHToTauTau_M-125',
                   #'ttH125' : 'ttHToTauTau',
                   #'ggHWW125' : 'GluGluHToWWTo2L2Nu_M-125',
                   #'qqHWW125' : 'VBFHToWWTo2L2Nu_M-125',
                   #'WminusHWW125' : 'HWminusJ_HToWW',
                   #'WplusHWW125' : 'HWplusJ_HToWW',
                 }
  elif options.era in ['legacy16','UL_16_preVFP','UL_16_postVFP']:
    sm_samples = { 'ggH125_SM' : 'GluGluToHToTauTau_M-125',
                   'qqH125': ['VBFHToTauTau_M-125','ZHToTauTau_M-125','WplusHToTauTau_M-125','WminusHToTauTau_M-125'],
                   'qqH95' : 'VBFHToTauTau_M-95',
                   #'qqH125' : 'VBFHToTauTau_M-125',
                   #'ZH125' : 'ZHToTauTau_M-125',
                   #'WplusH125' : 'WplusHToTauTau_M-125',
                   #'WminusH125' : 'WminusHToTauTau_M-125',
                   #'ttH125' : 'ttHJetToTT',
                   #'ggHWW125' : 'GluGluHToWWTo2L2Nu_M-125',
                   #'qqHWW125' : 'VBFHToWWTo2L2Nu_M-125',
                   #'WminusHWW125' : 'HWminusJ_HToWW',
                   #'WplusHWW125' : 'HWplusJ_HToWW',
                 }

if options.analysis == "vlq":
  vlq_samples = {
    "vlq_betaRd33_0_mU2_gU1":"VectorLQToTauTau_betaRd33_0_mU2_gU1",
    "vlq_betaRd33_0_mU2_gU2":"VectorLQToTauTau_betaRd33_0_mU2_gU2",
    "vlq_betaRd33_0_mU2_gU3":"VectorLQToTauTau_betaRd33_0_mU2_gU3",
    "vlq_betaRd33_0_mU3_gU1":"VectorLQToTauTau_betaRd33_0_mU3_gU1",
    "vlq_betaRd33_0_mU3_gU2":"VectorLQToTauTau_betaRd33_0_mU3_gU2",
    "vlq_betaRd33_0_mU3_gU3":"VectorLQToTauTau_betaRd33_0_mU3_gU3",
    "vlq_betaRd33_0_mU4_gU1":"VectorLQToTauTau_betaRd33_0_mU4_gU1",
    "vlq_betaRd33_0_mU4_gU2":"VectorLQToTauTau_betaRd33_0_mU4_gU2",
    "vlq_betaRd33_0_mU4_gU3":"VectorLQToTauTau_betaRd33_0_mU4_gU3",
    "vlq_betaRd33_minus1_mU2_gU1":"VectorLQToTauTau_betaRd33_minus1_mU2_gU1",
    "vlq_betaRd33_minus1_mU2_gU2":"VectorLQToTauTau_betaRd33_minus1_mU2_gU2",
    "vlq_betaRd33_minus1_mU2_gU3":"VectorLQToTauTau_betaRd33_minus1_mU2_gU3",
    "vlq_betaRd33_minus1_mU3_gU1":"VectorLQToTauTau_betaRd33_minus1_mU3_gU1",
    "vlq_betaRd33_minus1_mU3_gU2":"VectorLQToTauTau_betaRd33_minus1_mU3_gU2",
    "vlq_betaRd33_minus1_mU3_gU3":"VectorLQToTauTau_betaRd33_minus1_mU3_gU3",
    "vlq_betaRd33_minus1_mU4_gU1":"VectorLQToTauTau_betaRd33_minus1_mU4_gU1",
    "vlq_betaRd33_minus1_mU4_gU2":"VectorLQToTauTau_betaRd33_minus1_mU4_gU2",
    "vlq_betaRd33_minus1_mU4_gU3":"VectorLQToTauTau_betaRd33_minus1_mU4_gU3",
    "VLQ_betaRd33_0_M_2000":["VectorLQToTauTau_betaRd33_0_mU2_gU1","VectorLQToTauTau_betaRd33_0_mU2_gU2","VectorLQToTauTau_betaRd33_0_mU2_gU3"],
    "VLQ_betaRd33_0_M_3000":["VectorLQToTauTau_betaRd33_0_mU3_gU1","VectorLQToTauTau_betaRd33_0_mU3_gU2","VectorLQToTauTau_betaRd33_0_mU3_gU3"],
    "VLQ_betaRd33_0_M_4000":["VectorLQToTauTau_betaRd33_0_mU4_gU1","VectorLQToTauTau_betaRd33_0_mU4_gU2","VectorLQToTauTau_betaRd33_0_mU4_gU3"],
    "VLQ_betaRd33_minus1_M_2000":["VectorLQToTauTau_betaRd33_minus1_mU2_gU1","VectorLQToTauTau_betaRd33_minus1_mU2_gU2","VectorLQToTauTau_betaRd33_minus1_mU2_gU3"],
    "VLQ_betaRd33_minus1_M_3000":["VectorLQToTauTau_betaRd33_minus1_mU3_gU1","VectorLQToTauTau_betaRd33_minus1_mU3_gU2","VectorLQToTauTau_betaRd33_minus1_mU3_gU3"],
    "VLQ_betaRd33_minus1_M_4000":["VectorLQToTauTau_betaRd33_minus1_mU4_gU1","VectorLQToTauTau_betaRd33_minus1_mU4_gU2","VectorLQToTauTau_betaRd33_minus1_mU4_gU3"],
  }


if options.ggh_masses_powheg == "":
  mssm_samples = { 'ggH' : 'SUSYGluGluToHToTauTau_M-*', 'bbH' : 'SUSYGluGluToBBHToTauTau_M-*' }
else:
  mssm_samples = { 'ggH' : 'SUSYGluGluToHToTauTau_M-*_powheg', 'bbH' : 'SUSYGluGluToBBHToTauTau_M-*_powheg' }

if options.bbh_masses_powheg == "":
  mssm_nlo_samples = { 'bbH' : 'SUSYGluGluToBBHToTauTau_M-*-NLO' }
  mssm_lo_samples = { 'bbH-LO' : 'SUSYGluGluToBBHToTauTau_M-*' }
  mssm_nlo_qsh_samples = { 'bbH-QshUp' : 'SUSYGluGluToBBHToTauTau_M-*-NLO-QshUp', 'bbH-QshDown' : 'SUSYGluGluToBBHToTauTau_M-*-NLO-QshDown' }
else:
  mssm_lo_samples = { 'bbH' : 'SUSYGluGluToBBHToTauTau_M-*_powheg' }


if options.nlo_qsh and mssm_nlo_samples: mssm_nlo_samples.update(mssm_nlo_qsh_samples)
Hhh_samples = { 'ggH' : 'GluGluToRadionToHHTo2B2Tau_M-*' }

# set systematics: first index sets folder name contaning systematic samples, second index sets string to be appended to output histograms, third index specifies the weight to be applied , 4th lists samples that should be skipped
systematics = OrderedDict()
if not options.no_default: systematics['default'] = ('','', 'wt', [], False)

if options.syst_e_res != '':
    systematics['res_e_up'] = ('ERES_UP' , '_'+options.syst_e_res+'Up', 'wt', ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)
    systematics['res_e_down'] = ('ERES_DOWN' , '_'+options.syst_e_res+'Down', 'wt', ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)
if options.syst_mu_res != '':
    systematics['res_mu_up'] = ('MURES_UP' , '_'+options.syst_mu_res+'Up', 'wt', ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)
    systematics['res_mu_down'] = ('MURES_DOWN' , '_'+options.syst_mu_res+'Down', 'wt', ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)
if options.syst_tau_res != '':
    systematics['res_tau_up'] = ('TRES_UP' , '_'+options.syst_tau_res+'Up', 'wt', ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)
    systematics['res_tau_down'] = ('TRES_DOWN' , '_'+options.syst_tau_res+'Down', 'wt', ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)
if options.syst_e_scale != '':
    systematics['scale_e_up'] = ('ESCALE_UP' , '_'+options.syst_e_scale+'Up', 'wt', ['QCD','jetFakes'], False)
    systematics['scale_e_down'] = ('ESCALE_DOWN' , '_'+options.syst_e_scale+'Down', 'wt', ['QCD','jetFakes'], False)
if options.syst_mu_scale != '':
    systematics['scale_mu_up'] = ('MUSCALE_UP' , '_'+options.syst_mu_scale+'Up', 'wt',['QCD','jetFakes'], False)
    systematics['scale_mu_down'] = ('MUSCALE_DOWN' , '_'+options.syst_mu_scale+'Down', 'wt', ['QCD','jetFakes'], False)

if options.syst_tau_scale != '':
    systematics['scale_t_up'] = ('TSCALE_UP' , '_'+options.syst_tau_scale+'Up', 'wt', ['jetFakes'], False)
    systematics['scale_t_down'] = ('TSCALE_DOWN' , '_'+options.syst_tau_scale+'Down', 'wt', ['jetFakes'], False)
if options.syst_tau_scale_0pi != '':
    systematics['scale_t_0pi_up'] = ('TSCALE0PI_UP' , '_'+options.syst_tau_scale_0pi+'Up', 'wt', ['QCD','jetFakes'], False)
    systematics['scale_t_0pi_down'] = ('TSCALE0PI_DOWN' , '_'+options.syst_tau_scale_0pi+'Down', 'wt', ['QCD','jetFakes'], False)
if options.syst_tau_scale_1pi != '':
    systematics['scale_t_1pi_up'] = ('TSCALE1PI_UP' , '_'+options.syst_tau_scale_1pi+'Up', 'wt', ['QCD','jetFakes'], False)
    systematics['scale_t_1pi_down'] = ('TSCALE1PI_DOWN' , '_'+options.syst_tau_scale_1pi+'Down', 'wt', ['QCD','jetFakes'], False)
if options.syst_tau_scale_3prong != '':
    systematics['scale_t_3prong_up'] = ('TSCALE3PRONG_UP' , '_'+options.syst_tau_scale_3prong+'Up', 'wt', ['QCD','jetFakes'], False)
    systematics['scale_t_3prong_down'] = ('TSCALE3PRONG_DOWN' , '_'+options.syst_tau_scale_3prong+'Down', 'wt', ['QCD','jetFakes'], False)
if options.syst_tau_scale_3prong1pi0 != '':
    systematics['scale_t_3prong1pi0_up'] = ('TSCALE3PRONG1PI0_UP' , '_'+options.syst_tau_scale_3prong1pi0+'Up', 'wt', ['QCD','jetFakes'], False)
    systematics['scale_t_3prong1pi0_down'] = ('TSCALE3PRONG1PI0_DOWN' , '_'+options.syst_tau_scale_3prong1pi0+'Down', 'wt', ['QCD','jetFakes'], False)
### group TES
if options.syst_tau_scale_grouped != "":
    hist_name = options.syst_tau_scale_grouped
    syst_dict = {
        '0pi': ['1prong','TSCALE0PI'],
        '1pi': ['1prong1pizero','TSCALE1PI'],
        '3prong': ['3prong','TSCALE3PRONG'],
        '3prong1pi0': ['3prong1pizero','TSCALE3PRONG1PI0'],
    }
    for name, values in syst_dict.iteritems():
        syst_name = hist_name.replace("*group", values[0])
        systematics["scale_t_{}_up".format(name)] = (
            "{}_UP".format(values[1]), "_{}Up".format(syst_name), 
            "wt", ["QCD","jetFakes"], False
        )
        systematics["scale_t_{}_down".format(name)] = (
            "{}_DOWN".format(values[1]), "_{}Down".format(syst_name), 
            "wt", ["QCD","jetFakes"], False
        )
##
if options.syst_efake_0pi_scale != '':
    systematics['scale_efake_0pi_up'] = ('EFAKE0PI_UP' , '_'+options.syst_efake_0pi_scale+'Up', 'wt', ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W','jetFakes','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
    systematics['scale_efake_0pi_down'] = ('EFAKE0PI_DOWN' , '_'+options.syst_efake_0pi_scale+'Down', 'wt', ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W','jetFakes','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
if options.syst_efake_1pi_scale != '':
    systematics['scale_efake_1pi_up'] = ('EFAKE1PI_UP' , '_'+options.syst_efake_1pi_scale+'Up', 'wt', ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W','jetFakes','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
    systematics['scale_efake_1pi_down'] = ('EFAKE1PI_DOWN' , '_'+options.syst_efake_1pi_scale+'Down', 'wt', ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W','jetFakes','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
if options.syst_mufake_0pi_scale != '':
    systematics['scale_mufake_0pi_up'] = ('MUFAKE0PI_UP' , '_'+options.syst_mufake_0pi_scale+'Up', 'wt', ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W','jetFakes','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
    systematics['scale_mufake_0pi_down'] = ('MUFAKE0PI_DOWN' , '_'+options.syst_mufake_0pi_scale+'Down', 'wt', ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W','jetFakes','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
if options.syst_mufake_1pi_scale != '':
    systematics['scale_mufake_1pi_up'] = ('MUFAKE1PI_UP' , '_'+options.syst_mufake_1pi_scale+'Up', 'wt', ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W','jetFakes','ggH_hww125','qqH_hww125','EmbedZTT'], False)
    systematics['scale_mufake_1pi_down'] = ('MUFAKE1PI_DOWN' , '_'+options.syst_mufake_1pi_scale+'Down', 'wt', ['ZTT','VVT','VVJ','TTT','TTJ','QCD','signal','W','jetFakes','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
if options.syst_eff_t != '':
    systematics['syst_eff_t_up'] = ('' , '_'+options.syst_eff_t+'Up', 'wt*wt_tau_id_up', ['ZL','ZJ','VVJ','TTJ','QCD','W'], False)
    systematics['syst_eff_t_down'] = ('' , '_'+options.syst_eff_t+'Down', 'wt*wt_tau_id_down', ['ZL','ZJ','VVJ','TTJ','QCD','W'], False)
if options.syst_tquark != '':
    systematics['syst_tquark_up'] = ('' , '_'+options.syst_tquark+'Up', 'wt*wt_tquark_up', ['ZTT','ZL','ZJ','VVT','VVJ','QCD','W','signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
    systematics['syst_tquark_down'] = ('' , '_'+options.syst_tquark+'Down', 'wt*wt_tquark_down', ['ZTT','ZL','ZJ','VVJ','VVT','QCD','W', 'signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)  
#    syst_alt_name='CMS_htt_ttbarShape_TheoryVsData'
#    if 'ttbarShape' in options.syst_tquark:
#      syst_alt_name = options.syst_tquark.replace('ttbarShape','ttbarShape_TheoryVsData') 
#    systematics['syst_tquark_alt_up'] = ('' , '_'+syst_alt_name+'Up', 'wt*wt_tquark_alt', ['ZTT','ZL','ZJ','VVT','VVJ','QCD','W','signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
#    systematics['syst_tquark_alt_down'] = ('' , '_'+syst_alt_name+'Down', 'wt*(2.-wt_tquark_alt)', ['ZTT','ZL','ZJ','VVJ','VVT','QCD','W', 'signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False) 
if options.syst_zwt != '':
    systematics['syst_zwt_up'] = ('' , '_'+options.syst_zwt+'Up', 'wt*wt_zpt_up', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
    systematics['syst_zwt_down'] = ('' , '_'+options.syst_zwt+'Down', 'wt*wt_zpt_down', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
if options.syst_w_fake_rate != '':
    to_skip = ['ZTT','ZL','ZJ','VVT','VVJ','TTT','TTJ','QCD','signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT']
    if options.era in ["smsummer16",'cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','tauid2016','tauid2017','mvadm2016']: to_skip = ['ZTT','ZL','VVT','TTT','QCD','signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT']
    systematics['syst_w_fake_rate_up'] = ('' , '_'+options.syst_w_fake_rate+'Up', 'wt*wt_tau_fake_up', to_skip, False)
    systematics['syst_w_fake_rate_down'] = ('' , '_'+options.syst_w_fake_rate+'Down', 'wt*wt_tau_fake_down', to_skip, False)
if options.syst_jfake_m != '':
    systematics['syst_jfake_m_up'] = ('' , '_'+options.syst_jfake_m+'Up', 'wt*idisoweight_up_2', ['ZTT','QCD','signal','TT','EmbedZTT'], False)
    systematics['syst_jfake_m_down'] = ('' , '_'+options.syst_jfake_m+'Down', 'wt*idisoweight_down_2', ['ZTT','QCD','signal','TT','EmbedZTT'], False)
if options.syst_jfake_e != '':
    systematics['syst_jfake_e_up'] = ('' , '_'+options.syst_jfake_e+'Up', 'wt*idisoweight_up_1', ['ZTT','QCD','signal','TT','EmbedZTT'], False)
    systematics['syst_jfake_e_down'] = ('' , '_'+options.syst_jfake_e+'Down', 'wt*idisoweight_down_1', ['ZTT','QCD','signal','TT','EmbedZTT'], False)
if options.syst_scale_j != '':
    systematics['syst_scale_j_up'] = ('JES_UP' , '_'+options.syst_scale_j+'Up', 'wt', ['EmbedZTT','jetFakes'], False)
    systematics['syst_scale_j_down'] = ('JES_DOWN' , '_'+options.syst_scale_j+'Down', 'wt', ['EmbedZTT','jetFakes'], False)
if options.syst_res_j != '':
    systematics['syst_res_j_up'] = ('JER_UP' , '_'+options.syst_res_j+'Up', 'wt', ['EmbedZTT','jetFakes'], False)
    systematics['syst_res_j_down'] = ('JER_DOWN' , '_'+options.syst_res_j+'Down', 'wt', ['EmbedZTT','jetFakes'], False)
if options.syst_scale_j_corr != '':
    systematics['syst_scale_j_corr_up'] = ('JES_CORR_UP' , '_'+options.syst_scale_j_corr+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_corr_down'] = ('JES_CORR_DOWN' , '_'+options.syst_scale_j_corr+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_uncorr != '':
    systematics['syst_scale_j_uncorr_up'] = ('JES_UNCORR_UP' , '_'+options.syst_scale_j_uncorr+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_uncorr_down'] = ('JES_UNCORR_DOWN' , '_'+options.syst_scale_j_uncorr+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_rbal != '':
    systematics['syst_scale_j_rbal_up'] = ('JESRBAL_UP' , '_'+options.syst_scale_j_rbal+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_rbal_down'] = ('JESRBAL_DOWN' , '_'+options.syst_scale_j_rbal+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_rsamp != '':
    systematics['syst_scale_j_rsamp_up'] = ('JESRSAMP_UP' , '_'+options.syst_scale_j_rsamp+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_rsamp_down'] = ('JESRSAMP_DOWN' , '_'+options.syst_scale_j_rsamp+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_full != '':
    systematics['syst_scale_j_full_up'] = ('JESFULL_UP' , '_'+options.syst_scale_j_full+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_full_down'] = ('JESFULL_DOWN' , '_'+options.syst_scale_j_full+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_cent != '':
    systematics['syst_scale_j_cent_up'] = ('JESCENT_UP' , '_'+options.syst_scale_j_cent+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_cent_down'] = ('JESCENT_DOWN' , '_'+options.syst_scale_j_cent+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_hf != '':
    systematics['syst_scale_j_hf_up'] = ('JESHF_UP' , '_'+options.syst_scale_j_hf+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_hf_down'] = ('JESHF_DOWN' , '_'+options.syst_scale_j_hf+'Down', 'wt', ['EmbedZTT'], False)    
if options.syst_scale_j_full_corr != '':
    systematics['syst_scale_j_full_corr_up'] = ('JESFULL_CORR_UP' , '_'+options.syst_scale_j_full_corr+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_full_corr_down'] = ('JESFULL_CORR_DOWN' , '_'+options.syst_scale_j_full_corr+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_cent_corr != '':
    systematics['syst_scale_j_cent_corr_up'] = ('JESCENT_CORR_UP' , '_'+options.syst_scale_j_cent_corr+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_cent_corr_down'] = ('JESCENT_CORR_DOWN' , '_'+options.syst_scale_j_cent_corr+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_hf_corr != '':
    systematics['syst_scale_j_hf_corr_up'] = ('JESHF_CORR_UP' , '_'+options.syst_scale_j_hf_corr+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_hf_corr_down'] = ('JESHF_CORR_DOWN' , '_'+options.syst_scale_j_hf_corr+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_full_uncorr != '':
    systematics['syst_scale_j_full_uncorr_up'] = ('JESFULL_UNCORR_UP' , '_'+options.syst_scale_j_full_uncorr+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_full_uncorr_down'] = ('JESFULL_UNCORR_DOWN' , '_'+options.syst_scale_j_full_uncorr+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_cent_uncorr != '':
    systematics['syst_scale_j_cent_uncorr_up'] = ('JESCENT_UNCORR_UP' , '_'+options.syst_scale_j_cent_uncorr+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_cent_uncorr_down'] = ('JESCENT_UNCORR_DOWN' , '_'+options.syst_scale_j_cent_uncorr+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_scale_j_hf_uncorr != '':
    systematics['syst_scale_j_hf_uncorr_up'] = ('JESHF_UNCORR_UP' , '_'+options.syst_scale_j_hf_uncorr+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_scale_j_hf_uncorr_down'] = ('JESHF_UNCORR_DOWN' , '_'+options.syst_scale_j_hf_uncorr+'Down', 'wt', ['EmbedZTT'], False)
############ Regrouped JEC for full run2
if options.syst_scale_j_regrouped != "":
    hist_name = options.syst_scale_j_regrouped
    # need dict of syst names and folders of where the shifted trees are found
    names = ["Absolute", "Absolute_year", "BBEC1", "BBEC1_year", 
            "EC2", "EC2_year", "FlavorQCD", "HF", "HF_year", 
            "RelativeBal", "RelativeSample_year"]
    folders = ["JESABS", "JESABS_YEAR", "JESBBEC1", "JESBBEC1_YEAR", 
            "JESEC2", "JESEC2_YEAR", "JESFLAV", "JESHF", "JESHF_YEAR", 
            "JESRBAL", "JESRELSAMP_YEAR"]
    syst_dict = dict(zip(names, folders))

    replaceYear = ""
    if options.era in ["legacy16",'UL_16_preVFP','UL_16_postVFP']: replaceYear = "2016"
    elif (options.era == "cpsummer17" or options.era == 'UL_17'): replaceYear = "2017"
    elif (options.era == "cp18" or options.era == 'UL_18'): replaceYear = "2018"
    else: assert ValueError("Regrouped JES only works for full RunII analyses")

    if "*group" in hist_name:
        for name, folder in syst_dict.iteritems():
            if "year" in name: name=name.replace("year", replaceYear)
            syst_name = hist_name.replace("*group", name)
            systematics['syst_scale_j_{}_up'.format(syst_name)] = (
                "{}_UP".format(folder), "_{}Up".format(syst_name), 
                "wt", ["EmbedZTT",'jetFakes','QCD'], False)
            systematics['syst_scale_j_{}_down'.format(syst_name)] = (
                "{}_DOWN".format(folder), "_{}Down".format(syst_name), 
                "wt", ["EmbedZTT",'jetFakes','QCD'], False)

if options.syst_eff_b_weights != '':
    systematics['syst_b_weights_up'] = ('' , '_'+options.syst_eff_b_weights+'Up', 'wt*wt_btag_up/wt_btag', ['EmbedZTT','ZTT','ZL','ZLL','ZJ','EWKZ','signal','jetFakes','W','QCD','qqH_hww','ggH_hww'], False)
    systematics['syst_b_weights_down'] = ('' , '_'+options.syst_eff_b_weights+'Down', 'wt*wt_btag_down/wt_btag', ['EmbedZTT','ZTT','ZL','ZLL','ZJ','EWKZ','signal','jetFakes','W','QCD','qqH_hww','ggH_hww'], False)
if options.syst_eff_b != '':
    systematics['syst_b_up'] = ('../../../../../../../vols/cms/dw515/Offline/output/MSSM/mssm_{}_btag/BTAG_UP'.format(options.year), '_'+options.syst_eff_b+'Up', 'wt', ['EmbedZTT','ZJ','EWKZ','jetFakes','W','QCD','qqH_hww','ggH_hww'], False)
    systematics['syst_b_down'] = ('../../../../../../../vols/cms/dw515/Offline/output/MSSM/mssm_{}_btag/BTAG_DOWN'.format(options.year) , '_'+options.syst_eff_b+'Down', 'wt', ['EmbedZTT','ZJ','EWKZ','jetFakes','W','QCD','qqH_hww','ggH_hww'], False)
if options.syst_fake_b != '':
    systematics['syst_fake_b_up'] = ('../../../../../../../vols/cms/dw515/Offline/output/MSSM/mssm_{}_btag/BFAKE_UP'.format(options.year), '_'+options.syst_fake_b+'Up', 'wt', ['EmbedZTT'], False)
    systematics['syst_fake_b_down'] = ('../../../../../../../vols/cms/dw515/Offline/output/MSSM/mssm_{}_btag/BFAKE_DOWN'.format(options.year) , '_'+options.syst_fake_b+'Down', 'wt', ['EmbedZTT'], False)
if options.syst_zpt_es != '':
    systematics['syst_zpt_es_up'] = ('' , '_'+options.syst_zpt_es+'Up', 'wt*wt_zpt_esup', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
    systematics['syst_zpt_es_down'] = ('' , '_'+options.syst_zpt_es+'Down', 'wt*wt_zpt_esdown', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
if options.syst_zpt_tt != '':
    systematics['syst_zpt_tt_up'] = ('' , '_'+options.syst_zpt_tt+'Up', 'wt*wt_zpt_ttup', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
    systematics['syst_zpt_tt_down'] = ('' , '_'+options.syst_zpt_tt+'Down', 'wt*wt_zpt_ttdown', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
if options.syst_zpt_statpt0 != '':
    systematics['syst_zpt_statpt0_up'] = ('' , '_'+options.syst_zpt_statpt0+'Up', 'wt*wt_zpt_stat_m400pt0_up', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
    systematics['syst_zpt_statpt0_down'] = ('' , '_'+options.syst_zpt_statpt0+'Down', 'wt*wt_zpt_stat_m400pt0_down', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
if options.syst_zpt_statpt40 != '':
    systematics['syst_zpt_statpt40_up'] = ('' , '_'+options.syst_zpt_statpt40+'Up', 'wt*wt_zpt_stat_m400pt40_up', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
    systematics['syst_zpt_statpt40_down'] = ('' , '_'+options.syst_zpt_statpt40+'Down', 'wt*wt_zpt_stat_m400pt40_down', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
if options.syst_zpt_statpt80 != '':
    systematics['syst_zpt_statpt80_up'] = ('' , '_'+options.syst_zpt_statpt80+'Up', 'wt*wt_zpt_stat_m400pt80_up', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
    systematics['syst_zpt_statpt80_down'] = ('' , '_'+options.syst_zpt_statpt80+'Down', 'wt*wt_zpt_stat_m400pt80_down', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
if options.syst_z_mjj != '' and options.cat in ['vbf','dijet','dijet_lowboost','dijet_boosted',
        # 'ggh_lowMjj','qqh_lowMjj','misc_lowMjj','qcd_lowMjj','qqh_lowMjj','tt_lowMjj','zll_lowMjj','ztt_lowMjj','fake_lowMjj','jetFakes_lowMjj','zttEmbed_lowMjj',
        # 'ggh_highMjj','ggh_boosted_highMjj','qqh_boosted_highMjj','ggh_tight_highMjj','ggh_loose_highMjj','ggh_loose_boost_highMjj','ggh_tight_boost_highMjj','qqh_highMjj','misc_highMjj','qcd_highMjj','qqh_highMjj','tt_highMjj','zll_highMjj','ztt_highMjj','fake_highMjj','jetFakes_highMjj','zttEmbed_highMjj','jetFakes_tight_highMjj','jetFakes_loose_highMjj'
        ]:
    systematics['syst_z_mjj_up'] = ('' , '_'+options.syst_z_mjj+'Up', 'wt*wt_z_mjj_up', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','ggH_hww125','qqH_hww125', 'ggH_hww', 'qqH_hww','EmbedZTT'], False)
    systematics['syst_z_mjj_down'] = ('' , '_'+options.syst_z_mjj+'Down', 'wt*wt_z_mjj_down', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','ggH_hww125','qqH_hww125', 'ggH_hww', 'qqH_hww','EmbedZTT'], False)
if options.syst_qcd_scale != '' and options.cat in ['0jet','boosted','vbf','dijet','dijet_lowboost','dijet_boosted',
            ] and options.channel in ['em','et','mt','tt'] and options.era == 'smsummer16': 
    if "Mjj" not in options.cat:
        weight_up = 'wt*wt_scale_%s_%s' % (options.channel, options.cat)
        weight_down = 'wt*(2-wt_scale_%s_%s)' % (options.channel, options.cat)
    elif "low" in options.cat:
        weight_up = 'wt*wt_scale_%s_boosted' % (options.channel)
        weight_down = 'wt*(2-wt_scale_%s_boosted)' % (options.channel)
    if options.cat in ['dijet','dijet_lowboost','dijet_boosted',
            ]:
      weight_up = 'wt*wt_scale_%s_vbf' % (options.channel)
      weight_down = 'wt*(2-wt_scale_%s_vbf)' % (options.channel)

    systematics['syst_qcd_scale_up'] = ('' , '_'+options.syst_qcd_scale+'Up', weight_up, ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','qqH','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
    systematics['syst_qcd_scale_down'] = ('' , '_'+options.syst_qcd_scale+'Down', weight_down, ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','qqH','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
if options.syst_qcd_scale != '' and options.era != 'smsummer16':
    systematics['syst_qcd_scale_up'] = ('' , '_'+options.syst_qcd_scale+'Up', 'wt*wt_qcdscale_up', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*'], False)
    systematics['syst_qcd_scale_down'] = ('' , '_'+options.syst_qcd_scale+'Down', 'wt*wt_qcdscale_down', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*'], False)
if options.syst_quarkmass != '':
    systematics['syst_quarkmass_up'] = ('' , '_'+options.syst_quarkmass+'Up', 'wt*wt_quarkmass', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','qqH','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*'], False)
    systematics['syst_quarkmass_down'] = ('' , '_'+options.syst_quarkmass+'Down', 'wt*((1./wt_quarkmass-1.)*0.05 + 1.)', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','qqH','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*'], False)
if options.syst_ps != '':
    if '*PS' in options.syst_ps:
      hist_name_fsr = options.syst_ps.replace('*PS','PS_FSR')
      hist_name_isr = options.syst_ps.replace('*PS','PS_ISR')
      systematics['syst_ps_fsr_up'] = ('' , '_'+hist_name_fsr+'Up', 'wt*min(wt_ps_fsr_up,10)', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*','ggHsm_jhu_htt*','ggHps_jhu_htt*','ggHmm_jhu_htt*'], False)
      systematics['syst_ps_fsr_down'] = ('' , '_'+hist_name_fsr+'Down', 'wt*min(wt_ps_fsr_down,10)', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*','ggHsm_jhu_htt*','ggHps_jhu_htt*','ggHmm_jhu_htt*'], False)
      systematics['syst_ps_isr_up'] = ('' , '_'+hist_name_isr+'Up', 'wt*min(wt_ps_isr_up,10)', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*','ggHsm_jhu_htt*','ggHps_jhu_htt*','ggHmm_jhu_htt*'], False)
      systematics['syst_ps_isr_down'] = ('' , '_'+hist_name_isr+'Down', 'wt*min(wt_ps_isr_down,10)', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*','ggHsm_jhu_htt*','ggHps_jhu_htt*','ggHmm_jhu_htt*'], False)
    else: 
      systematics['syst_ps_up'] = ('' , '_'+options.syst_ps+'Up', 'wt*wt_ps_up', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','qqH','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*','ggHsm_jhu_htt*','ggHps_jhu_htt*','ggHmm_jhu_htt*'], False)
      systematics['syst_ps_down'] = ('' , '_'+options.syst_ps+'Down', 'wt*wt_ps_down', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','qqH','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*','ggHsm_jhu_htt*','ggHps_jhu_htt*','ggHmm_jhu_htt*'], False)
if options.syst_ue != '':
    systematics['syst_ue_up'] = ('' , '_'+options.syst_ue+'Up', 'wt*wt_ue_up', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','qqH','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*','ggHsm_jhu_htt*','ggHps_jhu_htt*','ggHmm_jhu_htt*'], False)
    systematics['syst_ue_down'] = ('' , '_'+options.syst_ue+'Down', 'wt*wt_ue_down', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','QCD','W','jetFakes','qqH','WminusH','WplusH','ZH','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT','ggH_ph_htt*','ggHsm_jhu_htt*','ggHps_jhu_htt*','ggHmm_jhu_htt*'], False)
if options.syst_tau_id_dm0 != '':
    systematics['syst_tau_id_dm0_up'] = ('' , '_'+options.syst_tau_id_dm0+'Up', 'wt*wt_tau_id_dm0_up', ['ZL','ZJ','VVJ','TTJ','QCD','W'], False)
    systematics['syst_tau_id_dm0_down'] = ('' , '_'+options.syst_tau_id_dm0+'Down', 'wt*wt_tau_id_dm0_down', ['ZL','ZJ','VVJ','TTJ','QCD','W'], False)
if options.syst_tau_id_dm1 != '':
    systematics['syst_tau_id_dm1_up'] = ('' , '_'+options.syst_tau_id_dm1+'Up', 'wt*wt_tau_id_dm1_up', ['ZL','ZJ','VVJ','TTJ','QCD','W'], False)
    systematics['syst_tau_id_dm1_down'] = ('' , '_'+options.syst_tau_id_dm1+'Down', 'wt*wt_tau_id_dm1_down', ['ZL','ZJ','VVJ','TTJ','QCD','W'], False)    
if options.syst_tau_id_dm10 != '':
    systematics['syst_tau_id_dm10_up'] = ('' , '_'+options.syst_tau_id_dm10+'Up', 'wt*wt_tau_id_dm10_up', ['ZL','ZJ','VVJ','TTJ','QCD','W'], False)
    systematics['syst_tau_id_dm10_down'] = ('' , '_'+options.syst_tau_id_dm10+'Down', 'wt*wt_tau_id_dm10_down', ['ZL','ZJ','VVJ','TTJ','QCD','W'], False)   
if options.syst_lfake_dm0 != '':
    systematics['syst_lfake_dm0_up'] = ('' , '_'+options.syst_lfake_dm0+'Up', 'wt*wt_lfake_dm0_up', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
    systematics['syst_lfake_dm0_down'] = ('' , '_'+options.syst_lfake_dm0+'Down', 'wt*wt_lfake_dm0_down', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)    
if options.syst_lfake_dm1 != '':
    systematics['syst_lfake_dm1_up'] = ('' , '_'+options.syst_lfake_dm1+'Up', 'wt*wt_lfake_dm1_up', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
    systematics['syst_lfake_dm1_down'] = ('' , '_'+options.syst_lfake_dm1+'Down', 'wt*wt_lfake_dm1_down', ['VVT','VVJ','TTT','TTJ','QCD','W','signal','jetFakes','EmbedZTT'], False)
if options.syst_qcd_shape_wsf != '':
    systematics['syst_qcd_shape_wsf_up'] = ('' , '_'+options.syst_qcd_shape_wsf.replace('cat',options.cat)+'Up', 'wt', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','jetFakes','signal','W','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
    systematics['syst_qcd_shape_wsf_down'] = ('' , '_'+options.syst_qcd_shape_wsf.replace('cat',options.cat)+'Down', 'wt', ['ZTT','ZL','ZJ','ZLL','VVT','VVJ','TTT','TTJ','jetFakes','signal','W','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
    if options.cat in ["0jet","boosted",
            # 'ggh_lowMjj','qqh_lowMjj','misc_lowMjj','qcd_lowMjj','qqh_lowMjj','tt_lowMjj','zll_lowMjj','ztt_lowMjj','fake_lowMjj','jetFakes_lowMjj','zttEmbed_lowMjj',
            ]:
        w_abs_shift=0.1
    if options.cat in ["vbf",'dijet','dijet_lowboost','dijet_boosted',
            # 'ggh_highMjj','ggh_boosted_highMjj','qqh_boosted_highMjj','ggh_tight_highMjj','ggh_loose_highMjj','ggh_loose_boost_highMjj','ggh_tight_boost_highMjj','qqh_highMjj','misc_highMjj','qcd_highMjj','qqh_highMjj','tt_highMjj','zll_highMjj','ztt_highMjj','fake_highMjj','jetFakes_highMjj','zttEmbed_highMjj','jetFakes_tight_highMjj','jetFakes_loose_highMjj'
            ]:
        w_abs_shift=0.3
if options.syst_scale_met_unclustered != '':
    systematics['syst_scale_met_unclustered_up'] = ('METUNCL_UP' , '_'+options.syst_scale_met_unclustered+'Up', 'wt', ['EWKZ','ZLL','ZL','ZJ','ZTT','W','signal','QCD','jetFakes','EmbedZTT','ggH_hww','qqH_hww'], False)
    systematics['syst_scale_met_unclustered_down'] = ('METUNCL_DOWN' , '_'+options.syst_scale_met_unclustered+'Down', 'wt', ['EWKZ','ZLL','ZL','ZJ','ZTT','W','signal','QCD','jetFakes','EmbedZTT','ggH_hww','qqH_hww'], False)
if options.syst_scale_met_clustered != '':
    systematics['syst_scale_met_clustered_up'] = ('METCL_UP' , '_'+options.syst_scale_met_clustered+'Up', 'wt', ['QCD','jetFakes','EmbedZTT'], False)
    systematics['syst_scale_met_clustered_down'] = ('METCL_DOWN' , '_'+options.syst_scale_met_clustered+'Down', 'wt', ['QCD','jetFakes','EmbedZTT'], False)
if options.syst_res_met != '':
    hist_name = options.syst_res_met
    if '$NJET' in hist_name: 
      hist_name_0jet = hist_name.replace('$NJET','0Jet')
      hist_name_1jet = hist_name.replace('$NJET','1Jet')
      hist_name_2jet = hist_name.replace('$NJET','2Jet')
      
      systematics['syst_res_met_0jet_up'] = ('MET_RES_NJETS0_UP' , '_'+hist_name_0jet+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      systematics['syst_res_met_0jet_down'] = ('MET_RES_NJETS0_DOWN' , '_'+hist_name_0jet+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      
      systematics['syst_res_met_1jet_up'] = ('MET_RES_NJETS1_UP' , '_'+hist_name_1jet+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      systematics['syst_res_met_1jet_down'] = ('MET_RES_NJETS1_DOWN' , '_'+hist_name_1jet+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      
      systematics['syst_res_met_2jet_up'] = ('MET_RES_NJETS2_UP' , '_'+hist_name_2jet+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      systematics['syst_res_met_2jet_down'] = ('MET_RES_NJETS2_DOWN' , '_'+hist_name_2jet+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
    elif options.analysis == "vlq":
      systematics['syst_res_met_up'] = ('MET_RES_UP' , '_'+hist_name+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','signal'], False)
      systematics['syst_res_met_down'] = ('MET_RES_DOWN' , '_'+hist_name+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','signal'], False)
    else:
      systematics['syst_res_met_up'] = ('MET_RES_UP' , '_'+hist_name+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww'], False)
      systematics['syst_res_met_down'] = ('MET_RES_DOWN' , '_'+hist_name+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww'], False) 
if options.syst_scale_met != '':
    hist_name = options.syst_scale_met
    if '$NJET' in hist_name: 
      hist_name_0jet = hist_name.replace('$NJET','0Jet')
      hist_name_1jet = hist_name.replace('$NJET','1Jet')
      hist_name_2jet = hist_name.replace('$NJET','2Jet')
      systematics['syst_scale_met_0jet_up'] = ('MET_SCALE_NJETS0_UP' , '_'+hist_name_0jet+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      systematics['syst_scale_met_0jet_down'] = ('MET_SCALE_NJETS0_DOWN' , '_'+hist_name_0jet+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      
      systematics['syst_scale_met_1jet_up'] = ('MET_SCALE_NJETS1_UP' , '_'+hist_name_1jet+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      systematics['syst_scale_met_1jet_down'] = ('MET_SCALE_NJETS1_DOWN' , '_'+hist_name_1jet+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      
      systematics['syst_scale_met_2jet_up'] = ('MET_SCALE_NJETS2_UP' , '_'+hist_name_2jet+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
      systematics['syst_scale_met_2jet_down'] = ('MET_SCALE_NJETS2_DOWN' , '_'+hist_name_2jet+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ'], False)
    elif options.analysis == "vlq":
      systematics['syst_scale_met_up'] = ('MET_SCALE_UP' , '_'+hist_name+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','signal'], False)
      systematics['syst_scale_met_down'] = ('MET_SCALE_DOWN' , '_'+hist_name+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','signal'], False)
    else:  
      systematics['syst_scale_met_up'] = ('MET_SCALE_UP' , '_'+hist_name+'Up', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww'], False)
      systematics['syst_scale_met_down'] = ('MET_SCALE_DOWN' , '_'+hist_name+'Down', 'wt', ['QCD','jetFakes','EmbedZTT','TT','TTJ','TTT','VV','VVT','VVJ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww'], False)  
    if options.embedding and options.analysis in ['mssmrun2','vlq']: 
      embed_hist_name='scale_embed_met_%s' % options.year 
      systematics['syst_scale_embed_met_up'] = ('MET_SCALE_UP' , '_'+embed_hist_name+'Up', 'wt', ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)
      systematics['syst_scale_embed_met_down'] = ('MET_SCALE_DOWN' , '_'+embed_hist_name+'Down', 'wt', ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False) 
if options.syst_scale_j_by_source != '':
    jes_sources={"AbsoluteFlavMap":1,"AbsoluteMPFBias":2,"AbsoluteScale":3,"AbsoluteStat":4,"FlavorQCD":5,"Fragmentation":6,"PileUpDataMC":7,"PileUpPtBB":8,"PileUpPtEC1":9,"PileUpPtEC2":10,"PileUpPtHF":11,"PileUpPtRef":12,"RelativeBal":13,"RelativeFSR":14,"RelativeJEREC1":15,"RelativeJEREC2":16,"RelativeJERHF":17,"RelativePtBB":18,"RelativePtEC1":19,"RelativePtEC2":20,"RelativePtHF":21,"RelativeStatEC":22,"RelativeStatFSR":23,"RelativeStatHF":24,"SinglePionECAL":25,"SinglePionHCAL":26,"TimePtEta":27}
    jes_to_process=[]
    for i in options.jes_sources.split(','):
      if ':' in i: jes_to_process+=range(int(i.split(':')[0]),int(i.split(':')[1])+1)
      else: jes_to_process.append(int(i)) 
    jes_to_process = list(set(jes_to_process))  
    for source in jes_sources:
      jes_num = jes_sources[source]  
      if jes_num not in jes_to_process: continue  
      replace_dict = {'n_jets':'n_jets_%i'%jes_num, 'n_bjets':'n_bjets_%i'%jes_num, 'mjj':'mjj_%i'%jes_num, 'jdeta':'jdeta_%i'%jes_num, 'jdphi':'jdphi_%i'%jes_num, 'jpt_1':'jpt_1_%i'%jes_num, 'jpt_2':'jpt_2_%i'%jes_num}
      syst_name = 'syst_scale_j_by_source_'+source
      hist_name = options.syst_scale_j_by_source.replace('SOURCE', source)
      systematics[syst_name+'_up'] = ('JES_UP' , '_'+hist_name+'Up', 'wt', ['jetFakes','EmbedZTT'], False,replace_dict)
      systematics[syst_name+'_down'] = ('JES_DOWN' , '_'+hist_name+'Down', 'wt', ['jetFakes','EmbedZTT'], False,replace_dict)

## em QCD uncertainties
if options.syst_em_qcd != '' and options.channel == 'em':
    hist_name = options.syst_em_qcd
    if '*BIN' in hist_name:
      hist_name_bini = hist_name.replace('*BIN', 'IsoExtrap')
      systematics['syst_em_qcd_extrap_up'] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_em_qcd_extrapup', ['ZLL','TT','TTJ','TTT','ZTT','ZL','ZJ','VVT','VVJ','W','signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
      systematics['syst_em_qcd_extrap_down'] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_em_qcd_extrapdown', ['ZLL','TT','TTJ','TTT','ZTT','ZL','ZJ','VVT','VVJ','W','signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
      for j in range(0,3):
        for i in range(1,3):
          hist_name_bini = hist_name.replace('*BIN', 'stat_njets%(j)i_unc%(i)i' % vars())
          systematics['syst_em_qcd_njets%(j)i_unc%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_em_qcd_njets%(j)s_unc%(i)i_up' % vars(), ['ZLL','TT','TTJ','TTT','ZTT','ZL','ZJ','VVT','VVJ','W','signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)
          systematics['syst_em_qcd_njets%(j)i_unc%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_em_qcd_njets%(j)i_unc%(i)i_down' % vars(), ['ZLL','TT','TTJ','TTT','ZTT','ZL','ZJ','VVT','VVJ','W','signal','jetFakes','EWKZ','ggH_hww125','qqH_hww125','ggH_hww','qqH_hww','EmbedZTT'], False)  

if options.syst_prefire != '':
    systematics['syst_prefire_up'] = ('' , '_'+options.syst_prefire+'Up', 'wt*wt_prefire_up/wt_prefire', ['QCD','jetFakes','EmbedZTT'], False)
    systematics['syst_prefire_down'] = ('' , '_'+options.syst_prefire+'Down', 'wt*wt_prefire_down/wt_prefire', ['QCD','jetFakes','EmbedZTT'], False)

if options.syst_tau_id_diff != '':

    hist_name = options.syst_tau_id_diff

    if options.analysis in ['mssmrun2','vlq']:
      if '*' in hist_name:
        # add the usual pT binned uncerts for et and mt with pT<100 
        if options.channel in ['et','mt']:

          pt_bins = ["20-25", "25-30", "30-35", "35-40", "40-500", "500-1000", "1000-inf"]

          for i in range(3,6):
            bin_name = pt_bins[i-1]
            hist_name_bini = hist_name.replace('*','%(bin_name)s' % vars())
            systematics['syst_tau_id_diff_bin%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*(wt_tau_id_pt_bin%(i)i_up*(pt_2<100) + (pt_2>=100))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
            systematics['syst_tau_id_diff_bin%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*(wt_tau_id_pt_bin%(i)i_down*(pt_2<100) + (pt_2>=100))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)

        if options.channel == 'tt':
          # add the usual dm binned uncerts for tt with pt<100
          for i in [0,1,10,11]:
            hist_name_bini = hist_name.replace('*','dm%(i)i' % vars())
            systematics['syst_tau_id_diff_dm%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_id_dm%(i)i_up' % vars(), ['QCD','jetFakes'], False)
            systematics['syst_tau_id_diff_dm%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_id_dm%(i)i_down'% vars(), ['QCD','jetFakes'], False)

        # add the dedicated high pT uncerts which are common for both channels
   
        pt_bins = {5: "highpT_100-500",  6: "highpT_500-inf"}

        # weights stored incorrectly, doing by hand instead
        #for i in range(5,7):
          #bin_name = pt_bins[i]
          #hist_name_bini = hist_name.replace('*','%(bin_name)s' % vars())
          #systematics['syst_tau_id_diff_highpt_bin%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_id_mssm_bin%(i)i_up' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
          #systematics['syst_tau_id_diff_highpt_bin%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_id_mssm_bin%(i)i_down' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
       
        if options.year == "2018":
          nom_40_plus = "1.00458784525"
          up_40_500, up_500_1000_const,up_500_1000_scale, up_1000_plus = "1.03510518813", "1.00458784525", "0.0305173428802", "0.0610346857604"
          down_40_500, down_500_1000_const,down_500_1000_scale, down_1000_plus = "0.952971628935", "1.00458784525", "-0.0516162163153", "-0.103232432631"
        elif options.year == "2017":
          nom_40_plus = "0.86966080714"
          up_40_500, up_500_1000_const,up_500_1000_scale, up_1000_plus = "0.903623575349", "0.86966080714", "0.0339627682083", "0.0679255364167"
          down_40_500, down_500_1000_const,down_500_1000_scale, down_1000_plus = "0.835685230027", "0.86966080714", "-0.0339755771129", "-0.0679511542259"
        elif options.year == "2016":
          nom_40_plus = "0.949771847264"
          up_40_500, up_500_1000_const,up_500_1000_scale, up_1000_plus = "0.982473629072", "0.949771847264", "0.0327017818078", "0.0654035636156"
          down_40_500, down_500_1000_const,down_500_1000_scale, down_1000_plus = "0.908265091132", "0.949771847264", "-0.0415067561324", "-0.0830135122647"

        bin5_up_formula = "(((gen_match_X==5&&pt_X>100&&pt_X<=500)*(%(up_40_500)s/%(nom_40_plus)s)) + ((gen_match_X==5&&pt_X>100&&pt_X<=500)==0))" % vars()
        bin5_down_formula = "(((gen_match_X==5&&pt_X>100&&pt_X<=500)*(%(down_40_500)s/%(nom_40_plus)s)) + ((gen_match_X==5&&pt_X>100&&pt_X<=500)==0))" % vars()
        bin6_up_formula = "((((gen_match_X==5&&pt_X>500&&pt_X<=1000)*((%(up_500_1000_const)s + %(up_500_1000_scale)s*(pt_X/500.))/%(nom_40_plus)s))+((gen_match_X==5&&pt_X>1000)*((%(up_500_1000_const)s + %(up_1000_plus)s)/%(nom_40_plus)s)))+((gen_match_X==5&&pt_X>500)==0))" % vars()
        bin6_down_formula = "((((gen_match_X==5&&pt_X>500&&pt_X<=1000)*((%(down_500_1000_const)s + %(down_500_1000_scale)s*(pt_X/500.))/%(nom_40_plus)s))+((gen_match_X==5&&pt_X>1000)*((%(down_500_1000_const)s + %(down_1000_plus)s)/%(nom_40_plus)s)))+((gen_match_X==5&&pt_X>500)==0))" % vars()


        if options.channel in ["mt","et"]:
          systematics['syst_tau_id_diff_highpt_bin5_up' % vars()] = ('' , '_'+hist_name.replace('*','highpT_100-500')+'Up', 'wt*{}'.format(bin5_up_formula.replace("X","2")), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
          systematics['syst_tau_id_diff_highpt_bin5_down' % vars()] = ('' , '_'+hist_name.replace('*','highpT_100-500')+'Down', 'wt*{}'.format(bin5_down_formula.replace("X","2")), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
          systematics['syst_tau_id_diff_highpt_bin6_up' % vars()] = ('' , '_'+hist_name.replace('*','highpT_500-inf')+'Up', 'wt*{}'.format(bin6_up_formula.replace("X","2")), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
          systematics['syst_tau_id_diff_highpt_bin6_down' % vars()] = ('' , '_'+hist_name.replace('*','highpT_500-inf')+'Down', 'wt*{}'.format(bin6_down_formula.replace("X","2")), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
        elif options.channel == "tt":
          systematics['syst_tau_id_diff_highpt_bin5_up' % vars()] = ('' , '_'+hist_name.replace('*','highpT_100-500')+'Up', 'wt*{}*{}'.format(bin5_up_formula.replace("X","1"),bin5_up_formula.replace("X","2")), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
          systematics['syst_tau_id_diff_highpt_bin5_down' % vars()] = ('' , '_'+hist_name.replace('*','highpT_100-500')+'Down', 'wt*{}*{}'.format(bin5_down_formula.replace("X","1"),bin5_down_formula.replace("X","2")), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
          systematics['syst_tau_id_diff_highpt_bin6_up' % vars()] = ('' , '_'+hist_name.replace('*','highpT_500-inf')+'Up', 'wt*{}*{}'.format(bin6_up_formula.replace("X","1"),bin6_up_formula.replace("X","2")), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
          systematics['syst_tau_id_diff_highpt_bin6_down' % vars()] = ('' , '_'+hist_name.replace('*','highpT_500-inf')+'Down', 'wt*{}*{}'.format(bin6_down_formula.replace("X","1"),bin6_down_formula.replace("X","2")), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)

        


        hist_name_pi0 = hist_name.replace('_eff_t','').replace('*','1ProngPi0Eff' % vars())
        hist_name_pi = hist_name.replace('_eff_t','').replace('*','3ProngEff' % vars())
        if options.channel in ['et','mt']:
          systematics['syst_tau_id_diff_emb_trk_pi_up' % vars()] = ('' , '_'+hist_name_pi+'Up', 'wt*((pt_2>=100)+(pt_2<100)*((tau_decay_mode_2<3)*1.008+(tau_decay_mode_2>9)*1.024))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes','ZTT','TTT','VVT','signal'], False)
          systematics['syst_tau_id_diff_emb_trk_pi_down' % vars()] = ('' , '_'+hist_name_pi+'Down', 'wt*((pt_2>=100)+(pt_2<100)*((tau_decay_mode_2<3)*0.992+(tau_decay_mode_2>9)*0.976))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes','ZTT','TTT','VVT','signal'], False)
          systematics['syst_tau_id_diff_emb_trk_pi0_up' % vars()] = ('' , '_'+hist_name_pi0+'Up', 'wt*((pt_2>=100)+(pt_2<100)*((tau_decay_mode_2==0 || tau_decay_mode_2==10) + ((tau_decay_mode_2==0 || tau_decay_mode_2==10)==0)*1.014))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes','ZTT','TTT','VVT','signal'], False) 
          systematics['syst_tau_id_diff_emb_trk_pi0_down' % vars()] = ('' , '_'+hist_name_pi0+'Down', 'wt*((pt_2>=100)+(pt_2<100)*((tau_decay_mode_2==0 || tau_decay_mode_2==10) + ((tau_decay_mode_2==0 || tau_decay_mode_2==10)==0)*0.986))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes','ZTT','TTT','VVT','signal'], False) 
        if options.channel in ['tt']:
          systematics['syst_tau_id_diff_emb_trk_pi_up' % vars()] = ('' , '_'+hist_name_pi+'Up', 'wt*((pt_2>=100)+(pt_2<100)*((tau_decay_mode_2<3)*1.008+(tau_decay_mode_2>9)*1.024))*((pt_1>=100)+(pt_1<100)*((tau_decay_mode_1<3)*1.008+(tau_decay_mode_1>9)*1.024))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes','ZTT','TTT','VVT','signal'], False)
          systematics['syst_tau_id_diff_emb_trk_pi_down' % vars()] = ('' , '_'+hist_name_pi+'Down', 'wt*((pt_2>=100)+(pt_2<100)*((tau_decay_mode_2<3)*0.992+(tau_decay_mode_2>9)*0.976))*((pt_1>=100)+(pt_1<100)*((tau_decay_mode_1<3)*0.992+(tau_decay_mode_1>9)*0.976))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes','ZTT','TTT','VVT','signal'], False)
          systematics['syst_tau_id_diff_emb_trk_pi0_up' % vars()] = ('' , '_'+hist_name_pi0+'Up', 'wt*((pt_2>=100)+(pt_2<100)*((tau_decay_mode_2==0 || tau_decay_mode_2==10) + ((tau_decay_mode_2==0 || tau_decay_mode_2==10)==0)*1.014))*((pt_1>=100)+(pt_1<100)*((tau_decay_mode_1==0 || tau_decay_mode_1==10) + ((tau_decay_mode_1==0 || tau_decay_mode_1==10)==0)*1.014))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes','ZTT','TTT','VVT','signal'], False)
          systematics['syst_tau_id_diff_emb_trk_pi0_down' % vars()] = ('' , '_'+hist_name_pi0+'Down', 'wt*((pt_2>=100)+(pt_2<100)*((tau_decay_mode_2==0 || tau_decay_mode_2==10) + ((tau_decay_mode_2==0 || tau_decay_mode_2==10)==0)*0.986))*((pt_1>=100)+(pt_1<100)*((tau_decay_mode_1==0 || tau_decay_mode_1==10) + ((tau_decay_mode_1==0 || tau_decay_mode_1==10)==0)*0.986))' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes','ZTT','TTT','VVT','signal'], False) 

 
    elif '*PT' in hist_name:
      for i in range(1,6):
        hist_name_bini = hist_name.replace('*PT','bin%(i)i' % vars())
        systematics['syst_tau_id_diff_bin%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_id_pt_bin%(i)i_up' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
        systematics['syst_tau_id_diff_bin%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_id_pt_bin%(i)i_down' % vars(), ['ZLL','TTJ','ZL','ZJ','VVJ','W','jetFakes'], False)
    elif '*DM' in hist_name:
      for i in [0,1,10,11]:
        hist_name_bini = hist_name.replace('*DM','DM%(i)i' % vars())
        systematics['syst_tau_id_diff_dm%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_id_dm%(i)i_up' % vars(), ['QCD','jetFakes'], False)
        systematics['syst_tau_id_diff_dm%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_id_dm%(i)i_down'% vars(), ['QCD','jetFakes'], False) 
    elif '*MVADM' in hist_name:
      for i in [0,1,2,10,11]:
        if options.channel != 'tt':
          hist_name_bini = hist_name.replace('*MVADM','pTlow_MVADM%(i)i' % vars())
          systematics['syst_tau_id_diff_lowpt_mvadm%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_id_lowpt_mvadm%(i)i_up' % vars(), ['QCD','jetFakes'], False)
          systematics['syst_tau_id_diff_lowpt_mvadm%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_id_lowpt_mvadm%(i)i_down'% vars(), ['QCD','jetFakes'], False)

        hist_name_bini = hist_name.replace('*MVADM','pThigh_MVADM%(i)i' % vars())
        systematics['syst_tau_id_diff_highpt_mvadm%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_id_highpt_mvadm%(i)i_up' % vars(), ['QCD','jetFakes'], False)
        systematics['syst_tau_id_diff_highpt_mvadm%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_id_highpt_mvadm%(i)i_down'% vars(), ['QCD','jetFakes'], False)

if options.syst_lep_trg_diff != '':
  hist_name = options.syst_lep_trg_diff
  chan = options.channel
  if options.analysis in ['mssmrun2','vlq']:
    hist_name_bini_slt = hist_name.replace('*','trigger_%(chan)s' % vars())
    hist_name_bini_ct = hist_name.replace('*','xtrigger_l_%(chan)s' % vars())
    if chan in ["et"]:
      systematics['syst_lep_trg_diff_singlelep_up' % vars()] = ('' , '_'+hist_name_bini_slt+'Up', 'wt*((1.02*(trg_singleelectron==1)) + (trg_singleelectron==0))' % vars(), ['QCD','jetFakes'], False)
      systematics['syst_lep_trg_diff_singlelep_down' % vars()] = ('' , '_'+hist_name_bini_slt+'Down', 'wt*((1.02*(trg_singleelectron==1)) + (trg_singleelectron==0))'% vars(), ['QCD','jetFakes'], False)
      systematics['syst_lep_trg_diff_crosstrg_up' % vars()] = ('' , '_'+hist_name_bini_ct+'Up', 'wt*((1.02*(trg_etaucross==1)) + (trg_etaucross==0))' % vars(), ['QCD','jetFakes'], False)
      systematics['syst_lep_trg_diff_crosstrg_down' % vars()] = ('' , '_'+hist_name_bini_ct+'Down', 'wt*((1.02*(trg_etaucross==1)) + (trg_etaucross==0))'% vars(), ['QCD','jetFakes'], False)
    elif chan in ["mt"]:
      systematics['syst_lep_trg_diff_singlelep_up' % vars()] = ('' , '_'+hist_name_bini_slt+'Up', 'wt*((1.02*(trg_singlemuon==1)) + (trg_singlemuon==0))' % vars(), ['QCD','jetFakes'], False)
      systematics['syst_lep_trg_diff_singlelep_down' % vars()] = ('' , '_'+hist_name_bini_slt+'Down', 'wt*((1.02*(trg_singlemuon==1)) + (trg_singlemuon==0))'% vars(), ['QCD','jetFakes'], False)
      systematics['syst_lep_trg_diff_crosstrg_up' % vars()] = ('' , '_'+hist_name_bini_ct+'Up', 'wt*((1.02*(trg_mutaucross==1)) + (trg_mutaucross==0))' % vars(), ['QCD','jetFakes'], False)
      systematics['syst_lep_trg_diff_crosstrg_down' % vars()] = ('' , '_'+hist_name_bini_ct+'Down', 'wt*((1.02*(trg_mutaucross==1)) + (trg_mutaucross==0))'% vars(), ['QCD','jetFakes'], False)


if options.syst_tau_trg_diff != '':
    hist_name = options.syst_tau_trg_diff
    chan = options.channel

    if options.analysis in ['mssmrun2','vlq']:
      if '*'in hist_name:

        # single tau uncerts
        hist_name_bini = hist_name.replace('*','trigger_single_t' % vars())
        systematics['syst_tau_trg_diff_singletau_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_trg_mssm_singletau_up' % vars(), ['QCD','jetFakes'], False)
        systematics['syst_tau_trg_diff_singletau_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_trg_mssm_singletau_down'% vars(), ['QCD','jetFakes'], False)

        if chan in ['mt','et']:

          # single lepton uncerts
          hist_name_bini = hist_name.replace('*','trigger_%(chan)s' % vars())
          systematics['syst_tau_trg_diff_singlelep_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_trg_mssm_singlelep_up' % vars(), ['QCD','jetFakes'], False)
          systematics['syst_tau_trg_diff_singlelep_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_trg_mssm_singlelep_down'% vars(), ['QCD','jetFakes'], False)

          # cross lepton uncerts
          hist_name_bini = hist_name.replace('*','xtrigger_l_%(chan)s' % vars())
          systematics['syst_tau_trg_diff_crosslep_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_trg_mssm_crosslep_up' % vars(), ['QCD','jetFakes'], False)
          systematics['syst_tau_trg_diff_crosslep_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_trg_mssm_crosslep_down'% vars(), ['QCD','jetFakes'], False)

        # double-tau / cross-tau uncerts
        for i in [0,1,10,11]:

          hist_name_bini = hist_name.replace('*','xtrigger_t_%(chan)s_dm%(i)i' % vars())

          systematics['syst_tau_trg_diff_dm%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_trg_mssm_dm%(i)i_up' % vars(), ['QCD','jetFakes'], False)
          systematics['syst_tau_trg_diff_dm%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_trg_mssm_dm%(i)i_down'% vars(), ['QCD','jetFakes'], False)

          if chan == 'tt':
            # for tt channels we decouple low and high pT so add additional uncertainty for high pT here
            hist_name_bini = hist_name.replace('*','xtrigger_t_%(chan)s_dm%(i)i_highpT' % vars())
            systematics['syst_tau_trg_diff_dm%(i)i_highpt_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_trg_mssm_highpt_dm%(i)i_up' % vars(), ['QCD','jetFakes'], False)
            systematics['syst_tau_trg_diff_dm%(i)i_highpt_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_trg_mssm_highpt_dm%(i)i_down'% vars(), ['QCD','jetFakes'], False)



    elif '*DM' in hist_name:
      for i in [0,1,10,11]:
        hist_name_bini = hist_name.replace('*DM','DM%(i)i' % vars())
        if options.channel != 'tt':
          systematics['syst_tau_trg_diff_dm%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*((wt_tau_trg_dm%(i)i_up==1) + (wt_tau_trg_dm%(i)i_up!=1)*(max(1.05,wt_tau_trg_dm%(i)i_up)))' % vars(), ['QCD','jetFakes'], False)
          systematics['syst_tau_trg_diff_dm%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*((wt_tau_trg_dm%(i)i_down==1) + (wt_tau_trg_dm%(i)i_down!=1)*(min(0.95,wt_tau_trg_dm%(i)i_down)))' % vars(), ['QCD','jetFakes'], False)
        else: 
          systematics['syst_tau_trg_diff_dm%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_trg_dm%(i)i_up' % vars(), ['QCD','jetFakes'], False)
          systematics['syst_tau_trg_diff_dm%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_trg_dm%(i)i_down'% vars(), ['QCD','jetFakes'], False)
    elif '*MVADM' in hist_name:

      for i in [0,1,2,10,11]:
        hist_name_bini = hist_name.replace('*MVADM','MVADM%(i)i' % vars())
        systematics['syst_tau_trg_diff_mvadm%(i)i_up' % vars()] = ('' , '_'+hist_name_bini+'Up', 'wt*wt_tau_trg_mvadm%(i)i_up' % vars(), ['QCD','jetFakes'], False)
        systematics['syst_tau_trg_diff_mvadm%(i)i_down' % vars()] = ('' , '_'+hist_name_bini+'Down', 'wt*wt_tau_trg_mvadm%(i)i_down'% vars(), ['QCD','jetFakes'], False)

if options.syst_embed_pt != '':
    hist_name = options.syst_embed_pt
    systematics['syst_embed_pt_up' % vars()] = ('' , '_'+hist_name+'Up', 'wt*wt_zpt_embed_ic' % vars(), ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)
    systematics['syst_embed_pt_down' % vars()] = ('' , '_'+hist_name+'Down', 'wt*(2-wt_zpt_embed_ic)'% vars(), ['QCD','jetFakes','EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal'], False)

if options.method in [17,18] and options.do_ff_systs and options.channel in ['et','mt','tt'] and options.era == 'mssmsummer16':
    processes = ['tt','w','qcd']
    dms = ['dm0', 'dm1']
    njets = ['njet0','njet1']
    for process in processes:
      template_name = 'ff_'+process+'_syst'
      if process is 'qcd' or options.channel == 'tt': template_name = 'ff_'+process+'_'+options.channel+'_syst'
      weight_name = 'wt_ff_'+options.cat+'_'+process+'_syst_'
      systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
      systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
      if options.channel == 'tt' and process in ['w','tt']: continue
      for dm in dms: 
        for njet in njets:
          template_name = 'ff_'+process+'_'+dm+'_'+njet
          if process != 'tt': template_name+='_'+options.channel
          template_name+='_stat'
          weight_name = 'wt_ff_'+options.cat+'_'+process+'_'+dm+'_'+njet+'_stat_'
          systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
          systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
    if options.channel == "tt":
      processes = ['dy', 'w', 'tt']
      for process in processes:
        template_name = 'ff_'+process+'_frac_tt_syst'
        weight_name = 'wt_ff_'+options.cat+'_'+process+'_frac_syst_'
        systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
        systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)

if options.method in [17,18] and options.channel in ['et','mt','tt'] and options.analysis=='cpdecay' and options.do_ff_systs:
  if options.channel in ['et','mt']:
    chan = options.channel
    lt_systs = {}
    lt_systs['ff_%(chan)s_qcd_syst' % vars()] = 'wt_ff_qcd_syst_' % vars()
    lt_systs['ff_%(chan)s_wjets_syst' % vars()] = 'wt_ff_wjets_syst_' % vars()
    lt_systs['ff_%(chan)s_ttbar_syst' % vars()] = 'wt_ff_ttbar_syst_' % vars()

    tt_systs['ff_%(chan)s_qcd_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_1+(n_jets==0)*wt_ff_qcd_syst_' % vars()
    tt_systs['ff_%(chan)s_qcd_syst_njets1' % vars()] = '(n_jets==0)*wt_ff_1+(n_jets>=1)*wt_ff_qcd_syst_' % vars()

    tt_systs['ff_%(chan)s_wjets_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_1+(n_jets==0)*wt_ff_wjets_syst_' % vars()
    tt_systs['ff_%(chan)s_wjets_syst_njets1' % vars()] = '(n_jets!=1)*wt_ff_1+(n_jets==1)*wt_ff_wjets_syst_' % vars()
    tt_systs['ff_%(chan)s_wjets_syst_njets2' % vars()] = '(n_jets<2)*wt_ff_1+(n_jets>=2)*wt_ff_wjets_syst_' % vars()


    tt_systs['ff_%(chan)s_wjets_met_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_1+(n_jets==0)*wt_ff_wjets_met_' % vars()
    tt_systs['ff_%(chan)s_wjets_met_closure_syst_njets1' % vars()] = '(n_jets!=1)*wt_ff_1+(n_jets==1)*wt_ff_wjets_met_' % vars()
    tt_systs['ff_%(chan)s_wjets_met_closure_syst_njets2' % vars()] = '(n_jets<2)*wt_ff_1+(n_jets>=2)*wt_ff_wjets_met_' % vars()

    tt_systs['ff_%(chan)s_qcd_met_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_1+(n_jets==0)*wt_ff_qcd_met_' % vars()
    tt_systs['ff_%(chan)s_qcd_met_closure_syst_njets1' % vars()] = '(n_jets==0)*wt_ff_1+(n_jets>=1)*wt_ff_qcd_met_' % vars()

    tt_systs['ff_%(chan)s_wjets_l_pt_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_1+(n_jets==0)*wt_ff_wjets_l_pt_' % vars()
    tt_systs['ff_%(chan)s_wjets_l_pt_closure_syst_njets1' % vars()] = '(n_jets!=1)*wt_ff_1+(n_jets==1)*wt_ff_wjets_l_pt_' % vars()
    tt_systs['ff_%(chan)s_wjets_l_pt_closure_syst_njets2' % vars()] = '(n_jets<2)*wt_ff_1+(n_jets>=2)*wt_ff_wjets_l_pt_' % vars()

    tt_systs['ff_%(chan)s_qcd_l_pt_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_1+(n_jets==0)*wt_ff_qcd_l_pt_' % vars()
    tt_systs['ff_%(chan)s_qcd_l_pt_closure_syst_njets1' % vars()] = '(n_jets==0)*wt_ff_1+(n_jets>=1)*wt_ff_qcd_l_pt_' % vars()

    for proc in ['qcd','wjets']:
      lt_systs['ff_%(chan)s_%(proc)s_met_closure_syst' % vars()] = 'wt_ff_%(proc)s_met_' % vars()
      lt_systs['ff_%(chan)s_%(proc)s_l_pt_closure_syst' % vars()] = 'wt_ff_%(proc)s_l_pt_' % vars()

      for njet in ['0','1','2']:
        for dm in ['0_sig_lt3','0_sig_gt3','1','2','10','11']:
          #lt_systs[('ff_%(chan)s_%(proc)s_stat_njets%(njet)s_mvadm%(dm)s' % vars()).replace('lt3','lt').replace('gt3','gt')] = 'wt_ff_%(proc)s_stat_njet%(njet)s_mvadm%(dm)s_' % vars()
          lt_systs[('ff_%(chan)s_%(proc)s_stat_unc1_njets%(njet)s_mvadm%(dm)s' % vars()).replace('lt3','lt').replace('gt3','gt')] = 'wt_ff_%(proc)s_stat_unc1_njet%(njet)s_mvadm%(dm)s_' % vars()
          lt_systs[('ff_%(chan)s_%(proc)s_stat_unc2_njets%(njet)s_mvadm%(dm)s' % vars()).replace('lt3','lt').replace('gt3','gt')] = 'wt_ff_%(proc)s_stat_unc2_njet%(njet)s_mvadm%(dm)s_' % vars()
    for template_name in lt_systs:
      weight_name = lt_systs[template_name]
      systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
      systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)

  if options.channel in ['tt']:
    tt_systs={}
    tt_systs['ff_tt_qcd_met_closure_syst' % vars()] = 'wt_ff_qcd_met_' % vars()
    tt_systs['ff_tt_qcd_met_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_1+(n_jets==0)*wt_ff_qcd_met_' % vars()
    tt_systs['ff_tt_qcd_met_closure_syst_njets1' % vars()] = '(n_jets!=1)*wt_ff_1+(n_jets==1)*wt_ff_qcd_met_' % vars()
    tt_systs['ff_tt_qcd_met_closure_syst_njets2' % vars()] = '(n_jets<2)*wt_ff_1+(n_jets>=2)*wt_ff_qcd_met_' % vars()
    tt_systs['ff_tt_qcd_syst' % vars()] = 'wt_ff_qcd_syst_' % vars()
    tt_systs['ff_tt_qcd_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_1+(n_jets==0)*wt_ff_qcd_syst_' % vars()
    tt_systs['ff_tt_qcd_syst_njets1' % vars()] = '(n_jets==0)*wt_ff_1+(n_jets>=1)*wt_ff_qcd_syst_' % vars()

    tt_systs['ff_tt_wjets_syst' % vars()] = 'wt_ff_wjets_syst_' % vars()
    for njet in ['0','1','2']:
      for dm in ['0_sig_lt3','0_sig_gt3','1','2','10','11']:
        #tt_systs[('ff_tt_qcd_stat_njets%(njet)s_mvadm%(dm)s' % vars()).replace('lt3','lt').replace('gt3','gt')] = 'wt_ff_qcd_stat_njet%(njet)s_mvadm%(dm)s_' % vars()
        tt_systs[('ff_tt_qcd_stat_unc1_njets%(njet)s_mvadm%(dm)s' % vars()).replace('lt3','lt').replace('gt3','gt')] = 'wt_ff_qcd_stat_unc1_njet%(njet)s_mvadm%(dm)s_' % vars()
        tt_systs[('ff_tt_qcd_stat_unc2_njets%(njet)s_mvadm%(dm)s' % vars()).replace('lt3','lt').replace('gt3','gt')] = 'wt_ff_qcd_stat_unc2_njet%(njet)s_mvadm%(dm)s_' % vars()

    for template_name in tt_systs:
      weight_name = tt_systs[template_name]
      systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
      systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)

  template_name = 'ff_%s_sub_syst' % (options.channel)
  systematics['ff_sub_up']   = ('' , '_'+template_name+'Up',   'wt_ff',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
  systematics['ff_sub_down'] = ('' , '_'+template_name+'Down', 'wt_ff', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)

if options.method in [17,18] and options.channel in ['et','mt','tt'] and options.analysis=='cpprod' and options.do_ff_systs:
  if options.channel in ['et','mt']:
    chan = options.channel
    lt_systs = {}
    lt_systs['ff_%(chan)s_qcd_syst' % vars()] = 'wt_ff_dmbins_qcd_syst_' % vars()
    lt_systs['ff_%(chan)s_wjets_syst' % vars()] = 'wt_ff_dmbins_wjets_syst_' % vars()
    lt_systs['ff_%(chan)s_ttbar_syst' % vars()] = 'wt_ff_dmbins_ttbar_syst_' % vars()

    lt_systs['ff_%(chan)s_qcd_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_dmbins_1+(n_jets==0)*wt_ff_dmbins_qcd_syst_' % vars()
    lt_systs['ff_%(chan)s_qcd_syst_njets1' % vars()] = '(n_jets==0)*wt_ff_dmbins_1+(n_jets>=1)*wt_ff_dmbins_qcd_syst_' % vars()

    lt_systs['ff_%(chan)s_wjets_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_dmbins_1+(n_jets==0)*wt_ff_dmbins_wjets_syst_' % vars()
    lt_systs['ff_%(chan)s_wjets_syst_njets1' % vars()] = '(n_jets!=1)*wt_ff_dmbins_1+(n_jets==1)*wt_ff_dmbins_wjets_syst_' % vars()
    lt_systs['ff_%(chan)s_wjets_syst_njets2' % vars()] = '(n_jets<2)*wt_ff_dmbins_1+(n_jets>=2)*wt_ff_dmbins_wjets_syst_' % vars()


    lt_systs['ff_%(chan)s_wjets_met_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_dmbins_1+(n_jets==0)*wt_ff_dmbins_wjets_met_' % vars()
    lt_systs['ff_%(chan)s_wjets_met_closure_syst_njets1' % vars()] = '(n_jets!=1)*wt_ff_dmbins_1+(n_jets==1)*wt_ff_dmbins_wjets_met_' % vars()
    lt_systs['ff_%(chan)s_wjets_met_closure_syst_njets2' % vars()] = '(n_jets<2)*wt_ff_dmbins_1+(n_jets>=2)*wt_ff_dmbins_wjets_met_' % vars()

    lt_systs['ff_%(chan)s_qcd_met_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_dmbins_1+(n_jets==0)*wt_ff_dmbins_qcd_met_' % vars()
    lt_systs['ff_%(chan)s_qcd_met_closure_syst_njets1' % vars()] = '(n_jets==0)*wt_ff_dmbins_1+(n_jets>=1)*wt_ff_dmbins_qcd_met_' % vars()

    lt_systs['ff_%(chan)s_wjets_l_pt_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_dmbins_1+(n_jets==0)*wt_ff_dmbins_wjets_l_pt_' % vars()
    lt_systs['ff_%(chan)s_wjets_l_pt_closure_syst_njets1' % vars()] = '(n_jets!=1)*wt_ff_dmbins_1+(n_jets==1)*wt_ff_dmbins_wjets_l_pt_' % vars()
    lt_systs['ff_%(chan)s_wjets_l_pt_closure_syst_njets2' % vars()] = '(n_jets<2)*wt_ff_dmbins_1+(n_jets>=2)*wt_ff_dmbins_wjets_l_pt_' % vars()

    lt_systs['ff_%(chan)s_qcd_l_pt_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_dmbins_1+(n_jets==0)*wt_ff_dmbins_qcd_l_pt_' % vars()
    lt_systs['ff_%(chan)s_qcd_l_pt_closure_syst_njets1' % vars()] = '(n_jets==0)*wt_ff_dmbins_1+(n_jets>=1)*wt_ff_dmbins_qcd_l_pt_' % vars()

    for proc in ['qcd','wjets']:
      lt_systs['ff_%(chan)s_%(proc)s_met_closure_syst' % vars()] = 'wt_ff_dmbins_%(proc)s_met_' % vars()
      lt_systs['ff_%(chan)s_%(proc)s_l_pt_closure_syst' % vars()] = 'wt_ff_dmbins_%(proc)s_l_pt_' % vars()
      for njet in ['0','1','2']:
        for dm in ['0','1','10','11']:
          lt_systs[('ff_%(chan)s_%(proc)s_stat_unc1_njets%(njet)s_dm%(dm)s' % vars())] = 'wt_ff_dmbins_%(proc)s_stat_unc1_njet%(njet)s_dm%(dm)s_' % vars()
          lt_systs[('ff_%(chan)s_%(proc)s_stat_unc2_njets%(njet)s_dm%(dm)s' % vars())] = 'wt_ff_dmbins_%(proc)s_stat_unc2_njet%(njet)s_dm%(dm)s_' % vars()
    for template_name in lt_systs:
      weight_name = lt_systs[template_name]
      systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
      systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)

  if options.channel in ['tt']:
    tt_systs={}
    tt_systs['ff_tt_qcd_met_closure_syst' % vars()] = 'wt_ff_dmbins_qcd_met_' % vars()
    tt_systs['ff_tt_qcd_met_closure_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_dmbins_1+(n_jets==0)*wt_ff_dmbins_qcd_met_' % vars()
    tt_systs['ff_tt_qcd_met_closure_syst_njets1' % vars()] = '(n_jets!=1)*wt_ff_dmbins_1+(n_jets==1)*wt_ff_dmbins_qcd_met_' % vars()
    tt_systs['ff_tt_qcd_met_closure_syst_njets2' % vars()] = '(n_jets<2)*wt_ff_dmbins_1+(n_jets>=2)*wt_ff_dmbins_qcd_met_' % vars()

#    tt_systs['ff_tt_qcd_syst' % vars()] = 'wt_ff_dmbins_qcd_syst_' % vars()
    tt_systs['ff_tt_qcd_syst_njets0' % vars()] = '(n_jets!=0)*wt_ff_dmbins_1+(n_jets==0)*wt_ff_dmbins_qcd_syst_' % vars()
    tt_systs['ff_tt_qcd_syst_njets1' % vars()] = '(n_jets==0)*wt_ff_dmbins_1+(n_jets>=1)*wt_ff_dmbins_qcd_syst_' % vars()

    tt_systs['ff_tt_wjets_syst' % vars()] = 'wt_ff_dmbins_wjets_syst_' % vars()
    for njet in ['0','1','2']:
      for dm in ['0','1','10','11']:
        tt_systs[('ff_tt_qcd_stat_unc1_njets%(njet)s_dm%(dm)s' % vars())] = 'wt_ff_dmbins_qcd_stat_unc1_njet%(njet)s_dm%(dm)s_' % vars()
        tt_systs[('ff_tt_qcd_stat_unc2_njets%(njet)s_dm%(dm)s' % vars())] = 'wt_ff_dmbins_qcd_stat_unc2_njet%(njet)s_dm%(dm)s_' % vars()

    for template_name in tt_systs:
      weight_name = tt_systs[template_name]
      systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
      systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)

  template_name = 'ff_%s_sub_syst' % (options.channel)
  systematics['ff_sub_up']   = ('' , '_'+template_name+'Up',   'wt_ff',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)
  systematics['ff_sub_down'] = ('' , '_'+template_name+'Down', 'wt_ff', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT'], True)

if options.method in [17,18] and options.channel in ['et','mt','tt'] and options.analysis in ['mssmrun2','vlq'] and options.do_ff_systs:
  ch = options.channel
  yr = options.year
  if options.channel in ['tt']:
    tt_systs={}
    
    for i in [1,2]:
      tt_systs['CMS_ff_total_qcd_stat_dR_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_stat_dR_unc%(i)i_' % vars()
      tt_systs['CMS_ff_total_qcd_stat_pt_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_stat_pt_unc%(i)i_' % vars()
    tt_systs['CMS_ff_total_qcd_syst_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_syst_' % vars()
    tt_systs['CMS_ff_total_wjets_syst_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_syst_' % vars()
    tt_systs['CMS_ff_total_ttbar_syst_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_ttbar_syst_' % vars()

    # new additional uncertainties
    tt_systs['CMS_ff_total_qcd_syst_dr_closure_%(ch)s_%(yr)s' % vars()]  = 'wt_ff_mssm_qcd_syst_dR_closure_' % vars()
    tt_systs['CMS_ff_total_qcd_syst_met_closure_%(ch)s_%(yr)s' % vars()]  = 'wt_ff_mssm_qcd_syst_met_closure_' % vars()
    tt_systs['CMS_ff_total_qcd_syst_pt_2_closure_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_syst_pt_2_closure_' % vars()
    tt_systs['CMS_ff_total_syst_alt_func_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_syst_alt_func_' % vars()

    for njet in [0,1]:
      for jetpt in ['low','med','high']:
        for i in [1,2,3]:
          tt_systs[('CMS_ff_total_qcd_stat_njet%(njet)i_jet_pt_%(jetpt)s_unc%(i)i_%(ch)s_%(yr)s' % vars())] = 'wt_ff_mssm_qcd_stat_njet%(njet)i_jet_pt_%(jetpt)s_unc%(i)i_' % vars()


    for template_name in tt_systs:
      weight_name = tt_systs[template_name]
      systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT','ZLL','TT','VV'], True)
      systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT','ZLL','TT','VV'], True)

  elif options.channel in ['et','mt']:
    lt_systs={}

    for njet in [0,1]:
      for i in [1,2]:
        lt_systs['CMS_ff_total_qcd_stat_ss_njets%(njet)i_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_stat_ss_njets%(njet)i_unc%(i)i_' % vars()
        lt_systs['CMS_ff_total_wjets_stat_met_njets%(njet)i_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_stat_met_njets%(njet)i_unc%(i)i_' % vars()
        lt_systs['CMS_ff_total_wjets_stat_l_pt_njets%(njet)i_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_stat_l_pt_njets%(njet)i_unc%(i)i_' % vars()

    for i in [1,2]:
      lt_systs['CMS_ff_total_qcd_stat_os_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_stat_os_unc%(i)i_' % vars()
      lt_systs['CMS_ff_total_qcd_stat_l_pt_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_stat_l_pt_unc%(i)i_' % vars()
      lt_systs['CMS_ff_total_qcd_stat_iso_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_stat_iso_unc%(i)i_' % vars()
      lt_systs['CMS_ff_total_wjets_stat_extrap_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_stat_extrap_unc%(i)i_' % vars()
      lt_systs['CMS_ff_total_ttbar_stat_met_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_ttbar_stat_met_unc%(i)i_' % vars()
      lt_systs['CMS_ff_total_ttbar_stat_l_pt_unc%(i)i_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_ttbar_stat_l_pt_unc%(i)i_' % vars()

    lt_systs['CMS_ff_total_qcd_syst_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_syst_' % vars()
    lt_systs['CMS_ff_total_qcd_syst_iso_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_syst_iso_' % vars()
    lt_systs['CMS_ff_total_wjets_syst_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_syst_' % vars()
    lt_systs['CMS_ff_total_wjets_syst_extrap_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_syst_extrap_' % vars()
    lt_systs['CMS_ff_total_ttbar_syst_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_ttbar_syst_' % vars()

    lt_systs['CMS_ff_total_wjets_frac_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_frac_' % vars()
    lt_systs['CMS_ff_total_ttbar_frac_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_ttbar_frac_' % vars()
  
    lt_systs['CMS_ff_total_low_pt_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_low_pt_' % vars()

    # new additional uncertainties
    lt_systs['CMS_ff_total_qcd_syst_met_closure_%(ch)s_%(yr)s' % vars()]  = 'wt_ff_mssm_qcd_syst_met_closure_' % vars()
    lt_systs['CMS_ff_total_wjets_syst_met_closure_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_syst_met_closure_' % vars()
    lt_systs['CMS_ff_total_ttbar_syst_met_closure_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_ttbar_syst_met_closure_' % vars()
    lt_systs['CMS_ff_total_wjets_syst_l_pt_closure_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_syst_l_pt_closure_' % vars()
    lt_systs['CMS_ff_total_ttbar_syst_l_pt_closure_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_ttbar_syst_l_pt_closure_' % vars()
    lt_systs['CMS_ff_total_syst_alt_func_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_syst_alt_func_' % vars()
    lt_systs['CMS_ff_total_qcd_syst_bkg_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_qcd_syst_bkg_' % vars()
    lt_systs['CMS_ff_total_wjets_syst_bkg_%(ch)s_%(yr)s' % vars()] = 'wt_ff_mssm_wjets_syst_bkg_' % vars()

    for njet in [0,1]:
      for jetpt in ['low','med','high']:
        for i in [1,2,3]:
          lt_systs[('CMS_ff_total_qcd_stat_njet%(njet)i_jet_pt_%(jetpt)s_unc%(i)i_%(ch)s_%(yr)s' % vars())] = 'wt_ff_mssm_qcd_stat_njet%(njet)i_jet_pt_%(jetpt)s_unc%(i)i_' % vars()
        for i in [1,2,3,4]:
          lt_systs[('CMS_ff_total_wjets_stat_njet%(njet)i_jet_pt_%(jetpt)s_unc%(i)i_%(ch)s_%(yr)s' % vars())] = 'wt_ff_mssm_wjets_stat_njet%(njet)i_jet_pt_%(jetpt)s_unc%(i)i_' % vars()
        for i in [1,2,3]:
          lt_systs[('CMS_ff_total_ttbar_stat_jet_pt_%(jetpt)s_unc%(i)i_%(ch)s_%(yr)s' % vars())] = 'wt_ff_mssm_ttbar_stat_jet_pt_%(jetpt)s_unc%(i)i_' % vars()

    if options.qcd_ff_closure or options.w_ff_closure:
      for key,val in lt_systs.items():
        if options.qcd_ff_closure:
          lt_systs[key] = '(wt_ff_mssm_qcd_1/wt_ff_mssm_1)*%(val)s' % vars()
        elif options.w_ff_closure:
          lt_systs[key] = '(wt_ff_mssm_wjets_1/wt_ff_mssm_1)*%(val)s' % vars()


    for template_name in lt_systs:
      weight_name = lt_systs[template_name]
      systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   weight_name+'up',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT','ZLL','TT','VV'], True)
      systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', weight_name+'down', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT','ZLL','TT','VV'], True)

  template_name = 'CMS_ff_total_sub_syst_%(ch)s_%(yr)s' % vars()
  systematics[template_name+'_up']   = ('' , '_'+template_name+'Up',   'wt_ff',   ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT','ZLL','TT','VV'], True)
  systematics[template_name+'_down'] = ('' , '_'+template_name+'Down', 'wt_ff', ['EWKZ','ZTT','ZJ','ZL','VVT','VVJ','TTT','TTJ','QCD','W','signal','EmbedZTT','ZLL','TT','VV'], True)


if options.syst_qcd_bkg: 
    systematics['qcd_sub_up'] = ('','_'+options.syst_qcd_bkg+'Up', 'wt', ['EWKZ','ZTT','ZJ','ZL','ZLL','VV','VVT','VVJ','TT','TTT','TTJ','W','signal','EmbedZTT','jetFakes'], False)
    systematics['qcd_sub_down'] = ('','_'+options.syst_qcd_bkg+'Down', 'wt', ['EWKZ','ZTT','ZJ','ZL','ZLL','VV','VVT','VVJ','TT','TTT','TTJ','W','signal','EmbedZTT','jetFakes'], False)


if options.no_systs: systematics = {'default' : ('','', 'wt', [], False)}

# sort systematics by tree's input directory name        
systematics = OrderedDict(sorted(systematics.items(), key=lambda key: key[1]))

if options.qcd_os_ss_ratio > 0:
    qcd_os_ss_ratio = options.qcd_os_ss_ratio
else:
    if options.analysis in ['sm','cpprod','cpdecay','mssmrun2','vlq']:
      if options.channel == 'et':
          qcd_os_ss_ratio = 1.0
          qcd_os_ss_ratio = 1.13 
          if options.cat == '0jet':
            if options.era in ['cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','cpsummer17','UL_17','cp18','UL_18','mvadm2016']: qcd_os_ss_ratio = 0.97
            else: qcd_os_ss_ratio = 1.0
          elif options.cat == 'boosted':
            if options.era in ['cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','cpsummer17','UL_17','cp18','UL_18','mvadm2016']: qcd_os_ss_ratio = 1.61
            else: qcd_os_ss_ratio = 1.28
          elif options.cat in ['vbf']: qcd_os_ss_ratio = 1.0
          elif options.cat in ['dijet','dijet_lowboost','dijet_boosted',
                  # 'ggh_highMjj','ggh_boosted_highMjj','qqh_boosted_highMjj','ggh_tight_highMjj','ggh_loose_highMjj','ggh_loose_boost_highMjj','ggh_tight_boost_highMjj','qqh_highMjj','misc_highMjj','qcd_highMjj','qqh_highMjj','tt_highMjj','zll_highMjj','ztt_highMjj','fake_highMjj','jetFakes_highMjj','zttEmbed_highMjj','jetFakes_tight_highMjj','jetFakes_loose_highMjj'
                  ]: qcd_os_ss_ratio = 1.82
          elif options.cat in ['btag']: qcd_os_ss_ratio = 1.16

      elif options.channel in ['mt','mj']: 
          #qcd_os_ss_ratio = 1.07
          qcd_os_ss_ratio = 1.12
          if options.cat == '0jet':
            if options.era in ['cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','cpsummer17','UL_17','cp18','UL_18','mvadm2016']: qcd_os_ss_ratio = 1.15
            else: qcd_os_ss_ratio = 1.07
          elif options.cat == 'boosted':
            if options.era in ['cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','cpsummer17','UL_17','cp18','UL_18','mvadm2016']: qcd_os_ss_ratio = 1.18
            else: qcd_os_ss_ratio = 1.06
          elif options.cat in ['vbf']: qcd_os_ss_ratio = 1.0
          elif options.cat in ['dijet','dijet_lowboost','dijet_boosted',
                  # 'ggh_highMjj','ggh_boosted_highMjj','qqh_boosted_highMjj','ggh_tight_highMjj','ggh_loose_highMjj','ggh_loose_boost_highMjj','ggh_tight_boost_highMjj','qqh_highMjj','misc_highMjj','qcd_highMjj','qqh_highMjj','tt_highMjj','zll_highMjj','ztt_highMjj','fake_highMjj','jetFakes_highMjj','zttEmbed_highMjj','jetFakes_tight_highMjj','jetFakes_loose_highMjj'
                  ]: qcd_os_ss_ratio = 1.23
          elif options.cat in ['btag']: qcd_os_ss_ratio = 1.01
      elif options.channel == 'zmm' or options.channel == 'zee':
          qcd_os_ss_ratio = 1.07   
      elif options.channel == 'em': 
          qcd_os_ss_ratio = 1.0
      else:
          qcd_os_ss_ratio = 1.0  
    else:
      if options.channel == 'et':
          qcd_os_ss_ratio = 1.02
          if options.cat == 'inclusive': qcd_os_ss_ratio = 1.13
          elif options.cat in ['nobtag', 'nobtag_tightmt', 'nobtag_loosemt']: qcd_os_ss_ratio = 1.11
          elif options.cat in ['btag', 'btag_tightmt', 'btag_loosemt']: qcd_os_ss_ratio = 1.16
      elif options.channel in ['mt','mj']: 
          qcd_os_ss_ratio = 1.18
          if options.cat == 'inclusive': qcd_os_ss_ratio = 1.12
          elif options.cat in ['nobtag', 'nobtag_tightmt', 'nobtag_loosemt']: qcd_os_ss_ratio = 1.14
          elif options.cat in ['btag', 'btag_tightmt', 'btag_loosemt']: qcd_os_ss_ratio = 1.01
      elif options.channel == 'zmm' or options.channel == 'zee':
          qcd_os_ss_ratio = 1.06   
      else:
          qcd_os_ss_ratio = 1.0
#if options.do_ss:
#    qcd_os_ss_ratio = 1.0


# Get array of signal masses to process        
ggh_masses=None
bbh_masses=None
sm_masses=None
bbh_masses_powheg=None
ggh_masses_powheg=None
if options.sm_masses != "": sm_masses = options.sm_masses.split(',')
if options.ggh_masses != "": ggh_masses = options.ggh_masses.split(',')
if options.bbh_masses != "": bbh_masses = options.bbh_masses.split(',')
if options.bbh_nlo_masses != "": bbh_nlo_masses = options.bbh_nlo_masses.split(',')
if options.ggh_masses_powheg != "": ggh_masses = options.ggh_masses_powheg.split(',')
if options.bbh_masses_powheg != "": bbh_masses = options.bbh_masses_powheg.split(',')


ROOT.TH1.SetDefaultSumw2(True)

# All functions defined here

def BuildCutString(wt='', sel='', cat='', sign='os',bkg_sel=''):
    full_selection = '(1)'
    if wt != '':
        full_selection = '('+wt+')'
    if sel != '':
        full_selection += '*('+sel+')'
    if sign != '':
        full_selection += '*('+sign+')'
    if bkg_sel != '':
        full_selection += '*('+bkg_sel+')'
    if cat != '':
        full_selection += '*('+cat+')'
    return full_selection

def GetZTTNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['ztt_sel'])
    return ana.SummedFactory('ZTT'+add_name, samples, plot, full_selection)
def GetEmbeddedNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    wt_ = wt
    if options.analysis in ['mssmrun2','vlq']:
      if options.channel == 'tt': 
        wt_+='*(pt_1/gen_match_1_pt<1.5&&pt_2/gen_match_2_pt<1.5)'
        if options.year == "2016": wt_+='*1.008'
        else: wt_+='*1.01'
      if options.channel in ['et','mt']: wt_+='*(pt_2/gen_match_2_pt<1.5)*1.005'
      #wt_+='*wt_emb_sel_kit/trackingweight_1'
      wt_+='*1/trackingweight_1'
    if options.channel == 'em':
      #for em channel there are non-closures wrt data and MC which are corrected here with these additional correction factors
      if options.era in ['cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','mvadm2016']: wt_+='*1.106'
      if options.era in ['cpsummer17','UL_17']: wt_+='*1.080'
      if options.era in ['cp18','UL_18']: wt_+='*1.101'
    full_selection = BuildCutString(wt_+'*(wt<2)', sel, cat, OSSS, z_sels['ztt_sel'])
    return ana.SummedFactory('EmbedZTT'+add_name, samples, plot, full_selection)


def GetZLLNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zll_sel'])
    return ana.SummedFactory('ZLL'+add_name, samples, plot, full_selection)

def GetZLNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zl_sel'])
    return ana.SummedFactory('ZL'+add_name, samples, plot, full_selection)

def GetZLEmbeddedNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    wt_=wt
    if options.analysis in ['mssmrun2','vlq']:
      #wt_+='*wt_emb_sel_kit/(trackingweight_1*trackingweight_2)'
      #wt_+='*wt_zpt_embed_ic/(trackingweight_1*trackingweight_2)'
      wt_+='*1/(trackingweight_1*trackingweight_2)'
    full_selection = BuildCutString(wt_, sel, cat, OSSS, '1')
    return ana.SummedFactory('EmbedZL'+add_name, samples, plot, full_selection)

def GetZJNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, z_sels['zj_sel'])
    return ana.SummedFactory('ZJ'+add_name, samples, plot, full_selection)

def GenerateZLL(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True, doZL=True, doZJ=True):
    if options.channel == 'em':
        zll_node = GetZLLNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
        ana.nodes[nodename].AddNode(zll_node)
    else:
        if doZL:
            zl_node = GetZLNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
            ana.nodes[nodename].AddNode(zl_node)
        if doZJ:
            zj_node = GetZJNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
            ana.nodes[nodename].AddNode(zj_node)  

def GenerateZTT(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    ztt_node = GetZTTNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)  
    ana.nodes[nodename].AddNode(ztt_node)
    split_taus = False
    if split_taus:
      ztt_node_rho = GetZTTNode(ana, '_rho'+add_name, samples, plot, wt, sel+'&&tauFlag_2==1', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(ztt_node_rho) 
      ztt_node_a1 = GetZTTNode(ana, '_a1'+add_name, samples, plot, wt, sel+'&&tauFlag_2==2', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(ztt_node_a1)
      ztt_node_pi = GetZTTNode(ana, '_pi'+add_name, samples, plot, wt, sel+'&&tauFlag_2==0', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(ztt_node_pi)
      ztt_node_other = GetZTTNode(ana, '_other'+add_name, samples, plot, wt, sel+'&&(tauFlag_2<0 || tauFlag_2>2)', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(ztt_node_other)  
    
def GenerateEmbedded(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    embed_node = GetEmbeddedNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
    ana.nodes[nodename].AddNode(embed_node)   
    if split_taus:
      embed_node_rho = GetEmbeddedNode(ana, '_rho'+add_name, samples, plot, wt, sel+'&&tauFlag_2==1', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(embed_node_rho)
      embed_node_a1 = GetEmbeddedNode(ana, '_a1'+add_name, samples, plot, wt, sel+'&&tauFlag_2==2', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(embed_node_a1)
      embed_node_pi = GetEmbeddedNode(ana, '_pi'+add_name, samples, plot, wt, sel+'&&tauFlag_2==0', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(embed_node_pi)
      embed_node_other = GetEmbeddedNode(ana, '_other'+add_name, samples, plot, wt, sel+'&&(tauFlag_2<0 || tauFlag_2>2)', cat, z_sels, get_os)
      ana.nodes[nodename].AddNode(embed_node_other) 
def GenerateZLEmbedded(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    embed_node = GetZLEmbeddedNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os)
    ana.nodes[nodename].AddNode(embed_node)  
    
def GetEWKZNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True, doJFakes=True):
    extra_sel = '1'
    if not doJFakes:
      if options.channel in ['et','mt']: extra_sel = 'gen_match_2<6'
      if options.channel == 'tt': extra_sel = '!(gen_match_1==6 || gen_match_2==6)'
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    if options.embedding: full_selection = BuildCutString(wt, sel, cat, OSSS, extra_top_sel+'&&'+extra_sel) 
    else: full_selection = BuildCutString(wt, sel, cat, OSSS, extra_sel)
    return ana.SummedFactory('EWKZ'+add_name, samples, plot, full_selection)

def GenerateEWKZ(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', z_sels={}, get_os=True):
    ewkz_node = GetEWKZNode(ana, add_name, samples, plot, wt, sel, cat, z_sels, get_os, options.method != 17)  
    ana.nodes[nodename].AddNode(ewkz_node) 

def GetggHWWNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, '1')
    return ana.SummedFactory('ggH_hww125'+add_name, samples, plot, full_selection)

def GetqqHWWNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', get_os=True):
    if get_os: OSSS = 'os'
    else: OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS, '1')
    return ana.SummedFactory('qqH_hww125'+add_name, samples, plot, full_selection)

def GenerateHWW(ana, add_name='', ggh_samples=[], qqh_samples=[], plot='', wt='', sel='', cat='', get_os=True, doggH=True, doqqH=True):
  if doggH:
      gghww_node = GetggHWWNode(ana, add_name, ggh_samples, plot, wt, sel, cat, get_os)
      ana.nodes[nodename].AddNode(gghww_node)
  if doqqH:
      qqhww_node = GetqqHWWNode(ana, add_name, qqh_samples, plot, wt, sel, cat, get_os)
      ana.nodes[nodename].AddNode(qqhww_node)

def GetTTTNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True):
  if get_os: OSSS = 'os'
  else: OSSS = '!os'  
  full_selection = BuildCutString(wt, sel, cat, OSSS, top_sels['ttt_sel'])
  return ana.SummedFactory('TTT'+add_name, samples, plot, full_selection)
  
def GetTTTforZTTNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True):
  if get_os: OSSS = 'os'
  else: OSSS = '!os'
  full_selection = BuildCutString(wt, sel, cat, OSSS, top_sels['ttt_ztt_sel'])
  return ana.SummedFactory('TTT_for_ZTT'+add_name, samples, plot, full_selection)


def GetTTJNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True):
  if get_os: OSSS = 'os'
  else: OSSS = '!os'  
  full_selection = BuildCutString(wt, sel, cat, OSSS, top_sels['ttj_sel'])
  return ana.SummedFactory('TTJ'+add_name, samples, plot, full_selection)

def GenerateTop(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', top_sels={}, get_os=True, doTTT=True, doTTJ=True):
  wt_=wt#+"*wt_tquark_down"
  if doTTT:
      ttt_node = GetTTTNode(ana, add_name, samples, plot, wt_, sel, cat, top_sels, get_os)
      ana.nodes[nodename].AddNode(ttt_node)
      if options.embedding:
        ttt_for_ztt_node = GetTTTforZTTNode(ana, add_name, samples, plot, wt_, sel, cat, top_sels, get_os)
        ana.nodes[nodename].AddNode(ttt_for_ztt_node)

      if split_taus:
       ttt_node_rho = GetTTTNode(ana, '_rho'+add_name, samples, plot, wt_, sel+'&&tauFlag_2==1', cat, top_sels, get_os)
       ana.nodes[nodename].AddNode(ttt_node_rho)
       ttt_node_a1 = GetTTTNode(ana, '_a1'+add_name, samples, plot, wt_, sel+'&&tauFlag_2==2', cat, top_sels, get_os)
       ana.nodes[nodename].AddNode(ttt_node_a1)
       ttt_node_pi = GetTTTNode(ana, '_pi'+add_name, samples, plot, wt_, sel+'&&tauFlag_2==0', cat, top_sels, get_os)
       ana.nodes[nodename].AddNode(ttt_node_pi)
       ttt_node_other = GetTTTNode(ana, '_other'+add_name, samples, plot, wt_, sel+'&&(tauFlag_2<0 || tauFlag_2>2)', cat, top_sels, get_os)
       ana.nodes[nodename].AddNode(ttt_node_other)

  if doTTJ:
      ttj_node = GetTTJNode(ana, add_name, samples, plot, wt_, sel, cat, top_sels, get_os)
      ana.nodes[nodename].AddNode(ttj_node)

def GetVVTNode(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True): 
  if get_os: OSSS = 'os'
  else: OSSS = '!os'  
  full_selection = BuildCutString(wt, sel, cat, OSSS, vv_sels['vvt_sel'])
  return ana.SummedFactory('VVT'+add_name, samples, plot, full_selection)

def GetVVTforZTTNode(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True):
  if get_os: OSSS = 'os'
  else: OSSS = '!os'
  full_selection = BuildCutString(wt, sel, cat, OSSS, vv_sels['vvt_ztt_sel'])
  return ana.SummedFactory('VVT_for_ZTT'+add_name, samples, plot, full_selection)


def GetVVJNode(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True): 
  if get_os: OSSS = 'os'
  else: OSSS = '!os'  
  full_selection = BuildCutString(wt, sel, cat, OSSS, vv_sels['vvj_sel'])
  return ana.SummedFactory('VVJ'+add_name, samples, plot, full_selection)

def GenerateVV(ana, add_name ='', samples=[], plot='', wt='', sel='', cat='', vv_sels={}, get_os=True, doVVT=True, doVVJ=True): 
  if doVVT:
      vvt_node = GetVVTNode(ana, add_name, samples, plot, wt, sel, cat, vv_sels, get_os)
      ana.nodes[nodename].AddNode(vvt_node)
      if options.embedding:
        vvt_for_ztt_node = GetVVTforZTTNode(ana, add_name, samples, plot, wt, sel, cat, vv_sels, get_os)
        ana.nodes[nodename].AddNode(vvt_for_ztt_node)
      if split_taus:
       vvt_node_rho = GetVVTNode(ana, '_rho'+add_name, samples, plot, wt, sel+'&&tauFlag_2==1', cat, vv_sels, get_os)
       ana.nodes[nodename].AddNode(vvt_node_rho)
       vvt_node_a1 = GetVVTNode(ana, '_a1'+add_name, samples, plot, wt, sel+'&&tauFlag_2==2', cat, vv_sels, get_os)
       ana.nodes[nodename].AddNode(vvt_node_a1)
       vvt_node_pi = GetVVTNode(ana, '_pi'+add_name, samples, plot, wt, sel+'&&tauFlag_2==0', cat, vv_sels, get_os)
       ana.nodes[nodename].AddNode(vvt_node_pi)
       vvt_node_other = GetVVTNode(ana, '_other'+add_name, samples, plot, wt, sel+'&&(tauFlag_2<0 || tauFlag_2>2)', cat, vv_sels, get_os)
       ana.nodes[nodename].AddNode(vvt_node_other)

  if doVVJ:
      vvj_node = GetVVJNode(ana, add_name, samples, plot, wt, sel, cat, vv_sels, get_os)
      ana.nodes[nodename].AddNode(vvj_node)
      
def GetWGNode(ana, add_name='', samples=[], plot='', wt='', sel='', cat='', get_os=True):
  if get_os:
      OSSS = 'os'
  else:
      OSSS = '!os'   
  full_selection = BuildCutString(wt, sel, cat, OSSS)
  wg_node = ana.SummedFactory('WGam'+add_name, samples, plot, full_selection)
  return wg_node

def GetWNode(ana, name='W', samples=[], data=[], plot='',plot_unmodified='', wt='', sel='', cat='', cat_data='', method=8, qcd_factor=qcd_os_ss_ratio, get_os=True):
  if get_os: OSSS = 'os'
  else: OSSS = '!os'
  full_selection = BuildCutString(wt, sel, cat, OSSS, '')
  if cats['w_shape'] != '': shape_cat = cats['w_shape']
  else: shape_cat = cat
  if method == 14:
      shape_cat = '(n_jets<=1 && n_loose_bjets>=1)*('+cats['baseline']+')'
  shape_selection = BuildCutString(wt, sel, shape_cat, OSSS, '')
  
  if method in [0, 8, 9, 15, 19]:
      w_node = ana.SummedFactory(name, samples, plot, full_selection)
  elif method in [10, 11]:
      control_sel = cats['w_sdb']+' && '+ OSSS
      w_control_full_selection = BuildCutString(wt, control_sel, cat, OSSS)
      w_control_full_selection_data = BuildCutString(wt, control_sel, cat_data, OSSS)
      subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt,control_sel,cat,cat_data,method,qcd_os_ss_ratio,True,False) 
      if shape_selection == full_selection:
          w_shape = None
      else:    
          w_shape = ana.SummedFactory('w_shape', samples, plot, shape_selection)
      w_node = HttWNode(name,
        ana.SummedFactory('data_obs', data, plot_unmodified, w_control_full_selection_data),
        subtract_node,
        ana.SummedFactory('W_cr', samples, plot, w_control_full_selection),
        ana.SummedFactory('W_sr', samples, plot, full_selection),
        w_shape)
  elif method in [12, 13, 14, 16, 23]:
      if method == 16:
          cat_nobtag = '('+cats['btag_wnobtag']+')*('+cats['baseline']+')'
          cat_nobtag_data = '('+cats_unmodified['btag_wnobtag']+')*('+cats_unmodified['baseline']+')'
          full_selection = BuildCutString(wt, sel, cat_nobtag, OSSS)
          ss_selection = BuildCutString(wt, '', cat_nobtag, '!os', '')
          os_selection = BuildCutString(wt, '', cat_nobtag, 'os', '')
          control_sel = cats['w_sdb']
          w_control_full_selection = BuildCutString(wt, control_sel, cat_nobtag, OSSS)
          w_control_full_selection_os = BuildCutString(wt, control_sel, cat_nobtag)
          w_control_full_selection_ss = BuildCutString(wt, control_sel, cat_nobtag, '!os')
          w_control_full_selection_os_data = BuildCutString(wt, control_sel, cat_nobtag_data)
          w_control_full_selection_ss_data = BuildCutString(wt, control_sel, cat_nobtag_data, '!os')
          btag_extrap_sel_num = BuildCutString(wt, sel, cat, OSSS, '')
          btag_extrap_sel_den = BuildCutString(wt, sel, cat_nobtag, OSSS, '')
          btag_extrap_num_node = ana.SummedFactory('btag', samples, plot, btag_extrap_sel_num)
          btag_extrap_den_node = ana.SummedFactory('no_btag', samples, plot, btag_extrap_sel_den)
          subtract_node_os = GetSubtractNode(ana,'_os',plot,plot_unmodified,wt,control_sel,cat_nobtag,cat_nobtag_data,method,qcd_os_ss_ratio,True,False) 
          subtract_node_ss = GetSubtractNode(ana,'_ss',plot,plot_unmodified,wt,control_sel,cat_nobtag,cat_nobtag_data,method,qcd_os_ss_ratio,False,False)
      elif method == 23:
          cat_nopt = '('+cats['dijet_lowboost']+')*('+cats['baseline']+')'
          cat_nopt_data = '('+cats_unmodified['dijet_lowboost']+')*('+cats_unmodified['baseline']+')'
          full_selection = BuildCutString(wt, sel, cat_nopt, OSSS)
          ss_selection = BuildCutString(wt, '', cat_nopt, '!os', '')
          os_selection = BuildCutString(wt, '', cat_nopt, 'os', '')
          control_sel = cats['w_sdb']
          w_control_full_selection = BuildCutString(wt, control_sel, cat_nopt, OSSS)
          w_control_full_selection_os = BuildCutString(wt, control_sel, cat_nopt)
          w_control_full_selection_ss = BuildCutString(wt, control_sel, cat_nopt, '!os')
          w_control_full_selection_os_data = BuildCutString(wt, control_sel, cat_nopt_data)
          w_control_full_selection_ss_data = BuildCutString(wt, control_sel, cat_nopt_data, '!os')
          btag_extrap_sel_num = BuildCutString(wt, sel, cat, OSSS, '')
          btag_extrap_sel_den = BuildCutString(wt, sel, cat_nopt, OSSS, '')
          btag_extrap_num_node = ana.SummedFactory('btag', samples, plot, btag_extrap_sel_num)
          btag_extrap_den_node = ana.SummedFactory('no_btag', samples, plot, btag_extrap_sel_den)
          subtract_node_os = GetSubtractNode(ana,'_os',plot,plot_unmodified,wt,control_sel,cat_nopt,cat_nopt_data,method,qcd_os_ss_ratio,True,False) 
          subtract_node_ss = GetSubtractNode(ana,'_ss',plot,plot_unmodified,wt,control_sel,cat_nopt,cat_nopt_data,method,qcd_os_ss_ratio,False,False)    
      else:
          full_selection = BuildCutString(wt, sel, cat, OSSS)
          ss_selection = BuildCutString(wt, '', cat, '!os', '')
          os_selection = BuildCutString(wt, '', cat, 'os', '')
          control_sel = cats['w_sdb']
          w_control_full_selection = BuildCutString(wt, control_sel, cat, OSSS)
          w_control_full_selection_os = BuildCutString(wt, control_sel, cat)
          w_control_full_selection_ss = BuildCutString(wt, control_sel, cat, '!os')
          w_control_full_selection_os_data = BuildCutString(wt, control_sel, cat_data)
          w_control_full_selection_ss_data = BuildCutString(wt, control_sel, cat_data, '!os')
          btag_extrap_num_node = None
          btag_extrap_den_node = None
          subtract_node_os = GetSubtractNode(ana,'_os',plot,plot_unmodified,wt,control_sel,cat,cat_data,method,qcd_os_ss_ratio,True,False) 
          subtract_node_ss = GetSubtractNode(ana,'_ss',plot,plot_unmodified,wt,control_sel,cat,cat_data,method,qcd_os_ss_ratio,False,False)
          
      if shape_selection == full_selection:
          w_shape = None
      else:    
          w_shape = ana.SummedFactory('w_shape', samples, plot, shape_selection)
      w_node = HttWOSSSNode(name,
        ana.SummedFactory('data_os', data, plot_unmodified, w_control_full_selection_os_data),
        subtract_node_os,
        ana.SummedFactory('data_ss', data, plot_unmodified, w_control_full_selection_ss_data),
        subtract_node_ss,
        ana.SummedFactory('W_cr', samples, plot, w_control_full_selection),
        ana.SummedFactory('W_sr', samples, plot, full_selection),
        ana.SummedFactory('W_os', samples, plot, os_selection),
        ana.SummedFactory('W_ss', samples, plot, ss_selection),
        w_shape,
        qcd_factor,
        get_os,
        btag_extrap_num_node,
        btag_extrap_den_node) 
      
  elif method in [21,22,24,25]:
    control_sel = cats['w_sdb']+' && '+ OSSS
    w_control_full_selection = BuildCutString(wt, control_sel, cat, OSSS)
    w_control_full_selection_data = BuildCutString(wt, control_sel, cat_data, OSSS)
    data_node=ana.SummedFactory('data_obs', data, plot_unmodified, w_control_full_selection_data)
    subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt,control_sel,cat,cat_data,method,qcd_os_ss_ratio,True,False) 
    
    qcd_subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt,cats['w_sdb'],cat,cat_data,8,qcd_os_ss_ratio,False,True)
    qcd_control_full_selection = BuildCutString(wt, cats['w_sdb'], cat, '!os')
    qcd_control_full_selection_data = BuildCutString(wt, cats['w_sdb'], cat_data, '!os')
    qcd_node = HttQCDNode('QCD'+add_name,
      ana.SummedFactory('data_ss', data, plot_unmodified, qcd_control_full_selection_data),
      qcd_subtract_node,
      qcd_factor,
      None)
    
    subtract_node.AddNode(qcd_node)
    
    if shape_selection == full_selection: w_shape = None
    else: w_shape = ana.SummedFactory('w_shape', samples, plot, shape_selection)
    
    if method in [22]:
      data_node = None
      wsf_num = GetWNode(ana, name, samples, data, plot,plot_unmodified, wt, sel, cat, cat_data, 21, qcd_factor, True)
      wsf_denum = GetWNode(ana, name, samples, data, plot,plot_unmodified, wt, sel, cat, cat_data, 8, qcd_factor, True)
    if method in [24]:
      data_node = None
      cat_nopt = '('+cats['dijet']+')*('+cats['baseline']+')'
      cat_nopt_data = '('+cats_unmodified['dijet_lowboost']+')*('+cats_unmodified['baseline']+')'
      wsf_num = GetWNode(ana, name, samples, data, plot,plot_unmodified, wt, sel, cat, cat_data, 8, qcd_factor, True)
      wsf_denum = GetWNode(ana, name, samples, data, plot,plot_unmodified, wt, sel, cat_nopt, cat_nopt_data, 8, qcd_factor, True)  
    if method in [25]:
      data_node = None
      wsf_num = GetWNode(ana, name, samples, data, plot,plot_unmodified, wt, sel, cat, cat_data, 24, qcd_factor, True)
      wsf_denum = GetWNode(ana, name, samples, data, plot,plot_unmodified, wt, sel, cat, cat_data, 8, qcd_factor, True)   
    else:
      wsf_num = None
      wsf_denum = None    
    
    w_node = HttWNode(name,
      data_node,
      subtract_node,
      ana.SummedFactory('W_cr', samples, plot, w_control_full_selection),
      ana.SummedFactory('W_sr', samples, plot, full_selection),
      w_shape,
      wsf_num,
      wsf_denum)
  return w_node

def GenerateW(ana, add_name='', samples=[], data=[], wg_samples=[], plot='', plot_unmodified='', wt='', sel='', cat='', cat_data='', method=8, qcd_factor=qcd_os_ss_ratio, get_os=True):
  w_node_name = 'W'  
  if options.channel == 'em':
      w_total_node = SummedNode('W'+add_name)
      w_total_node.AddNode(GetWGNode(ana, add_name, wg_samples, plot, wt, sel, cat, get_os))
      ana.nodes[nodename].AddNode(GetWGNode(ana, add_name, wg_samples, plot, wt, sel, cat, get_os))
      w_node_name+='J'
  ana.nodes[nodename].AddNode(GetWNode(ana, w_node_name+add_name, samples, data, plot, plot_unmodified, wt, sel, cat, cat_data, method, qcd_factor, get_os))
  if options.channel == 'em':
      w_total_node.AddNode(GetWNode(ana, w_node_name+add_name, samples, data, plot, plot_unmodified, wt, sel, cat, cat_data, method, qcd_factor, get_os))
      ana.nodes[nodename].AddNode(w_total_node)

def GetSubtractNode(ana,add_name,plot,plot_unmodified,wt,sel,cat,cat_data,method,qcd_os_ss_ratio,OSSS,includeW=False,w_shift=None):
  subtract_node = SummedNode('total_bkg'+add_name)
  if includeW:
      if w_shift is not None: w_wt = '%s*%f' %(wt,w_shift)
      else: w_wt = wt
      w_node = GetWNode(ana, 'W', wjets_samples, data_samples, plot, plot_unmodified, w_wt, sel, cat, cat_data,method, qcd_os_ss_ratio, OSSS)
      subtract_node.AddNode(w_node)
  ttt_node = GetTTTNode(ana, "", top_samples, plot, wt, sel, cat, top_sels, OSSS)
  ttj_node = GetTTJNode(ana, "", top_samples, plot, wt, sel, cat, top_sels, OSSS)
  vvt_node = GetVVTNode(ana, "", vv_samples, plot, wt, sel, cat, vv_sels, OSSS)
  vvj_node = GetVVJNode(ana, "", vv_samples, plot, wt, sel, cat, vv_sels, OSSS)
  subtract_node.AddNode(ttt_node)
  subtract_node.AddNode(ttj_node)
  subtract_node.AddNode(vvt_node)
  subtract_node.AddNode(vvj_node)
  if options.embedding and options.channel != 'zmm': 
    embed_node = GetEmbeddedNode(ana, "", embed_samples, plot, wt, sel, cat, z_sels, OSSS)
    subtract_node.AddNode(embed_node)
  else:
    ztt_node = GetZTTNode(ana, "", ztt_samples, plot, wt, sel, cat, z_sels, OSSS)
    subtract_node.AddNode(ztt_node)
  if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','tauid2016','tauid2017','cpsummer17','UL_17','cp18','UL_18','mvadm2016'] and options.method!=0: 
    ewkz_node = GetEWKZNode(ana, "", ewkz_samples, plot, wt, sel, cat, z_sels, OSSS)
    subtract_node.AddNode(ewkz_node)
  if options.channel not in ["em"]:
      zl_node = GetZLNode(ana, "", ztt_samples, plot, wt, sel, cat, z_sels, OSSS)
      zj_node = GetZJNode(ana, "", ztt_samples, plot, wt, sel, cat, z_sels, OSSS)
      subtract_node.AddNode(zl_node)
      subtract_node.AddNode(zj_node)
  if options.channel in ["em"]:
      zll_node = GetZLLNode(ana, "", ztt_samples, plot, wt, sel, cat, z_sels, OSSS)
      subtract_node.AddNode(zll_node)
  if options.channel == "em":
      wg_node = GetWGNode(ana, "", wgam_samples, plot, wt, sel, cat, OSSS)
      subtract_node.AddNode(wg_node)
  #print(subtract_node.shape)
  return subtract_node
      
def GenerateQCD(ana, add_name='', data=[], plot='', plot_unmodified='', wt='', sel='', cat='', cat_data='', method=8, qcd_factor=qcd_os_ss_ratio, get_os=True,w_shift=None):
    shape_node = None
    OSSS = "!os"
    if get_os: OSSS = "os"
    if method == 0:
        subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt,sel,cat,cat_data,8,1.0,get_os,True)
        full_selection = BuildCutString(wt, sel, cat_data, OSSS)
        ana.nodes[nodename].AddNode(HttQCDNode('QCD'+add_name,
          ana.SummedFactory('data', data, plot_unmodified, full_selection),
          subtract_node,
          1.0))
        return

    if options.channel != 'tt':
        sub_shift='*1.0'
        if 'qcd_sub_up' in systematic: sub_shift = '*1.1'
        if 'qcd_sub_down' in systematic: sub_shift = '*0.9'
        if method in [9, 11, 13, 14]:
            if method in [9, 11, 13]: 
              shape_cat = '('+cats[options.cat]+')*('+cats['qcd_loose_shape']+')'
              shape_cat_data = '('+cats_unmodified[options.cat]+')*('+cats_unmodified['qcd_loose_shape']+')'
            elif method == 14:
                shape_cat = '(n_jets<=1 && n_loose_bjets>=1)*('+cats['baseline']+')'
                shape_cat_data = '(n_jets<=1 && n_loose_bjets>=1)*('+cats_unmodified['baseline']+')'
            shape_selection = BuildCutString(wt, sel, shape_cat_data, '!os')
            subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt,sel,shape_cat,shape_cat_data,method,qcd_os_ss_ratio,False,True)  
            shape_node = SubtractNode('shape', ana.SummedFactory('data_ss', data, plot_unmodified, shape_selection), subtract_node)
        
        #if options.channel == 'em': qcd_os_ss_factor = 1
        qcd_os_ss_factor = qcd_factor
        weight = wt
        if method in [15,19]:
            #qcd_os_ss_factor = 1
            if get_os and options.channel == "em":
                weight = wt+'*wt_em_qcd' 
                #if options.era in ['cpsummer16','cpdecay16',"legacy16",'mvadm2016']: weight+='*wt_em_qcd_extrapup'
            if options.em_qcd_weight != '':
                weight=wt+'*'+options.em_qcd_weight
            if method == 19:
                shape_selection = BuildCutString(weight, sel, cats_unmodified['em_shape_cat'], '!os')
                subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,weight+sub_shift,sel,cats['em_shape_cat'],cats_unmodified['em_shape_cat'],method,1,False,True)
                shape_node = SubtractNode('shape', ana.SummedFactory('data_ss',data, plot_unmodified, shape_selection), subtract_node)
             
        if cats['qcd_shape'] != "" or (w_shift is not None and w_shift!=1.0):
            add_shape = False
            if cats['qcd_shape'] == '': shape_cat = cat
            else: 
              shape_cat = cats['qcd_shape']
              add_shape = True
            if cats['qcd_shape'] == '': shape_cat_data = cat_data
            else: 
              shape_cat_data = cats_unmodified['qcd_shape']
              add_shape = True
            if add_shape:
              shape_selection = BuildCutString(weight, sel, shape_cat_data, '!os')
              if method == 21: subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,weight,sel,shape_cat,shape_cat_data,22,1,False,True,w_shift)
              else: subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,weight,sel,shape_cat,shape_cat_data,method,1,False,True,w_shift)
              shape_node = SubtractNode('shape', ana.SummedFactory('data_ss',data, plot_unmodified, shape_selection), subtract_node)
            else: shape_node = None
        
        full_selection = BuildCutString(weight, sel, cat_data, '!os')
        if method == 21: 
          subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,weight,sel,cat,cat_data,22,qcd_os_ss_ratio,False,True)
        elif method == 24: 
          subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,weight,sel,cat,cat_data,25,qcd_os_ss_ratio,False,True)  
        else: subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,weight+sub_shift,sel,cat,cat_data,method,qcd_os_ss_ratio,False,True)
        if get_os: qcd_ratio = qcd_os_ss_factor
        else: qcd_ratio = 1.0
        ana.nodes[nodename].AddNode(HttQCDNode('QCD'+add_name,
          ana.SummedFactory('data_ss', data, plot_unmodified, full_selection),
          subtract_node,
          qcd_ratio,
          shape_node))
        
    else:
        if method == 8:
          qcd_sdb_cat = cats[options.cat]+' && '+cats['tt_qcd_norm'] 
          qcd_sdb_cat_data = cats_unmodified[options.cat]+' && '+cats_unmodified['tt_qcd_norm']
          
          subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt,sel,cat,cat_data,method,qcd_os_ss_ratio,False,True)        
          num_selection = BuildCutString(wt, sel, cat_data, '!os')
          num_node = SubtractNode('ratio_num',
                       ana.SummedFactory('data', data, plot_unmodified, num_selection),
                       subtract_node)
          if options.analysis == 'mssmsummer16': tau_id_wt = 'wt_tau2_id_loose'
          else: tau_id_wt = '1'
          subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt+'*'+tau_id_wt,sel,qcd_sdb_cat,qcd_sdb_cat_data,method,qcd_os_ss_ratio,False,True)
          den_selection = BuildCutString(wt, sel, qcd_sdb_cat_data, '!os')
          den_node = SubtractNode('ratio_den',
                       ana.SummedFactory('data', data, plot_unmodified, den_selection),
                       subtract_node)
          shape_node = None   
          full_selection = BuildCutString(wt, sel, qcd_sdb_cat_data, OSSS)
          subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt+'*'+tau_id_wt,sel,qcd_sdb_cat,qcd_sdb_cat_data,method,qcd_os_ss_ratio,get_os,True)
  
          ana.nodes[nodename].AddNode(HttQCDNode('QCD'+add_name,
            ana.SummedFactory('data', data, plot_unmodified, full_selection),
            subtract_node,
            1,
            shape_node,
            num_node,
            den_node))
        if method == 9:
          #osss = '(3.68434 -1.8015*dR + 0.488955*dR*dR -0.0489159*dR*dR*dR)' # this line to be inclusive by decay mode
          osss = '((tau_decay_mode_1<2&&tau_decay_mode_2<2)*(3.555-1.689*dR+0.4741*dR*dR-0.04973*dR*dR*dR) + ((tau_decay_mode_1<2&&tau_decay_mode_2>9)||(tau_decay_mode_1>9&&tau_decay_mode_2<2))*(2.851-1.008*dR+0.2374*dR*dR-0.02331*dR*dR*dR) + (tau_decay_mode_1>9&&tau_decay_mode_2>9)*(3.391-2.21*dR+0.7623*dR*dR-0.08775*dR*dR*dR))'
          weight = wt
          if get_os:
              weight = wt+'*'+osss
          full_selection = BuildCutString(weight, sel, cat_data, '!os')
          subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,weight,sel,cat,cat_data,method,qcd_os_ss_ratio,False,True)
          ana.nodes[nodename].AddNode(HttQCDNode('QCD'+add_name,
          ana.SummedFactory('data_ss', data, plot_unmodified, full_selection),
          subtract_node,
          1.0,
          None))
            
        
def GenerateFakeTaus(ana, add_name='', data=[], plot='',plot_unmodified='', wt='', sel='', cat_name='',get_os=True,ff_syst_weight=None):

    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'

    sub_wt=''
    if 'sub_syst' in add_name and 'Up' in add_name: sub_wt='*1.1'
    if 'sub_syst' in add_name and 'Down' in add_name: sub_wt='*0.9'

    wp = options.wp
    # Select data from anti-isolated region
    if options.channel != "tt":
        if options.channel == 'mt':
            if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','cpsummer17','UL_17','mvadm2016','cp18','UL_18']:
              anti_isolated_sel = cats['baseline'].replace('deepTauVsJets_%(wp)s_2>0.5' % vars(),'deepTauVsJets_%(wp)s_2<0.5 && deepTauVsJets_vvvloose_2>0.5' % vars())
            else: 
              anti_isolated_sel = '(iso_1<0.15 && mva_olddm_tight_2<0.5 && mva_olddm_vloose_2>0.5 && antiele_2 && antimu_2 && !leptonveto && trg_singlemuon)'
        elif options.channel == 'et': 
            if options.era in ["smsummer16","cpsummer16","cpdecay16","legacy16",'UL_16_preVFP','UL_16_postVFP',"cpsummer17",'UL_17',"mvadm2016","cp18",'UL_18']:
              anti_isolated_sel = cats['baseline'].replace('deepTauVsJets_%(wp)s_2>0.5' % vars(),'deepTauVsJets_%(wp)s_2<0.5 && deepTauVsJets_vvvloose_2>0.5' % vars())
            else: 
              anti_isolated_sel = '(iso_1<0.1  && mva_olddm_tight_2<0.5 && mva_olddm_vloose_2>0.5 && antiele_2 && antimu_2 && !leptonveto && trg_singleelectron)'
        ff_cat = cats[cat_name] +" && "+ anti_isolated_sel
        ff_cat_data = cats_unmodified[cat_name] +" && "+ anti_isolated_sel
        if options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','cpsummer17','UL_17','mvadm2016','cp18','UL_18']:
          if ff_syst_weight is not None and 'sub_syst' not in add_name: fake_factor_wt_string = '('+ff_syst_weight+'_1)'
          else:
            if options.analysis in ['cpprod']: 
              fake_factor_wt_string = "wt_ff_us_1"
              fake_factor_wt_string = "wt_ff_dmbins_1"
            elif options.analysis in ['mssmrun2','vlq']:
              #if options.w_ff_closure or options.qcd_ff_closure:
                #json_name = 'scripts/ff_strings.json'
                #with open(json_name) as json_file:
                  #ff_dict = json.load(json_file)
              if options.w_ff_closure:
                fake_factor_wt_string = "wt_ff_mssm_wjets_1"
                if options.ml_ff: fake_factor_wt_string = "wt_ff_reweight_wjets_dr_1"
                #fake_factor_wt_string = RawFFFromString(ff_dict[options.channel][options.year]['wjets'])
                #fake_factor_wt_string = ff_dict[options.channel][options.year]['wjets']
              elif options.qcd_ff_closure:
                fake_factor_wt_string = "wt_ff_mssm_qcd_1"
                if options.ml_ff: fake_factor_wt_string = "wt_ff_reweight_qcd_raw_1"
                #fake_factor_wt_string = ff_dict[channel][year]['qcd']
              else:
                fake_factor_wt_string = "wt_ff_mssm_1"
                if options.ml_ff: fake_factor_wt_string = "((wt_ff_reweight_qcd_1*ff_frac_mssm_qcd) + (wt_ff_reweight_wjets_1*ff_frac_mssm_wjets) + (wt_ff_reweight_ttbar_1*ff_frac_mssm_ttbar))"

            else: fake_factor_wt_string = "wt_ff_1"
        else:
          if ff_syst_weight is not None: fake_factor_wt_string = ff_syst_weight
          else: fake_factor_wt_string = "wt_ff_"+options.cat
          fake_factor_wt_string+='*wt_tau_id_loose'
        if wt is not "": wt+="*"+fake_factor_wt_string
        else: wt=fake_factor_wt_string

        full_selection = BuildCutString(wt, sel, ff_cat_data, OSSS, '')
        # Calculate FF for anti-isolated data (f1) then subtract contributions from real taus (f2)
        f1 = ana.SummedFactory('data', data, plot_unmodified, full_selection)

        if (not options.w_ff_closure and not options.qcd_ff_closure) or options.ml_ff:
          f2 = GetSubtractNode(ana,'',plot,plot_unmodified,wt+sub_wt,sel+'&&(gen_match_2<6)',ff_cat,ff_cat_data,8,1.0,get_os,True)
        elif options.qcd_ff_closure and not options.ml_ff:
          f2 = GetSubtractNode(ana,'',plot,plot_unmodified,wt+sub_wt,sel,ff_cat,ff_cat_data,8,1,get_os,True)
        elif options.w_ff_closure and not options.ml_ff:
          f2 = GetSubtractNode(ana,'',plot,plot_unmodified,wt+sub_wt,sel,ff_cat,ff_cat_data,8,1,get_os,False)
          full_selection_ss = BuildCutString(wt, sel, ff_cat_data, '!os', '')
          qcd_node =  SubtractNode('qcd', ana.SummedFactory('data', data, plot_unmodified, full_selection_ss), GetSubtractNode(ana,'',plot,plot_unmodified,wt+sub_wt,sel,ff_cat,ff_cat_data,12,1.0,False,True))
          f2.AddNode(qcd_node)
        ana.nodes[nodename].AddNode(SubtractNode('jetFakes'+add_name, f1, f2))

    if options.channel == 'tt':
        anti_isolated_sel_1 = '(mva_olddm_medium_1<0.5 && mva_olddm_vloose_1>0.5 && mva_olddm_medium_2>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto)'
        anti_isolated_sel_2 = '(mva_olddm_medium_2<0.5 && mva_olddm_vloose_2>0.5 && mva_olddm_medium_1>0.5 && antiele_1 && antimu_1 && antiele_2 && antimu_2 && !leptonveto)'
        if options.era == "mssmsummer16": 
          anti_isolated_sel_1 +=" && trg_doubletau"
          anti_isolated_sel_2 +=" && trg_doubletau"
        if options.era in ["smsummer16","cpsummer16","cpdecay16","legacy16",'UL_16_preVFP','UL_16_postVFP',"mvadm2016"]:
          anti_isolated_sel_1 = cats['baseline'].replace('mva_olddm_tight_1>0.5','mva_olddm_tight_1<0.5 && mva_olddm_vloose_1>0.5')
          anti_isolated_sel_2 = cats['baseline'].replace('mva_olddm_tight_2>0.5','mva_olddm_tight_2<0.5 && mva_olddm_vloose_2>0.5')
        if options.era in ["cpsummer17",'UL_17',"cp18",'UL_18','cpdecay16','legacy16','UL_16_preVFP','UL_16_postVFP']: # need to do also for MVA cats for 2016!
          anti_isolated_sel_2 = cats['baseline'].replace('deepTauVsJets_%(wp)s_2>0.5' % vars(),'deepTauVsJets_%(wp)s_2<0.5 && deepTauVsJets_vvvloose_2>0.5' % vars())
          anti_isolated_sel_1 = cats['baseline'].replace('deepTauVsJets_%(wp)s_1>0.5' % vars(),'deepTauVsJets_%(wp)s_1<0.5 && deepTauVsJets_vvvloose_1>0.5' % vars())
        ff_cat_1 = cats[cat_name] +" && "+ anti_isolated_sel_1
        ff_cat_2 = cats[cat_name] +" && "+ anti_isolated_sel_2
        ff_cat_1_data = cats_unmodified[cat_name] +" && "+ anti_isolated_sel_1
        ff_cat_2_data = cats_unmodified[cat_name] +" && "+ anti_isolated_sel_2
        if ff_syst_weight is not None and 'sub_syst' not in add_name: 
            fake_factor_wt_string_1 = '('+ff_syst_weight+'_1)'
            fake_factor_wt_string_2 = ff_syst_weight+'_2'

#wt_ff_dmbins_qcd_syst_down
            if options.analysis in ['cpprod'] and 'ff_tt_qcd_syst' in add_name:
              flat = 1.04
              if (options.era=='cpsummer17' or options.era == 'UL_17'): flat=1.07
              if (options.era=='cp18'or options.era == 'UL_18'): flat=1.06

              if 'njets0' in add_name:
                if 'Up' in add_name:
                  ff_syst_weight = '((n_jets==0)*wt_ff_dmbins_qcd_syst_down_1*%(flat)s + (n_jets>0)*wt_ff_dmbins_1)' % vars()
                if 'Down' in add_name:
                  ff_syst_weight = '((n_jets==0)*(2*wt_ff_dmbins_1-wt_ff_dmbins_qcd_syst_down_1*%(flat)s) + (n_jets>0)*wt_ff_dmbins_1)' % vars()
              elif 'njets1' in add_name:
                if 'Up' in add_name:
                  ff_syst_weight = '((n_jets==0)*wt_ff_dmbins_1 + (n_jets>0)*wt_ff_dmbins_qcd_syst_down_1*%(flat)s)' % vars()
                if 'Down' in add_name:
                  ff_syst_weight = '((n_jets==0)*wt_ff_dmbins_1 + (n_jets>0)*(2*wt_ff_dmbins_1-wt_ff_dmbins_qcd_syst_down_1*%(flat)s))' % vars()

              fake_factor_wt_string_1 = ff_syst_weight


            if options.analysis in ['cpdecay'] and 'ff_tt_qcd_syst' in add_name:
              if(options.era=='cp18' or options.era == 'UL_18'):
                if 'njets0' in add_name:
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0)*(1.24845-0.0578994*dR)/1.06 + (n_jets>0))*wt_ff_1'
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0)/(1.24845-0.0578994*dR)*1.06 + (n_jets>0))*wt_ff_1' 
                elif 'njets1' in add_name:
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0) + (n_jets>0)*(1.15364-0.0317013*dR)/1.06)*wt_ff_1' 
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0) + (n_jets>0)/(1.15364-0.0317013*dR)*1.06)*wt_ff_1'
                else: 
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0)*(1.24845-0.0578994*dR)/1.06 + (n_jets>0)*(1.15364-0.0317013*dR)/1.06)*wt_ff_1'     
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0)/(1.24845-0.0578994*dR)*1.06 + (n_jets>0)/(1.15364-0.0317013*dR)*1.06)*wt_ff_1'

              if (options.era=='cpsummer17' or options.era == 'UL_17'):
                if 'njets0' in add_name:
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0)*(1.40415-0.101317*dR)/1.07 + (n_jets>0))*wt_ff_1'
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0)/(1.40415-0.101317*dR)*1.07 + (n_jets>0))*wt_ff_1'
                elif 'njets1' in add_name:
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0) + (n_jets>0)*(1.16253-0.031154*dR)/1.07)*wt_ff_1'
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0) + (n_jets>0)/(1.16253-0.031154*dR)*1.07)*wt_ff_1'
                else:
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0)*(1.40415-0.101317*dR)/1.07 + (n_jets>0)*(1.16253-0.031154*dR)/1.07)*wt_ff_1'
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0)/(1.40415-0.101317*dR)*1.07 + (n_jets>0)/(1.16253-0.031154*dR)*1.07)*wt_ff_1'

              if options.era in ['legacy16','UL_16_preVFP','UL_16_postVFP']:
                if 'njets0' in add_name:
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0)*(1.22618-0.0546679*dR)/1.04 + (n_jets>0))*wt_ff_1'
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0)/(1.22618-0.0546679*dR)*1.04 + (n_jets>0))*wt_ff_1'
                elif 'njets1' in add_name:
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0) + (n_jets>0)*(1.09292-0.0200048*dR)/1.04)*wt_ff_1'
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0) + (n_jets>0)/(1.09292-0.0200048*dR)*1.04)*wt_ff_1'
                else:
                  if 'Up' in add_name:
                    ff_syst_weight = '((n_jets==0)*(1.22618-0.0546679*dR)/1.04 + (n_jets>0)*(1.09292-0.0200048*dR)/1.04)*wt_ff_1'
                  if 'Down' in add_name:
                    ff_syst_weight = '((n_jets==0)/(1.22618-0.0546679*dR)*1.04 + (n_jets>0)/(1.09292-0.0200048*dR)*1.04)*wt_ff_1'
              fake_factor_wt_string_1 = ff_syst_weight

            if options.analysis in ['cpprod','cpdecay','mssmrun2','vlq']:
              fake_factor_wt_string_2='0'
        else:
          if options.era in ["smsummer16","cpsummer16","cpdecay16","legacy16",'UL_16_preVFP','UL_16_postVFP',"cpsummer17",'UL_17',"mvadm2016","cp18",'UL_18']:
            # deep tau tight 2018 anti isolating the subleading tau
            if options.analysis == 'cpprod':
              fake_factor_wt_string_2='0'
              fake_factor_wt_string_1 = "wt_ff_dmbins_1"
            elif options.analysis == 'cpdecay':
              fake_factor_wt_string_2='0'
              fake_factor_wt_string_1 = "wt_ff_1"
            elif options.analysis in ['mssmrun2','vlq']:
              fake_factor_wt_string_2='0'
              if options.wp == 'medium':
                fake_factor_wt_string_1 = "wt_ff_mssm_1"
                if options.ml_ff:
                  fake_factor_wt_string_1 = "0.5*wt_ff_reweight_qcd_1*((os==0) + (wt_ff_reweight_qcd_dr_to_ar_1*(os==1)))"
                  fake_factor_wt_string_2 = "0.5*wt_ff_reweight_qcd_2*((os==0) + (wt_ff_reweight_qcd_dr_to_ar_2*(os==1)))"
              elif options.wp == 'tight':
                fake_factor_wt_string_1 = "wt_ff_mssm_tight_1"
                #fake_factor_wt_string_1 = '((n_prebjets==0 && jet_pt_1<1.25*pt_1)*((pt_1<200)*(15.4087*TMath::Landau(min(pt_1,199.),-15.7496,4.82075)+0.0870211) + (pt_1>=200)*0.27557) + (n_prebjets==0 && jet_pt_1>=1.25*pt_1&&jet_pt_1<1.5*pt_1)*((pt_1<200)*(-411615*TMath::Landau(min(pt_1,199.),-110.218,-14.0548)+0.084284) + (pt_1>=200)*0.28273) + (n_prebjets==0 &&jet_pt_1>=1.5*pt_1)*(0.0392299) + (n_prebjets>0&&jet_pt_1<1.25*pt_1)*((pt_1<200)*(11.7652*TMath::Landau(min(pt_1,199.),-12.9921,5.06968)+0.077124) + (pt_1>=200)*0.13219) + (n_prebjets>0&&jet_pt_1>=1.25*pt_1&&jet_pt_1<1.5*pt_1)*((pt_1<200)*(144.787*TMath::Landau(min(pt_1,199.),14.249,0.467844)+0.0529324) + (pt_1>=200)*0.07516) + (n_prebjets>0&&jet_pt_1>=1.5*pt_1)*(-51.0159*TMath::Landau(min(pt_1,199.),-142.619,-346.505)+0.0294523))*((n_deepbjets==0)*((0.950911+-0.05705*min(dR,5.)+0.0159116*pow(min(dR,5.),2)+0.00199494*pow(min(dR,5.),3))) + (n_deepbjets>0)*((2.53701+-2.23664*min(dR,5.)+0.87535*pow(min(dR,5.),2)+-0.101159*pow(min(dR,5.),3))))'

          else:    
            fake_factor_wt_string_1 = "wt_ff_"+options.cat+"_1"
            fake_factor_wt_string_2 = "wt_ff_"+options.cat+"_2"
            fake_factor_wt_string_1+='*wt_tau1_id_loose'
            fake_factor_wt_string_2+='*wt_tau2_id_loose'
        if wt is not "": 
            wt_1=wt+"*"+fake_factor_wt_string_1
            wt_2=wt+"*"+fake_factor_wt_string_2
        else: 
            wt_1=fake_factor_wt_string_1
            wt_2=fake_factor_wt_string_2
        
        if options.analysis == 'cpdecay' and get_os:
          if (options.era=='cp18' or options.era == 'UL_18'):
            fake_factor_wt_string_1+='*wt_ff_qcd_syst_down_1/wt_ff_1*((n_jets==0)*(1.24845-0.0578994*dR) + (n_jets>0)*(1.15364-0.0317013*dR))'  
          if (options.era=='cpsummer17' or options.era == 'UL_17'):
            fake_factor_wt_string_1+='*wt_ff_qcd_syst_down_1/wt_ff_1*((n_jets==0)*(1.40415-0.101317*dR) + (n_jets>0)*(1.16253-0.031154*dR))'
          if options.era in ['legacy16','UL_16_preVFP','UL_16_postVFP']:
            fake_factor_wt_string_1+='*wt_ff_qcd_syst_down_1/wt_ff_1*((n_jets==0)*(1.22618-0.0546679*dR) + (n_jets>0)*(1.09292-0.0200048*dR))'
          #fake_factor_wt_string_1='wt_ff_test_lead_1'
          if wt is not "":
              wt_1=wt+"*"+fake_factor_wt_string_1
          else:
              wt_1=fake_factor_wt_string_1
        full_selection_1 = BuildCutString(wt_1, sel, ff_cat_1_data, OSSS, '')
        full_selection_2 = BuildCutString(wt_2, sel, ff_cat_2_data, OSSS, '')
     
        if options.ff_ss_closure:
          # usual OS FF


          ff_total_node = SummedNode('jetFakes_pre'+add_name)
          f1_total_node = SummedNode('data')
          f1_total_node.AddNode(ana.SummedFactory('data_1', data, plot_unmodified, full_selection_1))
          f1_total_node.AddNode(ana.SummedFactory('data_2', data, plot_unmodified, full_selection_2))
          f2_total_node = SummedNode('total_bkg')
          f2_total_node.AddNode(GetSubtractNode(ana,'_1',plot,plot_unmodified,wt_1+sub_wt,sel+'*(gen_match_1<6)',ff_cat_1,ff_cat_1_data,8,1.0,get_os,True))
          f2_total_node.AddNode(GetSubtractNode(ana,'_2',plot,plot_unmodified,wt_2+sub_wt,sel+'*(gen_match_2<6)',ff_cat_2,ff_cat_2_data,8,1.0,get_os,True))

          #ana.nodes[nodename].AddNode(SubtractNode('jetFakes_pre'+add_name, f1_total_node, f2_total_node))

          # FF for SS data
          full_selection_1_ss = BuildCutString(wt_1, sel, ff_cat_1_data, '!os', '')
          full_selection_2_ss = BuildCutString(wt_2, sel, ff_cat_2_data, '!os', '')
          ff_total_node_ss = SummedNode('jetFakes_ss'+add_name)
          f1_total_node_ss = SummedNode('data_ss')
          f1_total_node_ss.AddNode(ana.SummedFactory('data_1_ss', data, plot_unmodified, full_selection_1_ss))
          f1_total_node_ss.AddNode(ana.SummedFactory('data_2_ss', data, plot_unmodified, full_selection_2_ss))
          f2_total_node_ss = SummedNode('total_bkg_ss')
          f2_total_node_ss.AddNode(GetSubtractNode(ana,'_1_ss',plot,plot_unmodified,wt_1+sub_wt,sel+'*(gen_match_1<6)',ff_cat_1,ff_cat_1_data,8,1.0,False,True))
          f2_total_node_ss.AddNode(GetSubtractNode(ana,'_2_ss',plot,plot_unmodified,wt_2+sub_wt,sel+'*(gen_match_2<6)',ff_cat_2,ff_cat_2_data,8,1.0,False,True))
          den_node = SubtractNode('jetFakes_ss'+add_name, f1_total_node_ss, f2_total_node_ss)

          # SS data - bkg
          subtract_node = GetSubtractNode(ana,'',plot,plot_unmodified,wt,sel+'*(gen_match_1<6 && gen_match_2<6)',cats[cat_name]+'&&'+cats['baseline'],cats_unmodified[cat_name]+'&&'+cats['baseline'],8,1.0,False,True)
          num_selection = BuildCutString(wt, sel, cats_unmodified[cat_name]+'&&'+cats['baseline'], '!os')
          num_node = SubtractNode('ratio_num',
                     ana.SummedFactory('data', data, plot_unmodified, num_selection),
                     subtract_node)

          ana.nodes[nodename].AddNode(HttQCDNode('jetFakes'+add_name,
          f1_total_node,
          f2_total_node,
          1,
          None,
          num_node,
          den_node))
        else: 
          ff_total_node = SummedNode('jetFakes'+add_name)
          f1_total_node = SummedNode('data')
          f1_total_node.AddNode(ana.SummedFactory('data_1', data, plot_unmodified, full_selection_1))
          if options.ml_ff: f1_total_node.AddNode(ana.SummedFactory('data_2', data, plot_unmodified, full_selection_2))
          f2_total_node = SummedNode('total_bkg')
          f2_total_node.AddNode(GetSubtractNode(ana,'_1',plot,plot_unmodified,wt_1+sub_wt,sel+'*(gen_match_1<6)',ff_cat_1,ff_cat_1_data,8,1.0,get_os,True))
          if options.ml_ff: f2_total_node.AddNode(GetSubtractNode(ana,'_2',plot,plot_unmodified,wt_2+sub_wt,sel+'*(gen_match_2<6)',ff_cat_2,ff_cat_2_data,8,1.0,get_os,True))
          ana.nodes[nodename].AddNode(SubtractNode('jetFakes'+add_name, f1_total_node, f2_total_node))


        if options.channel=='tt':
          #if options.analysis == 'cpprod': full_selection_extra = BuildCutString(wt+'*wt_ff_us_2', sel+'*(gen_match_2==6)', ff_cat_2_data, OSSS, '')
          #else: full_selection_extra = BuildCutString(wt+'*wt_ff_2', sel+'*(gen_match_2==6)', ff_cat_2_data, OSSS, '')

          full_selection_extra = BuildCutString(wt, sel, cats_unmodified[cat_name]+'&&'+cats['baseline'], OSSS, 'gen_match_2==6')       
  

          #wnode = ana.SummedFactory('Wfakes'+add_name, ztt_samples+vv_samples+wjets_samples+ewkz_samples+top_samples, plot, full_selection_extra)
#          ana.nodes[nodename].AddNode(wnode) 
        
def GenerateSMSignal(ana, add_name='', plot='', masses=['125'], wt='', sel='', cat='', get_os=True, sm_bkg = '',processes=['ggH','qqH','ZH','WminusH','WplusH']):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    if masses is not None:
        for mass in masses:  
            if sm_bkg != '':
                add_str = '_SM'+sm_bkg
            else:
                add_str = mass
            for key in sm_samples:
                if 'ggH' in key: 
                  full_selection = BuildCutString(wt+'*wt_ph_nnlops', sel, cat, OSSS)
                else: full_selection = BuildCutString(wt, sel, cat, OSSS)
                if True not in [proc in key for proc in processes]: continue
                sample_names=[]
                if isinstance(sm_samples[key], (list,)):
                  for i in sm_samples[key]: 
                    sample_names.append(i.replace('*',mass))
                else: sample_names = [sm_samples[key].replace('*',mass)]
                ana.nodes[nodename].AddNode(ana.SummedFactory(key.replace('*',add_str)+add_name, sample_names, plot, full_selection))
            
def GenerateMSSMSignal(ana, add_name='', bbh_add_name='', plot='', ggh_masses = ['1000'], bbh_masses = ['1000'], wt='', sel='', cat='', get_os=True, do_ggH=True, do_bbH=True):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    if options.gen_signal: OSSS='1'
    full_selection = BuildCutString(wt, sel, cat, OSSS)
    for key in mssm_samples:
        masses = None
        if key == 'ggH':
            masses = ggh_masses
        elif key == 'bbH':
            masses = bbh_masses
        if masses is not None:    
            for mass in masses:
                if key == 'ggH' and not do_ggH:
                    continue
                if key == 'bbH' and not do_bbH:
                    continue
                sample_name = mssm_samples[key].replace('*',mass)
                add_name_2 = ''
                if bbh_add_name == '-LO' and key is 'bbH': 
                    sample_name = mssm_lo_samples['bbH-LO'].replace('*',mass)
                    add_name_2 = bbh_add_name
                ana.nodes[nodename].AddNode(ana.BasicFactory(key+add_name_2+mass+add_name, sample_name, plot, full_selection))

# for CP signals
# need to multiply by the weights wt_cp_sm/ps/mm
def GenerateReweightedCPSignal(ana, add_name='', plot='', wt='', sel='', cat='', get_os=True):
    weights = {"sm": "wt_cp_sm*wt_ph_nnlops", "ps": "wt_cp_ps*wt_ph_nnlops", "mm": "wt_cp_mm*wt_ph_nnlops", "flat": "wt_ph_nnlops"}
    if get_os: 
        OSSS = 'os'
    else: 
        OSSS = '!os'
    if options.gen_signal: 
        OSSS='1' 
    for key, sample in sm_samples.iteritems():
        non_cp = True 
        for name in weights:
            if key.split("_")[1] == name:
                non_cp=False
                weight=wt+"*"+weights[name]
                full_selection = BuildCutString(weight, sel, cat, OSSS)
                name = key

                sample_names=[]
                if isinstance(sm_samples[key], (list,)):
                  for i in sm_samples[key]:
                    sample_names.append(i.replace('*',mass))
                else: sample_names = [sm_samples[key].replace('*',mass)]
                ana.nodes[nodename].AddNode(ana.SummedFactory(key.replace('*',mass)+add_name, sample_names, plot, full_selection))
                #ana.nodes[nodename].AddNode(ana.BasicFactory(name+mass+add_name, sample, plot, full_selection))
        if non_cp:
             full_selection = BuildCutString(wt, sel, cat, OSSS)
             name = key
             #ana.nodes[nodename].AddNode(ana.BasicFactory(name+mass+add_name, sample, plot, full_selection))

             sample_names=[]
             if isinstance(sm_samples[key], (list,)):
               for i in sm_samples[key]:
                 sample_names.append(i.replace('*',mass))
             else: sample_names = [sm_samples[key].replace('*',mass)]
             ana.nodes[nodename].AddNode(ana.SummedFactory(key.replace('*',mass)+add_name, sample_names, plot, full_selection))

def GenerateReweightedCPProdSignal(ana, add_name='', plot='', wt='', sel='', cat='', get_os=True):
    weights = {"sm": "wt_cp_prod_sm", "ps": "wt_cp_prod_ps", "mm": "wt_cp_prod_mm"}
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    if options.gen_signal:
        OSSS='1'
    for key, sample in sm_samples.iteritems():
        non_cp = True
        if 'ggH_sm' in key or 'ggH_ps' in key or 'ggH_mm' in key: non_cp=False
        if not non_cp:
          for name in weights: 
            tname = key.replace('*',mass)+add_name
            if 'ggH_'+name not in tname: 
              tname=key.replace('*',mass)+'_reweightedto_'+name+add_name
              #continue
              weight=wt+"*"+weights[name]+'*wt_quarkmass*wt_mg_nnlops'
            else:
              weight=wt+'*wt_quarkmass*wt_mg_nnlops*(wt_cp_prod_sm!=0)'
            non_cp=False
            #weight=wt+"*"+weights[name]+'*wt_quarkmass*wt_mg_nnlops'
            full_selection = BuildCutString(weight, sel, cat, OSSS)
            name = key

            sample_names=[]
            if isinstance(sm_samples[key], (list,)):
              for i in sm_samples[key]:
                sample_names.append(i.replace('*',mass))
            else: sample_names = [sm_samples[key].replace('*',mass)]
            ana.nodes[nodename].AddNode(ana.SummedFactory(tname, sample_names, plot, full_selection))
        if non_cp:
             name = key
             weight=wt
             if 'ggh' in name.lower():
               if 'reweighted_' in name: weight+='*wt_quarkmass*wt_mg_nnlops'
               else: weight+='*wt_ph_nnlops'
             full_selection = BuildCutString(weight, sel, cat, OSSS)
             sample_names=[]
             if isinstance(sm_samples[key], (list,)):
               for i in sm_samples[key]:
                 sample_names.append(i.replace('*',mass))
             else: sample_names = [sm_samples[key].replace('*',mass)]
             ana.nodes[nodename].AddNode(ana.SummedFactory(key.replace('*',mass)+add_name, sample_names, plot, full_selection))


def GenerateReWeightedMSSMSignal(ana, add_name='', plot='', ggh_masses = ['1000'], wt='', sel='', cat='', get_os=True):
  weights = {'ggh_t_':'wt_ggh_t', 'ggh_b_':'wt_ggh_b', 'ggh_i_':'wt_ggh_i', 'ggH_t_':'wt_ggH_t', 'ggH_b_':'wt_ggH_b', 'ggH_i_':'wt_ggH_i', 'ggA_t_':'wt_ggA_t', 'ggA_b_':'wt_ggA_b', 'ggA_i_':'wt_ggA_i' }  
  #weights = {'ggH_t_':'wt_ggH_t', 'ggH_b_':'wt_ggH_b', 'ggH_i_':'wt_ggH_i', 'ggA_t_':'wt_ggA_t', 'ggA_b_':'wt_ggA_b', 'ggA_i_':'wt_ggA_i' } 
  if get_os: OSSS = 'os'
  else: OSSS = '!os'
  if options.gen_signal: OSSS='1' 
  for mass in ggh_masses:
    if int(mass) == 95:
      full_selection = BuildCutString(wt, sel, cat, OSSS)
      sample_name = mssm_samples['ggH'].replace('*',mass)
      name = 'ggH'
      ana.nodes[nodename].AddNode(ana.BasicFactory(name+mass+add_name, sample_name, plot, full_selection))
    for name in weights:
      weight=wt+"*"+weights[name]    
      full_selection = BuildCutString(weight, sel, cat, OSSS)
      sample_name = mssm_samples['ggH'].replace('*',mass)
      ana.nodes[nodename].AddNode(ana.BasicFactory(name+mass+add_name, sample_name, plot, full_selection))
      if options.syst_mssm_ggh and add_name=='':
        for u in ['_scale', '_hdamp']:
          weight_up=wt+"*"+weights[name].replace('ggH','ggh')+u+'_up'
          weight_down=wt+"*"+weights[name].replace('ggH','ggh')+u+'_down'
          full_selection_up = BuildCutString(weight_up, sel, cat, OSSS)
          full_selection_down = BuildCutString(weight_down, sel, cat, OSSS)
          if 'scale' in u:
            syst_name_up='_QCDscale_ggH_REWEIGHTUp'
            syst_name_down='_QCDscale_ggH_REWEIGHTDown'
          else:
            syst_name_up = '_Hdamp_%sREWEIGHTUp' % (name.replace('ggA','ggH').replace('ggh','ggH'))
            syst_name_down = '_Hdamp_%sREWEIGHTDown' % (name.replace('ggA','ggH').replace('ggh','ggH'))
          sample_name = mssm_samples['ggH'].replace('*',mass)
          ana.nodes[nodename].AddNode(ana.BasicFactory(name+mass+add_name+syst_name_up, sample_name, plot, full_selection_up))      
          ana.nodes[nodename].AddNode(ana.BasicFactory(name+mass+add_name+syst_name_down, sample_name, plot, full_selection_down))      
          
def GenerateNLOMSSMSignal(ana, add_name='', plot='', ggh_nlo_masses = ['1000'], bbh_nlo_masses = ['1000'],wt='wt', sel='', cat='', doScales=True, doPDF=False, get_os=True,do_ggH=True, do_bbH=True):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    if options.gen_signal: OSSS='1'
    weights = {'':'1'}
    wt_noscale = wt
    if doScales: weights = {'':'1','muR1muF2':'wt_mur1_muf2','muR1muF0.5':'wt_mur1_muf0p5','muR2muF1':'wt_mur2_muf1','muR2muF2':'wt_mur2_muf2','muR0.5muF1':'wt_mur0p5_muf1','muR0.5muF0.5':'wt_mur0p5_muf0p5'}
    if doPDF:
      for i in range(1,101): weights['PDF_'+str(i)] = 'wt_pdf_'+str(i)
      weights['AlphaS_Up'] = 'wt_alphasup'
      weights['AlphaS_Down'] = 'wt_alphasdown'
    for weight in weights:
      wt = weights[weight]+'*'+wt_noscale  
      full_selection = BuildCutString(wt, sel, cat, OSSS)
      for key in mssm_nlo_samples:
          if 'Qsh' in key and weight is not '': continue
          if 'ggH' in key:
              masses = ggh_nlo_masses
          elif 'bbH' in key:
              masses = bbh_nlo_masses
          if masses is not None:    
              for mass in masses:
                  if key == 'ggH' and not do_ggH:
                      continue
                  if key == 'bbH' and not do_bbH:
                      continue
                  sample_name = mssm_nlo_samples[key].replace('*',mass)
                  ana.nodes[nodename].AddNode(ana.BasicFactory(key+mass+add_name+weight, sample_name, plot, full_selection))
        
        
def GenerateHhhSignal(ana, add_name='', plot='', masses = ['700'], wt='', sel='', cat='', get_os=True):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS)
    if masses is not None:
        for mass in masses:
            for key in Hhh_samples:
                sample_name = Hhh_samples[key].replace('*',mass)
                ana.nodes[nodename].AddNode(ana.BasicFactory(key+mass+add_name, sample_name, plot, full_selection))

def GenerateVLQSignal(ana, add_name='', plot='', vlq_sig= [] ,wt='', sel='', cat='', get_os=True):
    if get_os:
        OSSS = 'os'
    else:
        OSSS = '!os'
    full_selection = BuildCutString(wt, sel, cat, OSSS)
    for key in vlq_samples:
      if key in vlq_sig:
        if isinstance(vlq_samples[key], (list,)):
          sample_names = []
          for i in vlq_samples[key]:
            sample_names.append(i)
        else: sample_names = [vlq_samples[key]]
        ana.nodes[nodename].AddNode(ana.SummedFactory(key+add_name, sample_names, plot, full_selection))

 
def PrintSummary(nodename='', data_strings=['data_obs'], add_names=''):
    print ''
    print '################### Summary ###################'
    nodes = ana.nodes[nodename].SubNodes()
    bkg_total = ufloat(0.000000001,0.000000001)
    sig_total = ufloat(0.000000001,0.000000001)
    for node in nodes:
        if options.channel == 'em' and node.name == 'W': continue
        if node.shape.rate.n == 0: per_err = 0
        else: per_err = node.shape.rate.s/node.shape.rate.n
        print node.name.ljust(10) , ("%.2f" % node.shape.rate.n).ljust(10), '+/-'.ljust(5), ("%.2f" % node.shape.rate.s).ljust(7), "(%.4f)" % per_err 
        if True in [node.name.find(add_name) != -1 and add_name is not '' for add_name in add_names]: continue
        if len(signal_samples) != 0: sig_samp_cond = [node.name.find(sig) != -1 for sig in signal_samples.keys()]
        else: sig_samp_cond = []
        if True in sig_samp_cond and node.name.find("_SM"+options.add_sm_background) ==-1:
            sig_total += node.shape.rate
        elif node.name not in data_strings or (options.method == 18 and 'jetFakes' == node.name):
            bkg_total += node.shape.rate      
    if bkg_total.n == 0: per_err = 0        
    else: per_err = bkg_total.s/bkg_total.n
    print 'Total bkg'.ljust(10) , ("%.2f" % bkg_total.n).ljust(10), '+/-'.ljust(5), ("%.2f" % bkg_total.s).ljust(7), "(%.4f)" % per_err
    if sig_total.n == 0: per_err = 0
    else: per_err = sig_total.s/sig_total.n
    print 'Total sig'.ljust(10) , ("%.2f" % sig_total.n).ljust(10), '+/-'.ljust(5), ("%.2f" % sig_total.s).ljust(7), "(%.4f)" % per_err
    print '###############################################'
    print ''
    

def FixBins(ana,outfile='output.root'):
    #Fix empty histograms
    nodes = ana.nodes[nodename].SubNodes()
    for node in nodes:
        if 'data_obs' in node.name: continue
        hist = outfile.Get(nodename+'/'+node.name)
        outfile.cd(nodename)
        #Fix empty histogram
        if hist.Integral() == 0.0:
            hist.SetBinContent(hist.GetNbinsX()/2, 0.00001)
            hist.SetBinError(hist.GetNbinsX()/2, 0.00001)
            hist.Write(hist.GetName(),ROOT.TObject.kWriteDelete)
        outfile.cd()

                                                                                                                                
def NormFFSysts(ana,outfile='output.root'):
    nominal_hist = outfile.Get(nodename+'/jetFakes')
    if isinstance(nominal_hist,ROOT.TH2): nominal_scale = nominal_hist.Integral(-1,-1,-1,-1)
    else: nominal_scale = nominal_hist.Integral(0,nominal_hist.GetNbinsX()+2)
    directory = outfile.Get(nodename)
    outfile.cd(nodename)
    hists_to_add=[]
    for key in directory.GetListOfKeys():
        hist_name = key.GetName()
        hist = directory.Get(hist_name).Clone()
        if not isinstance(hist,ROOT.TDirectory):
           if 'jetFakes' not in hist_name: continue
           if hist_name == 'jetFakes': continue
           if isinstance(hist,ROOT.TH2): norm = nominal_scale/hist.Integral(-1,-1,-1,-1)
           else: norm = nominal_scale/hist.Integral(0,hist.GetNbinsX()+2)
           hist.Scale(norm)
           norm_hist_name = hist_name
           norm_hist_name = norm_hist_name.replace('jetFakes','jetFakes_norm')
           hist.SetName(norm_hist_name)
           hists_to_add.append(hist)
    for hist in hists_to_add: hist.Write()
    
def NormWFakeSysts(ana,outfile='output.root'):
    nominal_hist = outfile.Get(nodename+'/W')
    if isinstance(nominal_hist,ROOT.TH2): nominal_scale = nominal_hist.Integral(-1,-1,-1,-1)
    else: nominal_scale = nominal_hist.Integral(-1,-1)
    directory = outfile.Get(nodename)
    outfile.cd(nodename)
    hists_to_add=[]
    for key in directory.GetListOfKeys():
        hist_name = key.GetName()
        hist = directory.Get(hist_name)
        if not isinstance(hist,ROOT.TDirectory):
           if 'W' not in hist_name or options.syst_w_fake_rate not in hist_name or 'Wplus' in hist_name or 'Wminus' in hist_name or 'WH' in hist_name: continue
           if isinstance(hist,ROOT.TH2): norm = nominal_scale/hist.Integral(-1,-1,-1,-1)
           else: norm = nominal_scale/hist.Integral(-1,-1)
           hist.Scale(norm)
           hists_to_add.append(hist)
    for hist in hists_to_add: hist.Write("",ROOT.TObject.kOverwrite)

def NormEmbedToMC(ana,outfile='output.root'):
    nominal_hist = outfile.Get(nodename+'/ZTT')
    nominal_hist_embed = outfile.Get(nodename+'/EmbedZTT')
    if isinstance(nominal_hist,ROOT.TH2): scale = nominal_hist.Integral(-1,-1,-1,-1)/nominal_hist_embed.Integral(-1,-1,-1,-1)
    else: scale = nominal_hist.Integral(-1,-1)/nominal_hist_embed.Integral(-1,-1)
    directory = outfile.Get(nodename)
    outfile.cd(nodename)
    hists_to_add=[]
    for key in directory.GetListOfKeys():
        hist_name = key.GetName()
        hist = directory.Get(hist_name)
        if not isinstance(hist,ROOT.TDirectory):
           if 'EmbedZTT' not in hist_name: continue
           hist.Scale(scale)
           hists_to_add.append(hist)
    for hist in hists_to_add: hist.Write("",ROOT.TObject.kOverwrite)

    
def TTBarEmbeddingSyst(ana,outfile,template_name):
    nominal_hist = outfile.Get(nodename+'/EmbedZTT')    
    shift_hist = outfile.Get(nodename+'/TTT_embed_syst')
    shift_hist_2 = outfile.Get(nodename+'/VVT_embed_syst')
    shift_hist.Add(shift_hist_2)
    shift_hist.Scale(0.1)
    up_hist = nominal_hist.Clone()
    down_hist = nominal_hist.Clone()
    up_hist.Add(shift_hist)
    down_hist.Add(shift_hist,-1)
    up_hist.SetName('EmbedZTT_'+template_name+'Up')
    down_hist.SetName('EmbedZTT_'+template_name+'Down')
    outfile.cd(nodename)
    up_hist.Write()
    down_hist.Write()
    
def OverwriteNames(input_string,replace_dict={}):
    for i in replace_dict:
      input_string=input_string.replace(i,replace_dict[i])
    return input_string
    
def PDFUncerts(nodename, infile):
  def RMS(a):
    from numpy import mean, sqrt, square  
    rms = sqrt(mean(square(a-mean(a))))
    return rms
  outstring1=''
  outstring2=''
  
  for mass in bbh_nlo_masses:
    nominal_error=ROOT.Double()
    nominal = outfile.Get(nodename+'/bbH'+mass).IntegralAndError(-1, -1,nominal_error)
    sample_name='SUSYGluGluToBBHToTauTau_M-'+mass+'-NLO'
    evt_nom = ana.info[sample_name]['evt']
    pdf_variations_nosf=[]
    pdf_variations=[]
    pdf_variations_nosf.append(nominal)
    pdf_variations.append(nominal)
    for i in range(1,101):
      wt = 'wt_pdf_'+str(i)
      name = 'PDF_'+str(i)
      val = outfile.Get(nodename+'/bbH'+mass+name).Integral(-1, -1)
      pdf_variations_nosf.append(val)
      evt_var = ana.info[sample_name]['evt_'+wt]
      sf = evt_nom/evt_var
      pdf_variations.append(val*sf)
    pdf_uncert_nosf=RMS(pdf_variations_nosf)/nominal*100
    pdf_uncert=RMS(pdf_variations)/nominal*100
    #print pdf_uncert_nosf, pdf_uncert
    
    outstring1+=mass+','
    outstring2+=str(pdf_uncert)+','
    
    nominal_error=ROOT.Double()
    nominal = outfile.Get(nodename+'/bbH'+mass).IntegralAndError(-1, -1,nominal_error)
    alphas_down_error_nosf=ROOT.Double()
    alphas_down_nosf = outfile.Get(nodename+'/bbH'+mass+'AlphaS_Down').IntegralAndError(-1, -1,alphas_down_error_nosf)
    alphas_up_error_nosf=ROOT.Double()
    alphas_up_nosf = outfile.Get(nodename+'/bbH'+mass+'AlphaS_Up').IntegralAndError(-1, -1,alphas_up_error_nosf)
    evt_var = ana.info[sample_name]['evt_wt_alphasdown']
    sf = evt_nom/evt_var
    alphas_down=alphas_down_nosf*sf
    alphas_down_error=alphas_down_error_nosf*sf
    evt_var = ana.info[sample_name]['evt_wt_alphasup']
    sf = evt_nom/evt_var
    alphas_up=alphas_up_nosf*sf
    alphas_up_error=alphas_up_error_nosf*sf
    
    alphas_uncert = (alphas_up-alphas_down)/2/nominal
    alphas_uncert_error = math.sqrt(alphas_down_error**2+alphas_up_error**2)/(alphas_up-alphas_down)*alphas_uncert
    #(alphas_up_nosf-alphas_down_nosf)/2/nominal
    #print alphas_uncert*100, '\\% $\\pm$', alphas_uncert_error*100,'\\%'
  print outstring1
  print outstring2

def DONLOUncerts(nodename,infile):
    def LargestDiff(nominal,scales_shifted):
        largest_diff=0
        value = nominal
        for i in scales_shifted:
            diff = abs(scales_shifted[i] - nominal)
            if diff > largest_diff: 
              largest_diff = diff
              value = scales_shifted[i]
        return value
    if not options.bbh_nlo_masses: return
    outstring='%'+options.channel+' '+options.datacard+'\n'
    if options.nlo_qsh: outstring+='\\begin{table}[H]\n\\centering\n\\resizebox{\\textwidth}{!}{\n\\begin{tabular}{ |c|c|c| }\n\\hline\nSignal Mass (GeV) & Qsh Uncertainty &  Qsh Uncertainty (*)'
    else: outstring+='\\begin{table}[H]\n\\centering\n\\resizebox{\\textwidth}{!}{\n\\begin{tabular}{ |c|c|c| }\n\\hline\nSignal Mass (GeV) & Scale Uncertainty &  Scale Uncertainty (*)'
    outstring += '\\\\\n\\hline\n'
    #outstring2='{'
    #outstring3='{'
    #outstring4='{'
    for mass in bbh_nlo_masses:
      nominal_error=ROOT.Double()
      nominal = outfile.Get(nodename+'/bbH'+mass).IntegralAndError(-1, -1,nominal_error) 
      samples = {'bbH*':'', 'bbH*muR0.5muF0.5':'wt_mur0p5_muf0p5', 'bbH*muR1muF0.5':'wt_mur1_muf0p5', 'bbH*muR0.5muF1':'wt_mur0p5_muf1', 'bbH*muR2muF2':'wt_mur2_muf2', 'bbH*muR2muF1':'wt_mur2_muf1', 'bbH*muR1muF2':'wt_mur1_muf2'} 
      qsh_down_error=ROOT.Double()
      qsh_up_error=ROOT.Double()
      if options.nlo_qsh:
        qsh_down = outfile.Get(nodename+'/bbH-QshDown'+mass).IntegralAndError(-1, -1,qsh_down_error) 
        qsh_up = outfile.Get(nodename+'/bbH-QshUp'+mass).IntegralAndError(-1, -1,qsh_up_error)
        qsh_uncert_1=(max(nominal,qsh_down,qsh_up) - min(nominal,qsh_down,qsh_up))/2
        up_error = nominal_error 
        down_error = nominal_error
        if max(nominal,qsh_down,qsh_up) is qsh_up: up_error = qsh_up_error
        if max(nominal,qsh_down,qsh_up) is qsh_down: up_error = qsh_down_error
        if min(nominal,qsh_down,qsh_up) is qsh_up: down_error = qsh_up_error
        if min(nominal,qsh_down,qsh_up) is qsh_down: down_error = qsh_down_error
        qsh_error_1 = math.sqrt(up_error**2 + down_error**2)
        qsh_uncert_2 = (qsh_up - qsh_down)/2
        qsh_error_2 = math.sqrt(qsh_up_error**2 + qsh_down_error**2)
      scale_max = nominal
      scale_min = nominal
      scale_nosf_max = nominal
      scale_nosf_min = nominal
      up_dic = {}
      down_dic = {}
      for samp in samples:
        acceptance_error = ROOT.Double()  
        acceptance = outfile.Get(nodename+'/'+samp.replace('*',mass)).IntegralAndError(-1, -1,acceptance_error)
        if samp is 'bbH*': sf = 1.0 
        else: 
          sample_name='SUSYGluGluToBBHToTauTau_M-'+mass+'-NLO'
          evt_nom = ana.info[sample_name]['evt']
          evt_var = ana.info[sample_name]['evt_'+samples[samp]]
          sf = evt_nom/evt_var
        acceptance_nosf = acceptance
        acceptance*=sf
        #if samples[samp] in ['wt_mur0p5_muf0p5', 'wt_mur1_muf0p5', 'wt_mur0p5_muf1']: down_dic[samples[samp]] = acceptance
        #if samples[samp] in ['wt_mur2_muf2','wt_mur2_muf1','wt_mur1_muf2']: up_dic[samples[samp]] = acceptance
        if samples[samp] in ['wt_mur0p5_muf0p5']: 
            down_dic[samples[samp]] = [acceptance,acceptance_error]
        if samples[samp] in ['wt_mur2_muf2']: up_dic[samples[samp]] = [acceptance,acceptance_error]
        if acceptance > scale_max: 
            scale_max = acceptance
            up_error = acceptance_error
        if acceptance < scale_min:
            scale_min = acceptance
            down_error = acceptance_error
        if acceptance_nosf > scale_nosf_max: scale_nosf_max = acceptance_nosf
        if acceptance_nosf < scale_nosf_min: scale_nosf_min = acceptance_nosf
      #up_nom = LargestDiff(nominal,up_dic)
      #down_nom = LargestDiff(nominal,down_dic)
      up_nom = up_dic['wt_mur2_muf2'][0]
      down_nom = down_dic['wt_mur0p5_muf0p5'][0]
      uncert = (scale_max-scale_min)/2
      uncert_error = math.sqrt(up_error**2+down_error**2)/(scale_max-scale_min)*uncert
      uncert_nosf = (scale_nosf_max-scale_nosf_min)/2
      uncert_alt_method = (up_nom-down_nom)/2
      uncert_alt_error = math.sqrt(up_dic['wt_mur2_muf2'][1]**2 + down_dic['wt_mur0p5_muf0p5'][1]**2)/(up_nom-down_nom) *uncert_alt_method
      pythia_error=ROOT.Double()
      pythia_yield = outfile.Get(nodename+'/bbH'+mass).IntegralAndError(-1, -1,pythia_error) 
      outstring +=mass#+ ' & '+ str(round(pythia_yield,1))+' $\pm$ '+str(round(pythia_error,1))+ ' & '+ str(round(nominal,1))+' $\pm$ '+str(round(nominal_error,1))+ '('+str(round((pythia_yield-nominal)*100/pythia_yield,2))+' \%)' + ' & '+ str(round(uncert_nosf/nominal,2)) + ' & '+ str(round(uncert/nominal,2))+ ' & '+ str(round(uncert_alt_method/nominal,2))
      if options.nlo_qsh: 
        outstring+=' & '+ str(round(100*qsh_uncert_1/nominal,1))+' $\pm$ '+str(round(100*qsh_error_1/nominal,1))+' & '+ str(round(100*qsh_uncert_2/nominal,1))+' $\pm$ '+str(round(100*qsh_error_2/nominal,1))+'\\\\\n'  
      else:
        outstring+=' & '+ str(round(100*uncert/nominal,1))+'\\% $\\pm$ '+ str(round(100*uncert_error/nominal,1)) +'\\% & '+ str(round(100*uncert_alt_method/nominal,1))+'\\% $\\pm$ '+str(round(100*uncert_alt_error/nominal,1)) +'\\% \\\\\n'
      #outstring2 +=str(round(100*uncert/nominal,1))+','
      #outstring3 +=str(round(100*uncert_error/nominal,1))+','
      #outstring4+=mass+','
    outstring+='\\hline\n\\end{tabular}}\n\\end{table}'
    print outstring
    #print outstring2
    #print outstring3
    #print outstring4
    
def ScaleUncertBand(nodename='',outfile='output.root',NormScales=True):
    hist_names=['bbH*muR0.5muF0.5','bbH*muR1muF0.5','bbH*muR0.5muF1','bbH*muR2muF1','bbH*muR1muF2','bbH*muR2muF2']
    mass = '100'
    if options.draw_signal_mass: mass = options.draw_signal_mass
    hists=[]
    for hist_name in hist_names:
        sf = 1.0
        if NormScales:
            sample_name='SUSYGluGluToBBHToTauTau_M-'+mass+'-NLO'
            evt_nom = ana.info[sample_name]['evt']
            evt_var = ana.info[sample_name]['evt_'+hist_name.replace('bbH*','wt_').replace('muF','_muf').replace('muR','mur').replace('0.5','0p5')]
            sf = evt_nom/evt_var
        hist_name=hist_name.replace('*',mass)
        hist = outfile.Get(nodename+'/'+hist_name).Clone()
        if NormScales: hist.Scale(sf)
        hists.append(hist)
    nom_hist = outfile.Get(nodename+'/bbH'+mass)
    up_hist = nom_hist.Clone()
    down_hist = nom_hist.Clone()
    up_hist.SetName('ScaleUp')
    down_hist.SetName('ScaleDown')
    for i in range (1,nom_hist.GetNbinsX()+2):
        for hist in hists:
          max_content = up_hist.GetBinContent(i)
          min_content = down_hist.GetBinContent(i)
          content = hist.GetBinContent(i)
          if content > max_content: up_hist.SetBinContent(i,content)
          if content < min_content: down_hist.SetBinContent(i,content)
    outfile.cd(nodename)
    up_hist.Write()
    down_hist.Write()
    outfile.cd()
        
    
def DYUncertBand(outfile='output.root',ScaleToData=True):
    bkg_hist = outfile.Get(nodename+'/total_bkg')
    nominal_hist = outfile.Get(nodename+'/ZLL')
    up_hist = outfile.Get(nodename+'/total_bkg').Clone()
    down_hist = outfile.Get(nodename+'/total_bkg').Clone()
    up_hist.SetName('total_bkg_up')
    down_hist.SetName('total_bkg_down')
    shifts=['_ES', '_TT', '_Stat0', '_Stat40', '_Stat80']
    for i in range(1,nominal_hist.GetNbinsX()+2):
      nom_content = nominal_hist.GetBinContent(i)
      bkg_content = bkg_hist.GetBinContent(i)
      uncert=0
      for shift in shifts:
          shift_hist_up = outfile.Get(nodename+'/ZLL'+shift+'Up')
          shift_hist_down = outfile.Get(nodename+'/ZLL'+shift+'Down')
          up = abs(shift_hist_up.GetBinContent(i) - nom_content)
          down = abs(shift_hist_down.GetBinContent(i) - nom_content)
          uncert=math.sqrt(max([up,down])**2+uncert**2)
      up_hist.SetBinContent(i, bkg_content+uncert)
      down_hist.SetBinContent(i, bkg_content-uncert)
    outfile.cd(nodename)
    up_hist.Write()
    down_hist.Write()
    outfile.cd()
    if ScaleToData:
      data_hist=outfile.Get(nodename+'/data_obs')
      data_total=data_hist.Integral(-1,-1)
      bkg_total=bkg_hist.Integral(-1,-1)
      data_hist.Scale(bkg_total/data_total)
      outfile.cd(nodename)
      data_hist.Write()
      outfile.cd()
    
def GetTotals(ana,add_name="",outfile='outfile.root'):
    # add histograms to get totals for backgrounds split into real/fake taus and make a total backgrounds histogram
    outfile.cd(nodename)
    nodes = ana.nodes[nodename].SubNodes()
    nodenames=[]
    for node in nodes: nodenames.append(node.name)
   # for i in ['TT', 'VV', 'Z']:
   #     j = 'T'
   #     outname = i+add_name
   #     first_hist=True
   #     if options.channel == 'em' and i is 'Z':
   #         if first_hist and 'ZLL'+add_name in nodenames: 
   #             sum_hist = ana.nodes[nodename].nodes['ZLL'+add_name].shape.hist.Clone()
   #             first_hist=False
   #         elif 'ZLL'+add_name in nodenames: sum_hist.Add(ana.nodes[nodename].nodes['ZLL'+add_name].shape.hist.Clone())
   #         if first_hist and'ZTT'+add_name in nodenames: 
   #             sum_hist = ana.nodes[nodename].nodes['ZTT'+add_name].shape.hist.Clone()
   #             first_hist=False
   #         elif 'ZTT'+add_name in nodenames: sum_hist.Add(ana.nodes[nodename].nodes['ZTT'+add_name].shape.hist.Clone())
   #         if not first_hist:
   #             sum_hist.SetName(outname)
   #             sum_hist.Write()
   #     elif (options.channel == 'zee' or options.channel == 'zmm') and i is 'Z':
   #         if first_hist and 'ZLL'+add_name in nodenames:
   #             sum_hist = ana.nodes[nodename].nodes['ZLL'+add_name].shape.hist.Clone()
   #             first_hist=False
   #         elif 'ZLL'+add_name in nodenames: sum_hist.Add(ana.nodes[nodename].nodes['ZLL'+add_name].shape.hist.Clone())
   #         if not first_hist:
   #             sum_hist.SetName(outname)
   #             sum_hist.Write()
   #     else: 
   #         if i is 'Z':
   #             outname = 'ZLL'+add_name
   #             j = 'L'
   #         if i+'J' or i+j in [node.name for node in nodes]:
   #             if first_hist and i+'J'+add_name in nodenames: 
   #                 sum_hist = ana.nodes[nodename].nodes[i+'J'+add_name].shape.hist.Clone()
   #                 first_hist=False
   #             elif i+'J'+add_name in nodenames: sum_hist.Add(ana.nodes[nodename].nodes[i+'J'+add_name].shape.hist.Clone())
   #             if first_hist and i+j+add_name in nodenames:
   #                 sum_hist = ana.nodes[nodename].nodes[i+j+add_name].shape.hist.Clone()    
   #                 first_hist=False
   #             elif i+j+add_name in nodenames: sum_hist.Add(ana.nodes[nodename].nodes[i+j+add_name].shape.hist.Clone())
   #         if not first_hist:    
   #             sum_hist.SetName(outname)
   #             sum_hist.Write()
    first_hist=True
    for node in nodes:
        if len(signal_samples) != 0: sig_samp_cond = [node.name.find(sig) != -1 for sig in signal_samples.keys()]
        else: sig_samp_cond = []
        if True not in sig_samp_cond and node.name != 'data_obs' and node.name.find("_SM"+options.add_sm_background) ==-1 and not node.name.endswith('Up') and not node.name.endswith('Down'):
            if options.embedding and node.name.startswith('ZTT'): continue
            if 'embed_syst' in node.name: continue
            if options.method == 18 and 'jetFakes' == node.name: continue
            if options.channel == 'em' and node.name == 'W': continue
            if add_name not in node.name: continue
            if first_hist:
                total_bkg = ana.nodes[nodename].nodes[node.name].shape.hist.Clone()
                first_hist=False
            else: total_bkg.Add(ana.nodes[nodename].nodes[node.name].shape.hist.Clone())
    if not first_hist:        
        total_bkg.SetName('total_bkg'+add_name)
        total_bkg.Write()
    outfile.cd()
     
    #hist = outfile.Get(nodename+'/total_bkg')
    #qcd_hist = outfile.Get(nodename+'/QCD')
    #ztt_hist = outfile.Get(nodename+'/ZTT')
    #tot_sub = total_bkg.Integral(-1,-1) - qcd_hist.Integral(-1,-1) - ztt_hist.Integral(-1,-1)
    #out_str=''
    #hist_names = ['TT','VV','W','ZLL','EmbedZTT']
    #for name in hist_names:
    #  hist = outfile.Get(nodename+'/'+name)
    #  frac = hist.Integral(-1,-1)/tot_sub
    #  out_str+='%s = %.2f, ' % (name, frac)
    #print out_str

    
def CompareShapes(compare_w_shapes, compare_qcd_shapes):
    if compare_w_shapes:
      nominal_hist = outfile.Get(nodename+'/W')
      nominal_scale = nominal_hist.Integral(0,nominal_hist.GetNbinsX()+2)
      directory = outfile.Get(nodename)
      outfile.cd(nodename)
      shape_hist = outfile.Get(nodename+'/W_shape')
      shape_scale = shape_hist.Integral(0,shape_hist.GetNbinsX()+2)
      shape_hist.Scale(nominal_scale/shape_scale)
      shape_hist.Write()
    if compare_qcd_shapes:
      nominal_hist = outfile.Get(nodename+'/QCD')
      nominal_scale = nominal_hist.Integral()
      directory = outfile.Get(nodename)
      outfile.cd(nodename)
      shape_hist = outfile.Get(nodename+'/QCD_shape')
      shape_scale = shape_hist.Integral()
      shape_hist.Scale(nominal_scale/shape_scale)
      shape_hist.Write()
      
def AppendNameToSamples(samples=[],name_to_add=None):
  if name_to_add is None or name_to_add is '': return samples
  elif type(samples) is dict:
    new_samples = {}  
    for key in samples: new_samples[key] = samples[key]+name_to_add
    return new_samples  
  else:
    new_samples = []
    for sample in samples: new_samples.append(sample+name_to_add)    
    return new_samples
    
def RunPlotting(ana, cat='',cat_data='', sel='', add_name='', wt='wt', do_data=True, samples_to_skip=[], outfile='output.root',ff_syst_weight=None):
    doTTJ = 'TTJ' not in samples_to_skip
    doTTT = 'TTT' not in samples_to_skip
    doVVJ = 'VVJ' not in samples_to_skip
    doVVT = 'VVT' not in samples_to_skip
    doZL  = 'ZL'  not in samples_to_skip
    doZJ  = 'ZJ'  not in samples_to_skip
    
    zll_samples=list(ztt_samples)
    if options.analysis in ['cpdecay','cpprod','mssmrun2','vlq']: zll_samples+=ewkz_samples

    # produce template for observed data
    if do_data:
        if options.do_ss:
          OSSS = '!os'
        else:
            OSSS = 'os'
        weight='wt'
        if options.add_wt : weight+='*'+options.add_wt
        full_selection = BuildCutString(weight, sel, cat_data, OSSS)
        ana.nodes[nodename].AddNode(ana.SummedFactory('data_obs', data_samples, plot_unmodified, full_selection))
    
    # produce templates for backgrounds
    if options.method in [17] and options.channel != "em":
        doVVJ=False
        doTTJ=False

        if 'jetFakes' not in samples_to_skip:
            GenerateFakeTaus(ana, add_name, data_samples, plot, plot_unmodified, wt, sel, options.cat,not options.do_ss,ff_syst_weight)
            #if options.channel == 'tt':
        
        # use existing methods to calculate background due to non-fake taus
        add_fake_factor_selection = "gen_match_2<6"
        if options.channel == "tt": add_fake_factor_selection = "gen_match_1<6 && gen_match_2<6"
        residual_cat=cat+"&&"+add_fake_factor_selection
        if 'EmbedZTT' not in samples_to_skip and options.embedding:
           GenerateEmbedded(ana, add_name, embed_samples, plot, wt, sel, residual_cat, z_sels, not options.do_ss)  
            #if do_data: GenerateZTT(ana, add_name, ztt_samples+top_samples+vv_samples+ewkz_samples, plot, wt, sel, residual_cat, z_sels, not options.do_ss)
        #if 'ZTT' not in samples_to_skip and not options.embedding:
        if 'ZTT' not in samples_to_skip:
            GenerateZTT(ana, add_name, ztt_samples, plot, wt, sel, residual_cat, z_sels, not options.do_ss)
        if 'ZLL' not in samples_to_skip:
            GenerateZLL(ana, add_name, zll_samples, plot, wt, sel, residual_cat, z_sels, not options.do_ss,doZL,False)
        if 'TT' not in samples_to_skip:    
            GenerateTop(ana, add_name, top_samples, plot, wt, sel, residual_cat, top_sels, not options.do_ss, doTTT, doTTJ)  
        if 'VV' not in samples_to_skip:
            GenerateVV(ana, add_name, vv_samples, plot, wt, sel, residual_cat, vv_sels, not options.do_ss, doVVT, doVVJ)  
        #if 'EWKZ' not in samples_to_skip and options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'tauid2016','cpsummer17','tauid2017','cp18','mvadm2016']: 
            #GenerateEWKZ(ana, add_name, ewkz_samples, plot, wt, sel, residual_cat, z_sels, not options.do_ss)
        #if 'ggH_hww' not in samples_to_skip and 'qqH_hww' not in samples_to_skip and options.analysis == 'cpprod':
        #  GenerateHWW(ana, add_name, gghww_samples, qqhww_samples, plot, wt, sel, cat, not options.do_ss, True, True)
        if 'W' not in samples_to_skip and options.channel=='tt' and options.analysis in ['cpprod','cpdecay','mssmrun2','vlq'] and 'VV' not in samples_to_skip and 'ZTT' not in samples_to_skip and not options.ml_ff:
            GenerateW(ana, 'fakes'+add_name, ztt_samples+vv_samples+wjets_samples+ewkz_samples+top_samples, data_samples, wgam_samples, plot, plot_unmodified, wt, sel+'&&gen_match_1!=6&&gen_match_2==6', cat, cat_data, 8, qcd_os_ss_ratio, not options.do_ss)
        if options.channel in ['mt','et']:
          # need to add back the other fake components when testing the FF validations
          if options.w_ff_closure and not options.ml_ff:
            if 'ZLL' not in samples_to_skip:
                GenerateZLL(ana, add_name+'_res', ztt_samples, plot, wt, sel, cat, z_sels, not options.do_ss,False,True)
            if 'TT' not in samples_to_skip:
                GenerateTop(ana, add_name+'_res', top_samples, plot, wt, sel, cat, top_sels, not options.do_ss, False, True)
            if 'VV' not in samples_to_skip:
                GenerateVV(ana, add_name+'_res', vv_samples, plot, wt, sel, cat, vv_sels, not options.do_ss, False, True)
            if 'QCD' not in samples_to_skip:
                GenerateQCD(ana, add_name+'_res', data_samples, plot, plot_unmodified, wt, sel, cat, cat_data, 12, 1.1, not options.do_ss,wshift)
          if options.qcd_ff_closure and not options.ml_ff:
            if 'ZLL' not in samples_to_skip:
                GenerateZLL(ana, add_name+'_res', ztt_samples, plot, wt, sel, cat, z_sels, not options.do_ss,False,True)
            if 'TT' not in samples_to_skip:
                GenerateTop(ana, add_name+'_res', top_samples, plot, wt, sel, cat, top_sels, not options.do_ss, False, True)
            if 'VV' not in samples_to_skip:
                GenerateVV(ana, add_name+'_res', vv_samples, plot, wt, sel, cat, vv_sels, not options.do_ss, False, True)
            if 'W' not in samples_to_skip:
                GenerateW(ana, add_name+'_res', wjets_samples, data_samples, wgam_samples, plot, plot_unmodified, wt, sel, cat, cat_data, 8, qcd_os_ss_ratio, not options.do_ss)

    else:
        method = options.method
        if options.method == 18:

            if 'jetFakes' not in samples_to_skip:
                GenerateFakeTaus(ana, add_name, data_samples, plot, plot_unmodified, wt, sel, options.cat,not options.do_ss,ff_syst_weight)
            if options.channel == 'tt': method = 8
            elif options.cat == "btag_loosemt" or options.cat == "btag_tight": method = 16
            elif options.channel == 'et' or options.channel == 'mt': method = 12
        if 'EmbedZTT' not in samples_to_skip and options.embedding and options.channel != 'zmm' and options.channel != 'zee':    
            GenerateEmbedded(ana, add_name, embed_samples, plot, wt, sel, cat, z_sels, not options.do_ss)
            if options.cat in ['pass','fail']:
              cat_pass = '('+cats['pass']+')*('+cats['baseline']+')'  
              cat_fail = '('+cats['fail']+')*('+cats['baseline']+')'
              GenerateEmbedded(ana, '_pass'+add_name, embed_samples, plot, wt, sel, cat_pass, z_sels, not options.do_ss)
              GenerateEmbedded(ana, '_fail'+add_name, embed_samples, plot, wt, sel, cat_fail, z_sels, not options.do_ss)
        if 'ZTT' not in samples_to_skip and not options.embedding:
          GenerateZTT(ana, add_name, ztt_samples, plot, wt, sel, cat, z_sels, not options.do_ss)
        if 'ZTT' not in samples_to_skip and options.embedding and 'VV' not in samples_to_skip and 'TT' not in samples_to_skip:
          GenerateZTT(ana, add_name, ztt_samples+top_samples+vv_samples+ewkz_samples, plot, wt, sel, cat, z_sels, not options.do_ss)
        if 'ZLL' not in samples_to_skip:
            GenerateZLL(ana, add_name, zll_samples, plot, wt, sel, cat, z_sels, not options.do_ss,doZL,doZJ)
        if options.embedding and options.channel in ['zmm','zee'] and 'EmbedZLL' not in samples_to_skip: GenerateZLEmbedded(ana, add_name, embed_samples, plot, wt, sel, cat, z_sels, not options.do_ss)
        if 'TT' not in samples_to_skip:    
            GenerateTop(ana, add_name, top_samples, plot, wt, sel, cat, top_sels, not options.do_ss, doTTT, doTTJ) 
            if 'mvadm' in options.cat:
              cat_mvarho = '('+cats['mvadm_rho']+')*('+cats['baseline']+')'
              cat_mvaa1 = '('+cats['mvadm_a1']+')*('+cats['baseline']+')'
              cat_mvapi = '('+cats['mvadm_pi']+')*('+cats['baseline']+')'
              cat_mvanotrho = '('+cats['mvadm_notrho']+')*('+cats['baseline']+')'
              GenerateTop(ana, '_mvarho'+add_name, top_samples, plot, wt, sel, cat_mvarho, top_sels, not options.do_ss, doTTT, False)
              GenerateTop(ana, '_mvaa1'+add_name, top_samples, plot, wt, sel, cat_mvaa1, top_sels, not options.do_ss, doTTT, False)
              GenerateTop(ana, '_mvapi'+add_name, top_samples, plot, wt, sel, cat_mvapi, top_sels, not options.do_ss, doTTT, False)
              GenerateTop(ana, '_mvanotrho'+add_name, top_samples, plot, wt, sel, cat_mvanotrho, top_sels, not options.do_ss, doTTT, False)             
        if 'VV' not in samples_to_skip:
            GenerateVV(ana, add_name, vv_samples, plot, wt, sel, cat, vv_sels, not options.do_ss, doVVT, doVVJ) 
            if 'mvadm' in options.cat:
              cat_mvarho = '('+cats['mvadm_rho']+')*('+cats['baseline']+')'
              cat_mvaa1 = '('+cats['mvadm_a1']+')*('+cats['baseline']+')'
              cat_mvapi = '('+cats['mvadm_pi']+')*('+cats['baseline']+')'
              cat_mvanotrho = '('+cats['mvadm_notrho']+')*('+cats['baseline']+')'
              GenerateVV(ana, '_mvarho'+add_name, vv_samples, plot, wt, sel, cat_mvarho, vv_sels, not options.do_ss, doVVT, False)
              GenerateVV(ana, '_mvaa1'+add_name, vv_samples, plot, wt, sel, cat_mvaa1, vv_sels, not options.do_ss, doVVT, False)
              GenerateVV(ana, '_mvapi'+add_name, vv_samples, plot, wt, sel, cat_mvapi, vv_sels, not options.do_ss, doVVT, False)
              GenerateVV(ana, '_mvanotrho'+add_name, vv_samples, plot, wt, sel, cat_mvanotrho, vv_sels, not options.do_ss, doVVT, False) 
        if 'W' not in samples_to_skip:
            sel_mod=sel
            if options.method==0 and True in ['baseline_aisotau1' in x for x in options.set_alias] and options.channel=='tt': sel_mod =sel+'&&(gen_match_1==6)'
            if options.method==0 and  True in ['baseline_aisotau2' in x for x in options.set_alias] and options.channel=='tt': sel_mod =sel+'&&(gen_match_2==6)'
            GenerateW(ana, add_name, wjets_samples, data_samples, wgam_samples, plot, plot_unmodified, wt, sel_mod, cat, cat_data, method, qcd_os_ss_ratio, not options.do_ss)
        if 'QCD' not in samples_to_skip:
            GenerateQCD(ana, add_name, data_samples, plot, plot_unmodified, wt, sel, cat, cat_data, method, qcd_os_ss_ratio, not options.do_ss,wshift)
        #if 'EWKZ' not in samples_to_skip and options.era in ['smsummer16','cpsummer16','cpdecay16',"legacy16",'tauid2016','cpsummer17','tauid2017','cp18','mvadm2016'] and options.method!=0: 
        #    GenerateEWKZ(ana, add_name, ewkz_samples, plot, wt, sel, cat, z_sels, not options.do_ss) 
        #if 'ggH_hww' not in samples_to_skip and 'qqH_hww' not in samples_to_skip and options.analysis in ['cpprod','mssmrun2'] and options.channel == 'em':
        #  GenerateHWW(ana, add_name, gghww_samples, qqhww_samples, plot, wt, sel, cat, not options.do_ss, True, True)    
        if options.method==0 and options.channel=='tt':
            sel_mod = sel
            if True in ['baseline_aisotau1' in x for x in options.set_alias]: sel_mod =sel+'&&(gen_match_1!=6)'
            if True in ['baseline_aisotau2' in x for x in options.set_alias]: sel_mod =sel+'&&(gen_match_2!=6)' 
            GenerateW(ana, add_name+'_realextra', wjets_samples, data_samples, wgam_samples, plot, plot_unmodified, wt, sel_mod, cat, cat_data, method, qcd_os_ss_ratio, not options.do_ss)
        # used to plot em QCD uncerts
        #GenerateQCD(ana, '_shapeup', data_samples, plot, plot_unmodified, wt+'*wt_em_qcd_shapeup', sel, cat, cat_data, method, qcd_os_ss_ratio, not options.do_ss,wshift)
        #GenerateQCD(ana, '_shapedown', data_samples, plot, plot_unmodified, wt+'*wt_em_qcd_shapedown', sel, cat, cat_data, method, qcd_os_ss_ratio, not options.do_ss,wshift)
  
        if compare_w_shapes:
          cat_relax=cats['w_shape_comp']
          GenerateW(ana, '_shape', wjets_samples, data_samples, wgam_samples, plot, plot_unmodified, wt, sel, cat_relax, cats_unmodified['w_shape_comp'], 8, qcd_os_ss_ratio, not options.do_ss)    
        if compare_qcd_shapes:
          if options.channel == 'em':
            GenerateQCD(ana, '_shape', data_samples, plot, plot_unmodified, wt, sel, cat, cat_data, 19, qcd_os_ss_ratio, not options.do_ss,wshift)
          else:
            cat_relax=cats['qcd_shape_comp']
            GenerateQCD(ana, '_shape', data_samples, plot, plot_unmodified, wt, sel, cat_relax,cats_unmodified['qcd_shape_comp'], method, qcd_os_ss_ratio, not options.do_ss)
           
    if 'signal' not in samples_to_skip:
        if options.analysis in ['sm','cpprod','cpdecay']:
            procs=[]
            for proc in sm_samples:
                if True not in [samp in proc for samp in samples_to_skip]: procs.append(proc)   
            if options.analysis == 'cpdecay' and options.sm_masses!="": GenerateReweightedCPSignal(ana, add_name, plot, wt, sel, cat, not options.do_ss) 
            elif options.analysis == 'cpprod' and options.sm_masses!="": GenerateReweightedCPProdSignal(ana, add_name, plot, wt, sel, cat, not options.do_ss) 
            else: GenerateSMSignal(ana, add_name, plot, sm_masses, wt, sel, cat, not options.do_ss,processes=procs)
        elif options.analysis in ['mssm','mssmrun2'] and (options.ggh_masses != "" or options.bbh_masses != "" or options.ggh_masses_powheg != "" or options.bbh_masses_powheg != "") :
            bbh_add_name = ''
            if options.bbh_nlo_masses and not options.bbh_masses_powheg: bbh_add_name = '-LO'
            #GenerateMSSMSignal(ana, add_name, bbh_add_name, plot, ggh_masses, bbh_masses, wt, sel, cat, not options.do_ss)
            GenerateMSSMSignal(ana, add_name, bbh_add_name, plot, [], bbh_masses, wt, sel, cat, not options.do_ss)
            if options.add_sm_background:
                GenerateSMSignal(ana, add_name, plot, ['125'],  wt, sel, cat, not options.do_ss, options.add_sm_background)  
        elif options.analysis == 'Hhh':
            GenerateHhhSignal(ana, add_name, plot, ggh_masses, wt, sel, cat, not options.do_ss)
        elif options.analysis == 'vlq':
            GenerateVLQSignal(ana, add_name, plot, options.vlq_sig, wt, sel, cat, not options.do_ss)
        if options.analysis in ['mssm','mssmrun2'] and options.bbh_nlo_masses != "":
            GenerateNLOMSSMSignal(ana, add_name, plot, [''], bbh_nlo_masses, wt, sel, cat, options.doNLOScales, options.doPDF, not options.do_ss)
        if options.analysis in ['mssm','mssmrun2'] and options.doMSSMReWeighting:
          GenerateReWeightedMSSMSignal(ana, add_name, plot, ggh_masses, wt, sel, cat, not options.do_ss) 
          
    if options.syst_embedding_tt and options.embedding and systematic == 'default':
        GenerateTop(ana, '_embed_syst', top_samples, plot, wt, sel, cat, top_sels_embed, not options.do_ss, True, False)        
        GenerateVV(ana, '_embed_syst', vv_samples, plot, wt, sel, cat, top_sels_embed, not options.do_ss, True, False)

def Get1DBinNumFrom2D(h2d,xbin,ybin):
    Nxbins = h2d.GetNbinsX()
    return (ybin-1)*Nxbins + xbin -1

def Get1DBinNumFrom3D(h3d,xbin,ybin,zbin):
    Nxbins = h3d.GetNbinsX()
    Nybins = h3d.GetNbinsY()
    return (zbin-1)*Nxbins*Nybins + (ybin-1)*Nxbins + xbin -1

def UnrollHist2D(h2d,inc_y_of=True):
    # inc_y_of = True includes the y over-flow bins
    if inc_y_of: n = 1
    else: n = 0
    Nbins = (h2d.GetNbinsY()+n)*(h2d.GetNbinsX())
    h1d = ROOT.TH1D(h2d.GetName(), '', Nbins, 0, Nbins)
    for i in range(1,h2d.GetNbinsX()+1):
      for j in range(1,h2d.GetNbinsY()+1+n):
        glob_bin = Get1DBinNumFrom2D(h2d,i,j)
        content = h2d.GetBinContent(i,j)
        error = h2d.GetBinError(i,j)
        h1d.SetBinContent(glob_bin+1,content)
        h1d.SetBinError(glob_bin+1,error)
        #if 'sjdphi' in plot: h1d.GetXaxis().SetBinLabel(glob_bin+1,'%.1f-%.1f' % (h2d.GetXaxis().GetBinLowEdge(i),h2d.GetXaxis().GetBinLowEdge(i+1)))
        #else:
        #  h1d.GetXaxis().SetBinLabel(glob_bin+1,'%.0f-%.0f' % (h2d.GetXaxis().GetBinLowEdge(i),h2d.GetXaxis().GetBinLowEdge(i+1)))
        #if 'sdphi' in options.var: h1d.GetXaxis().SetBinLabel(glob_bin+1,'%.1f-%.1f' % (h2d.GetXaxis().GetBinLowEdge(i),h2d.GetXaxis().GetBinLowEdge(i+1)))
    #h1d.LabelsOption('v','X')
    return h1d

def UnrollHist3D(h3d,inc_y_of=False,inc_z_of=True):
    if inc_y_of: ny = 1
    else: ny = 0
    if inc_z_of: nz = 1
    else: nz = 0
    
    Nbins = (h3d.GetNbinsZ()+nz)*(h3d.GetNbinsY()+ny)*(h3d.GetNbinsX())
    h1d = ROOT.TH1D(h3d.GetName(), '', Nbins, 0, Nbins)
    for i in range(1,h3d.GetNbinsX()+1):
      for j in range(1,h3d.GetNbinsY()+1+ny):
        for k in range(1,h3d.GetNbinsZ()+1+nz):    
          glob_bin = Get1DBinNumFrom3D(h3d,i,j,k)
          content = h3d.GetBinContent(i,j,k)
          error = h3d.GetBinError(i,j,k)
          h1d.SetBinContent(glob_bin+1,content)
          h1d.SetBinError(glob_bin+1,error)
          #if 'sjdphi' in plot: h1d.GetXaxis().SetBinLabel(glob_bin+1,'%.1f-%.1f' % (h3d.GetXaxis().GetBinLowEdge(i),h3d.GetXaxis().GetBinLowEdge(i+1)))
          #else:
          #  h1d.GetXaxis().SetBinLabel(glob_bin+1,'%.0f-%.0f' % (h3d.GetXaxis().GetBinLowEdge(i),h3d.GetXaxis().GetBinLowEdge(i+1)))
         # if 'sdphi' in options.var: h1d.GetXaxis().SetBinLabel(glob_bin+1,'%.1f-%.1f' % (h3d.GetXaxis().GetBinLowEdge(i),h3d.GetXaxis().GetBinLowEdge(i+1)))
    #h1d.LabelsOption('v','X')
    return h1d

def NormSignals(outfile,add_name):
    # When adding signal samples to the data-card we want to scale all XS to 1pb - correct XS times BR is then applied at combine harvestor level
    samples_to_skip = add_names_dict[add_name] 
    if 'signal' not in samples_to_skip:
        outfile.cd(nodename)
        if options.analysis in ['sm','cpprod','cpdecay'] or options.add_sm_background:
            if options.analysis == "sm":
                masses = sm_masses
            else:
                masses = [options.add_sm_background]
            for samp in sm_samples:
                if options.analysis == "sm":
                    samp_name = samp
                else:
                    samp_name = samp+"_SM"
                if masses is not None:    
                    for mass in masses:        
                        if isinstance(sm_samples[samp], (list,)): xs = ana.info[sm_samples[samp][0].replace('*',mass)]['xs']
                        else: xs = ana.info[sm_samples[samp].replace('*',mass)]['xs']
                        if xs == 1.: continue
                        sf = 1.0/xs
                        if not outfile.GetDirectory(nodename).GetListOfKeys().Contains(samp_name+mass+add_name): continue
                        sm_hist = outfile.Get(nodename+'/'+samp_name+mass+add_name)
                        sm_hist.Scale(sf)
                        sm_hist.Write("",ROOT.TObject.kOverwrite)
        if options.analysis in ["mssm",'mssmrun2']:
            for samp in mssm_samples:
                if samp == 'ggH':
                    masses = ggh_masses
                elif samp == 'bbH' and not options.bbh_nlo_masses:
                    masses = bbh_masses
                elif 'bbH' in samp:
                    masses = bbh_nlo_masses
                if masses is not None:    
                    for mass in masses:
                        xs = ana.info[mssm_samples[samp].replace('*',mass)]['xs']
                        if xs == 1.: continue
                        sf = 1.0/xs
                        if sf == 1.0: continue
                        mssm_hist = outfile.Get(nodename+'/'+samp+mass+add_name)
                        mssm_hist.Scale(sf)
                        mssm_hist.Write("",ROOT.TObject.kOverwrite)
                        if options.doMSSMReWeighting and samp == 'ggH': 
                          re_weighted_names = ['ggh_t_','ggh_b_','ggh_i_','ggH_t_','ggH_b_','ggH_i_','ggA_t_','ggA_b_','ggA_i_']
                          #re_weighted_names = ['ggH_t_','ggH_b_','ggH_i_','ggA_t_','ggA_b_','ggA_i_']
                          for name in re_weighted_names:
                            mssm_hist = outfile.Get(nodename+'/'+name+mass+add_name)
                            #mssm_hist = ana.nodes[nodename].nodes[name+mass+add_name].shape.hist
                            mssm_hist.Scale(sf)
                            mssm_hist.Write("",ROOT.TObject.kOverwrite)
        if options.analysis == "Hhh":
            for samp in Hhh_samples:
                masses = ggh_masses
                if masses is not None:
                    for mass in masses:
                        xs = ana.info[Hhh_samples[samp].replace('*',mass)]['xs']
                        if xs == 1.: continue
                        sf = 1.0/xs
                        mssm_hist = outfile.Get(nodename+'/'+samp+mass+add_name)
                        mssm_hist.Scale(sf)
                        mssm_hist.Write("",ROOT.TObject.kOverwrite)
        if options.analysis == "vlq":
            for samp in vlq_samples:
              xs = ana.info[vlq_samples[samp]]['xs']
              if xs == 1.: continue
              sf = 1.0/xs
              vlq_hist = outfile.Get(nodename+'/'+samp+add_name)
              vlq_hist.Scale(sf)
              vlq_hist.Write("",ROOT.TObject.kOverwrite)
        outfile.cd()

def RenameMSSMrun2Datacards(outfile):
  chan = options.channel
  renames = {
   'CMS_eff_trigger_' : 'CMS_eff_trigger_emb_',
   'CMS_eff_xtrigger_l_' : 'CMS_eff_xtrigger_l_emb_',
   'CMS_eff_xtrigger_t_': 'CMS_eff_xtrigger_t_emb_',
   'CMS_eff_t_': 'CMS_eff_t_emb_',
   'CMS_scale_e' : 'CMS_scale_e_emb',
   'CMS_scale_t_' : 'CMS_scale_t_emb_',
   'CMS_eff_trigger_single_t_':'CMS_eff_trigger_single_t_emb_',
   'scale_embed_met_' : 'scale_embed_met_%s_' % chan,
  }  
  directory = outfile.Get(nodename)

  # count number of directories
  count = 0    
  for key in directory.GetListOfKeys(): count += 1

  i = 0
  for key in directory.GetListOfKeys():
    if i < count:
      name = key.GetName()
      histo = directory.Get(name)
      if not isinstance(histo,ROOT.TDirectory) and 'EmbedZTT' in name:
        new_name = name.replace('EmbedZTT', 'EMB')
        histo.SetName(new_name)
        directory.cd()
        #histo.Write(new_name,ROOT.TObject.kWriteDelete)
        histo.Write(new_name)
        directory.Delete(name+';1')
        for x in renames:
          if x in new_name:
            histo_clone = histo.Clone()
            y = renames[x]
            new_name_2 = new_name.replace(x,y)
            print new_name,new_name_2
            histo_clone.SetName(new_name_2)
            histo_clone.Write(new_name_2)
            break
      elif not isinstance(histo,ROOT.TDirectory) and 'VVT' in name and not 'VVT_for_ZTT' in name:
        new_name = name.replace('VVT', 'VVL')
        histo.SetName(new_name)
        directory.cd()
        histo.Write(new_name)
        directory.Delete(name+';1')
      elif not isinstance(histo,ROOT.TDirectory) and 'TTT' in name and not 'TTT_for_ZTT' in name:
        new_name = name.replace('TTT', 'TTL')
        histo.SetName(new_name)
        directory.cd()
        histo.Write(new_name)
        directory.Delete(name+';1')
      elif not isinstance(histo,ROOT.TDirectory) and 'for_ZTT' in name:
        new_name = name.replace('_for_ZTT', '')
        histo.SetName(new_name)
        directory.cd()
        histo.Write(new_name)
        directory.Delete(name+';1')
      elif not isinstance(histo,ROOT.TDirectory) and 'Wfakes' in name:
        new_name = name.replace('Wfakes', 'wFakes')
        histo.SetName(new_name)
        directory.cd()
        histo.Write(new_name)
        directory.Delete(name+';1')
      #elif not isinstance(histo,ROOT.TDirectory) and ('ggH_' in name or 'ggA_' in name) and name[:6].count('_') == 1:
        #new_name = name[:5]+'_'+name[5:]
        #histo.SetName(new_name)
        #directory.cd()
        #histo.Write(new_name)
        #directory.Delete(name+';1')
      elif not isinstance(histo,ROOT.TDirectory) and 'bbH' in name and name[:4].count('_') == 0:
        new_name = name[:3]+'_'+name[3:]
        histo.SetName(new_name)
        directory.cd()
        histo.Write(new_name)
        directory.Delete(name+';1')
      elif isinstance(histo,ROOT.TDirectory):
        directory.Delete(name+';1')
      #elif not isinstance(histo,ROOT.TDirectory) and ('WplusH' in name or 'WminusH' in name):
        #directory.Delete(name+';1')
      elif not isinstance(histo,ROOT.TDirectory) and 'ggH' in name and name[:4].count('_') == 0 and 'ggH125_SM' not in name and 'ggH95' not in name:
        directory.Delete(name+';1')
      elif not isinstance(histo,ROOT.TDirectory) and 'ggH125_SM' in name :
        new_name = name.replace('ggH125_SM', 'ggH125')
        histo.SetName(new_name)
        directory.cd()
        histo.Write(new_name)
        directory.Delete(name+';1')

    i += 1



def TotalUnc(h0, hists=[]):
  #sum in quadrature several systematic uncertainties to form total uncertainty band
  hout = h0.Clone()
  hup = h0.Clone()
  hdown = h0.Clone()
  hout.SetName(h0.GetName()+'_uncerts_total')
  hup.SetName(h0.GetName()+'_uncerts_total_up')
  hdown.SetName(h0.GetName()+'_uncerts_total_down')
  for i in range(1,h0.GetNbinsX()+2):
    x0 = h0.GetBinContent(i)
    uncerts_up = [0.]
    uncerts_down = [0.]
    for h in hists:
      x = h.GetBinContent(i)
      if x>x0: uncerts_up.append(x-x0)
      if x<x0: uncerts_down.append(x0-x)
    up = 0.
    down = 0.
    for u in uncerts_up: up+=u**2
    for u in uncerts_down: down+=u**2
    up = up**.5
    down = down**.5

    hup.SetBinContent(i,x0+up)
    hdown.SetBinContent(i,x0-down)
    c = (x0+up + x0-down)/2
    u = (up+down)/2
    hout.SetBinContent(i,c)
    hout.SetBinError(i,u)
  return (hout, hup, hdown)

def RawFFFromString(string):
  string = string[1:-1]
  bracket_count = 0
  new_string = ''
  for i in string:
    if i == '(': bracket_count += 1
    elif i == ')': bracket_count -= 1
    new_string += i
    if bracket_count == 0: break
  return new_string

# Create output file
is_2d=False
is_3d=False
var_name = options.var.split('[')[0]
var_name = var_name.split('(')[0]
if var_name.count(',') == 1:
    is_2d = True
    var_name = var_name.split(',')[0]+'_vs_'+var_name.split(',')[1]
if var_name.count(',') == 2:
    is_3d = True
    var_name = var_name.split(',')[0]+'_vs_'+var_name.split(',')[1]+'_vs_'+var_name.split(',')[2]    

if options.datacard != "": datacard_name = options.datacard
else: datacard_name = options.cat
if options.extra_name != "": 
#  output_name = options.outputfolder+'/datacard_'+options.extra_name+'_'+datacard_name+'_'+options.channel+'_'+options.year+'.root'
  datacard_name+='_'+options.extra_name
#else: 
output_name = options.outputfolder+'/datacard_'+var_name+'_'+datacard_name+'_'+options.channel+'_'+options.year+'.root'
outfile = ROOT.TFile(output_name, 'RECREATE')
    
cats['cat'] = '('+cats[options.cat]+')*('+cats['baseline']+')'
if options.channel=="em": cats['em_shape_cat'] = '('+cats[options.cat]+')*('+cats['loose_baseline']+')'
sel = options.sel
plot = options.var
plot_unmodified = plot
if options.datacard != "": nodename = options.channel+'_'+options.datacard
else: nodename = options.channel+'_'+options.cat   

add_names_dict = {}
add_names = []
cats_unmodified = copy.deepcopy(cats)

syst_names = {}

max_systs_per_pass = 30 # code uses too much memory if we try and process too many systematics at once so set the maximum number of systematics processed per loop here
while len(systematics) > 0:
  ana = Analysis()
  #if options.syst_scale_j_by_source != '': ana.writeSubnodes(False) # storing subnodes uses too much memory when doing JES uncertainties split by source
  
  ana.remaps = {}
  if options.channel == 'em':
      ana.remaps['MuonEG'] = 'data_obs'
  elif options.channel in ['mt','mj','zmm']:
      ana.remaps['SingleMuon'] = 'data_obs'
  elif (options.era != 'cp18' or options.era != 'UL_18') and (options.channel == 'et' or options.channel == 'zee'):
      ana.remaps['SingleElectron'] = 'data_obs'
  elif (options.era == 'cp18' or options.era == 'UL_18') and (options.channel == 'et' or options.channel == 'zee'):
      ana.remaps['EGamma'] = 'data_obs'
  elif options.channel == 'tt':
      ana.remaps['Tau'] = 'data_obs'  
      
  ana.nodes.AddNode(ListNode(nodename))

  prev_dir=None    
  for index, systematic in enumerate(list(systematics.keys())[:max_systs_per_pass]):
      syst_names[systematic] = systematics[systematic][1]
      if prev_dir is not None and systematics[systematic][0] is not prev_dir: continue # this ensures that we process the same trees from every call to ana.Run() - i.e trees in sub-directory systematics[systematic][0]
      prev_dir = systematics[systematic][0]
      print "Processing:", systematic
      print ""
      print systematics[systematic][0]

      plot = options.var
      cats=copy.deepcopy(cats_unmodified)
      wshift=1.0
      if systematic == 'syst_qcd_shape_wsf_up' and w_abs_shift is not None: wshift+=w_abs_shift
      if systematic == 'syst_qcd_shape_wsf_down' and w_abs_shift is not None: wshift-=w_abs_shift
      if options.syst_scale_j_by_source != '' and 'syst_scale_j_by_source' in systematic:
        # if JES systematic split by source then the category and plotting variable strings need to be modified to use the shifted variables
        replace_dict = systematics[systematic][5]  
        for cat in cats: cats[cat] = OverwriteNames(cats[cat], replace_dict)
        plot = OverwriteNames(plot, replace_dict)
      
      add_folder_name = systematics[systematic][0]
      add_name = systematics[systematic][1]
      isFFSyst = systematics[systematic][4]
      ff_syst_weight = None
      if not isFFSyst: weight = systematics[systematic][2]
      else:
          weight='wt'
          ff_syst_weight = systematics[systematic][2]
      if options.add_wt is not "": weight+="*"+options.add_wt
      if options.channel == "tt" and options.era == 'mssmsummer16': weight+='*wt_tau_id_medium'
      if options.channel == "tt" and options.era in ['smsummer16']: weight+='*wt_tau_id_tight'
      if options.cat == '0jet' and options.era in ['smsummer16']: weight+='*wt_lfake_rate'

      samples_to_skip = systematics[systematic][3]
      add_names.append(add_name)
      syst_add_name=add_folder_name
      

      mc_input_folder_name = options.folder
      if add_folder_name != '': mc_input_folder_name += '/'+add_folder_name
      
      if options.signal_folder: signal_mc_input_folder_name = options.signal_folder
      else: signal_mc_input_folder_name = options.folder
      if add_folder_name != '': signal_mc_input_folder_name += '/'+add_folder_name
      
      if options.embed_folder: embed_input_folder_name = options.embed_folder
      else: embed_input_folder_name = options.folder
      if add_folder_name != '' and 'EmbedZTT' not in samples_to_skip: embed_input_folder_name += '/'+add_folder_name
    
      # Add all data files
      for sample_name in data_samples:
          ana.AddSamples(options.folder+'/'+sample_name+'_'+options.channel+'_{}.root'.format(options.year), 'ntuple', None, sample_name)
      
      # Add all MC background files
      for sample_name in ztt_samples + vv_samples + wgam_samples + top_samples + ztt_shape_samples + wjets_samples+ewkz_samples+gghww_samples+qqhww_samples:
          ana.AddSamples(mc_input_folder_name+'/'+sample_name+'_'+options.channel+'_{}.root'.format(options.year), 'ntuple', None, sample_name)
          #ana.AddSamples(mc_input_folder_name+'/'+sample_name+'_'+options.channel+'_{}.root'.format(options.year), 'ntuple', options.folder+'/'+sample_name+'_'+options.channel+'_{}.root'.format(options.year), sample_name) # this fixes issues if a sample is not included in systematic sub directory (e.g because systematics doesn't affect it) but at the same time can make it easier to miss issues like a sample missing that should be there 
 
      # Add embedded samples if using
      if options.embedding: 
        for sample_name in embed_samples:
          ana.AddSamples(embed_input_folder_name+'/'+sample_name+'_'+options.channel+'_{}.root'.format(options.year), 'ntuple', None, sample_name)
          
      # Add all MC signal files
      
      if options.analysis in ['sm','cpprod','cpdecay']:
          signal_samples = sm_samples
      elif options.analysis in ['mssm','mssmrun2']:
          signal_samples = mssm_samples
          if options.bbh_nlo_masses and mssm_nlo_samples: signal_samples['bbH'] = mssm_nlo_samples['bbH']
          if options.nlo_qsh: signal_samples.update(mssm_nlo_qsh_samples)
          if options.bbh_nlo_masses and options.bbh_masses:  signal_samples.update(mssm_lo_samples)
      elif options.analysis == 'Hhh':
          signal_samples = Hhh_samples
      elif options.analysis == "vlq":
          signal_samples = vlq_samples

      for samp in signal_samples:
          if options.analysis in ['sm','cpprod','cpdecay']:
              masses=sm_masses
          elif samp == 'ggH':
              masses = ggh_masses
          elif (samp == 'bbH' and not options.bbh_nlo_masses) or samp == 'bbH-LO':
              masses = bbh_masses
          elif 'bbH' in samp:
              masses = bbh_nlo_masses
          elif options.analysis == "vlq":
              masses = None
          if masses is not None:    
              for mass in masses:
                  sample_names=[]
                  if isinstance(signal_samples[samp], (list,)): 
                     for i in signal_samples[samp]: sample_names.append(i.replace('*',mass)) 
                  else: sample_names = [signal_samples[samp].replace('*',mass)]
                  tree_name = 'ntuple'
                  if options.gen_signal: tree_name = 'gen_ntuple'
                  for sample_name in sample_names: 
                    #if 'amcatnloFXFX' in sample_name and False:
                    #  new_sig_folder = '/vols/cms/dw515/Offline/output/SM/Oct26_2016_newsig/'
                    #  if add_folder_name != '': new_sig_folder += '/'+add_folder_name
                    ana.AddSamples(signal_mc_input_folder_name+'/'+sample_name+'_'+options.channel+'_{}.root'.format(options.year), tree_name, None, sample_name)
          elif options.analysis == "vlq":
              if options.vlq_sig != "":
                if samp in options.vlq_sig.split(","):
                  sample_names=[]
                  if isinstance(vlq_samples[samp], (list,)):
                    for i in vlq_samples[samp]: sample_names.append(i)
                  else: sample_names = [vlq_samples[samp]]
                  for sample_name in sample_names:
                    ana.AddSamples(signal_mc_input_folder_name+'/'+sample_name+'_'+options.channel+'_{}.root'.format(options.year), 'ntuple', None, sample_name)
                   #ana.AddSamples(signal_mc_input_folder_name+'/'+signal_samples[samp]+'_'+options.channel+'_{}.root'.format(options.year), 'ntuple', None, signal_samples[samp])

      if options.add_sm_background and options.analysis in ['mssm','mssmrun2']:
          for samp in sm_samples:
            sample_names=[]
            if isinstance(sm_samples[samp], (list,)):
              for i in sm_samples[samp]: sample_names.append(i.replace('*',options.add_sm_background))
            else: sample_names = [sm_samples[samp].replace('*',options.add_sm_background)]
            for sample_name in sample_names:
              ana.AddSamples(signal_mc_input_folder_name+'/'+sample_name+'_'+options.channel+'_{}.root'.format(options.year), 'ntuple', None, sample_name)



              
      ana.AddInfo(options.paramfile, scaleTo='data_obs')
  
      # Add data only for default
      if systematic == 'default': do_data = True
      else: do_data = False
              
      #Run default plot 
      if options.scheme == 'signal': 
          samples_to_skip.extend(['TTT','TTJ','VVT','VVJ','W','QCD','jetFakes','ZLL','ZTT','ZL','EWKZ','ggH_hww'])
          do_data = False
      if options.scheme == "noTT":
          samples_to_skip.extend(["TTT","TTJ"])
      RunPlotting(ana, cats['cat'], cats_unmodified['cat'], sel, add_name, weight, do_data, samples_to_skip,outfile,ff_syst_weight)
      #if options.era == "tauid2016" and options.channel in ['et','mt']: 
      #    RunPlotting(ana, cats['pass']+'&&'+cats['baseline'], cats_unmodified['pass']+'&&'+cats_unmodified['baseline'], sel, "pass"+add_name, weight, False, samples_to_skip,outfile,ff_syst_weight)
      #    RunPlotting(ana, cats['fail']+'&&'+cats['baseline'], cats_unmodified['fail']+'&&'+cats_unmodified['baseline'], sel, "fail"+add_name, weight, False, samples_to_skip,outfile,ff_syst_weight)
     
      if options.do_custom_uncerts and options.custom_uncerts_wt_up != "" and options.custom_uncerts_wt_down !="":
          add_names.append("_custom_uncerts_up")
          add_names.append("_custom_uncerts_down")
          RunPlotting(ana, cats['cat'], cats_unmodified['cat'], sel, '_custom_uncerts_up', weight+'*'+options.custom_uncerts_wt_up, do_data, ['signal'],outfile,ff_syst_weight)
          RunPlotting(ana, cats['cat'], cats_unmodified['cat'], sel, '_custom_uncerts_down', weight+'*'+options.custom_uncerts_wt_down, do_data, ['signal'],outfile,ff_syst_weight)

      add_names_dict[add_name] = samples_to_skip
      
      del systematics[systematic]
  ana.Run()
  ana.nodes.Output(outfile)

  # fix negative bns,empty histograms etc.
  FixBins(ana,outfile)
  for n in add_names: 
    GetTotals(ana,n,outfile)
  PrintSummary(nodename, ['data_obs'], add_names)

if compare_w_shapes or compare_qcd_shapes: CompareShapes(compare_w_shapes, compare_qcd_shapes)
    
#if options.method in [17,18] and options.do_ff_systs: NormFFSysts(ana,outfile)
if (options.era in ["smsummer16"] and options.syst_w_fake_rate and options.method != 8) or options.era in ["tauid2016"]: NormWFakeSysts(ana,outfile)
#NormEmbedToMC(ana,outfile) # this is to check embedding sensitivity after scaling embedding to MC yields

if options.syst_embedding_tt and options.embedding and not options.no_default: TTBarEmbeddingSyst(ana,outfile,options.syst_embedding_tt)

if options.doNLOScales: 
    ScaleUncertBand(nodename,outfile)
    DONLOUncerts(nodename,outfile)
if options.doPDF:
    PDFUncerts(nodename,outfile)

def Symmetrise(hist):
  nbins = hist.GetNbinsX()
  if nbins % 2:
    print 'N X bins in 2D histogram is not even so cannot symmetrise!'
    return
  to_skip = ['data_obs','ggH','qqH','WH','WplusH','WminusH','ZH']
  to_skip = ['data_obs']
  if True in [x in hist.GetName() for x in to_skip]: return
  for i in range(1,nbins/2+1):
    lo_bin = i
    hi_bin = nbins-i+1
    for j in range(1,hist.GetNbinsY()+2):
      c1 = hist.GetBinContent(lo_bin,j)
      c2 = hist.GetBinContent(hi_bin,j)
      e1 = hist.GetBinError(lo_bin,j)
      e2 = hist.GetBinError(hi_bin,j)
      cnew = (c1+c2)/2
      enew = math.sqrt(e1**2 + e2**2)/2
      hist.SetBinContent(lo_bin,j,cnew)
      hist.SetBinContent(hi_bin,j,cnew)
      hist.SetBinError(lo_bin,j,enew)
      hist.SetBinError(hi_bin,j,enew)

def MergeXBins(hist):
  nxbins = hist.GetNbinsX()
  nybins = hist.GetNbinsY()
  nbins = hist.GetNbinsX()*hist.GetNbinsY()
  to_skip = ['data_obs','ggH','qqH','WH','WplusH','WminusH','ZH']
  if True in [x in hist.GetName() for x in to_skip]: return
  for i in range(1,nybins+2):
    tot_err = ROOT.Double()
    tot = hist.IntegralAndError(1,nxbins,i,i,tot_err)
    for j in range(1, nxbins+1):
      hist.SetBinContent(j,i,tot/nxbins)
      hist.SetBinError(j,i,tot_err/nxbins)

## print average weights
#x_lines = []
#y_labels = []
#directory = outfile.Get(nodename)
#outfile.cd(nodename)
#for key in directory.GetListOfKeys():
#  hist_name = key.GetName()
#  hist = directory.Get(hist_name).Clone()
#
#  if not isinstance(hist,ROOT.TDirectory):
#    print hist.GetName(), hist.GetEntries(), hist.Integral(-1,-1)/hist.GetEntries()

# sm 2D unrolling
if is_2d and options.do_unrolling:
  x_lines = []
  y_labels = []
  first_hist = True
  # loop over all TH2Ds and for each one unroll to produce TH1D and add to datacard
  directory = outfile.Get(nodename)  
  outfile.cd(nodename)
  hists_to_add = []
  for key in directory.GetListOfKeys():
    hist_name = key.GetName()
    hist = directory.Get(hist_name).Clone()
 
    if not isinstance(hist,ROOT.TDirectory):
      include_of = True
      if 'dijet' in options.cat: include_of = False
 
      if options.symmetrise: Symmetrise(hist)
      if options.mergeXbins: MergeXBins(hist)

      h1d = UnrollHist2D(hist,include_of)
      hists_to_add.append(h1d)
      if first_hist:
        first_hist=False
        Nxbins = hist.GetNbinsX()
        for i in range(1,hist.GetNbinsY()+1): x_lines.append(Nxbins*i)
        for j in range(1,hist.GetNbinsY()+1): y_labels.append([hist.GetYaxis().GetBinLowEdge(j),hist.GetYaxis().GetBinLowEdge(j+1)])
        if include_of: y_labels.append([hist.GetYaxis().GetBinLowEdge(hist.GetNbinsY()+1),-1])
  for hist in hists_to_add: hist.Write("",ROOT.TObject.kOverwrite)

# sm 3D unrolling
if is_3d and options.do_unrolling:
  x_lines = []
  y_labels = []
  z_labels = []
  first_hist = True
  # loop over all TH3Ds and for each one unroll to produce TH1D and add to datacard
  directory = outfile.Get(nodename)
  outfile.cd(nodename)
  hists_to_add = []
  for key in directory.GetListOfKeys():
    hist_name = key.GetName()
    hist = directory.Get(hist_name).Clone()
    if not isinstance(hist,ROOT.TDirectory):
      include_y_of = False
      include_z_of = True
      h1d = UnrollHist3D(hist,include_y_of,include_z_of)
      hists_to_add.append(h1d)
      if first_hist:
        first_hist=False
        Nxbins = hist.GetNbinsX()
        for i in range(1,hist.GetNbinsY()+1): x_lines.append(Nxbins*i)
        for j in range(1,hist.GetNbinsY()+1): y_labels.append([hist.GetYaxis().GetBinLowEdge(j),hist.GetYaxis().GetBinLowEdge(j+1)])
        if include_y_of: y_labels.append([hist.GetYaxis().GetBinLowEdge(hist.GetNbinsY()+1),-1])
        for j in range(1,hist.GetNbinsZ()+1): z_labels.append([hist.GetZaxis().GetBinLowEdge(j),hist.GetZaxis().GetBinLowEdge(j+1)])
        if include_z_of: z_labels.append([hist.GetZaxis().GetBinLowEdge(hist.GetNbinsZ()+1),-1])
  for hist in hists_to_add: hist.Write("",ROOT.TObject.kOverwrite)


# make systematic uncertainty histograms

custom_uncerts_up_name = options.custom_uncerts_up_name
custom_uncerts_down_name = options.custom_uncerts_down_name
#if len(syst_names)>1 and options.do_custom_uncerts:
if options.do_custom_uncerts:

  custom_uncerts_up_name = 'total_bkg_uncerts_total_up' 
  custom_uncerts_down_name = 'total_bkg_uncerts_total_down' 

  directory = outfile.Get(nodename)
  h0 = directory.Get('total_bkg')
  hists=[]

  for x in syst_names:
    if x == 'default': continue
    h=h0.Clone()
    syst = syst_names[x]
    print 'add syst', x, syst, h0.Integral()
    h.SetName(h0.GetName()+syst)
    for key in directory.GetListOfKeys():
      name = key.GetName()
      if name.startswith('ggH') or name.startswith('bbH') or name.startswith('qqH') or name.startswith('VH') or name.startswith('ZH') or name.startswith('WH') or name.startswith('WplusH') or name.startswith('WminusH'): continue
      histo = directory.Get(name)
      if not isinstance(histo,ROOT.TDirectory) and name.endswith(syst) and not name.startswith('jetFakes_norm'):
        #print name.replace(syst,''), name, histo.Integral() 
        histo_nom = directory.Get(name.replace(syst,''))
        histo.Add(histo_nom,-1)
        if 'scale_embed_met' in name:
          # this part scales the embed MET uncertainties to their actual values
          scales = {
            "mt_2018": 0.29,
            "mt_2017": 0.69,
            "mt_2016": 0.44,
            "et_2018": 0.29,
            "et_2017": 0.69,
            "et_2016": 0.44,
            "tt_2018": 0.29,
            "tt_2017": 0.69,
            "tt_2016": 0.44,
          }
          scale = (scales['%s_%s' % (options.channel,options.year)]**2+0.1**2)**.5
          print 'scaling ', name, ' by ', scale
          for i in range(1,histo.GetNbinsX()+1): 
            histo.SetBinContent(i,histo.GetBinContent(i)*scale)
        h.Add(histo)
    hists.append(h.Clone())

  norm_systs = {}
  if options.embedding: 
    if options.channel in ['em']:
      norm_systs['embed_yield'] = (0.04,['EmbedZTT'])
      norm_systs['dy_xs'] = (0.02,['ZLL'])
      norm_systs['lumi'] = (0.025,['TTT','VVT','ZLL'])
    else:
      norm_systs['embed_yield'] = (0.04,['EmbedZTT'])
      norm_systs['dy_xs'] = (0.02,['ZL'])
      norm_systs['lumi'] = (0.025,['TTT','VVT','ZL'])
  else: 
    if options.channel in ['em']:
      norm_systs['dy_xs'] = (0.02,['ZLL','ZTT'])
      norm_systs['lumi'] = (0.025,['TTT','VVT','ZLL','ZTT'])
    else:
      norm_systs['dy_xs'] = (0.02,['ZL','ZTT'])
      norm_systs['lumi'] = (0.025,['TTT','VVT','ZL','ZTT'])
  norm_systs['ttbar_xs'] = (0.04,['TTT'])
  norm_systs['vv_xs'] = (0.05,['VVT'])
  if options.channel == 'tt':
    if options.embedding: norm_systs['tau_id'] = (0.06,['EmbedZTT'])
    else: norm_systs['tau_id'] = (0.06,['ZTT'])
    norm_systs['wfakes_yield'] = (0.2,['Wfakes'])
    norm_systs['l_fakerate'] = (0.05,['VVT','ZL','TTT'])
  if options.channel in ['et','mt']:
    if options.embedding: 
      norm_systs['tau_id'] = (0.03,['EmbedZTT'])
      norm_systs['l_id'] = (0.02,['EmbedZTT','TTT','VVT','ZL'])
    else: 
      norm_systs['tau_id'] = (0.02,['ZTT'])
      norm_systs['l_id'] = (0.02,['ZTT','TTT','VVT','ZL'])
    norm_systs['l_fakerate'] = (0.2,['ZL'])
  if options.channel in ['em']:
    if options.embedding: 
      norm_systs['m_id'] = (0.02,['EmbedZTT','TTT','VVT','VVJ','ZLL','W'])
      norm_systs['e_id'] = (0.02,['EmbedZTT','TTT','VVT','VVJ','ZLL','W'])
      norm_systs['m_trg'] = (0.02,['EmbedZTT','TTT','VVT','VVJ','ZLL','W'])
      norm_systs['e_trg'] = (0.02,['EmbedZTT','TTT','VVT','VVJ','ZLL','W'])
    else:
      norm_systs['m_id'] = (0.02,['ZTT','TTT','TTJ','VVT','VVJ','ZLL','W'])
      norm_systs['e_id'] = (0.02,['ZTT','TTT','TTJ','VVT','VVJ','ZLL','W'])
      norm_systs['m_trg'] = (0.02,['ZTT','TTT','TTJ','VVT','VVJ','ZLL','W'])
      norm_systs['e_trg'] = (0.02,['ZTT','TTT','TTJ','VVT','VVJ','ZLL','W'])
    norm_systs['l_fakerate'] = (0.2,['ZLL','W'])
  
  for syst in norm_systs:
  
    val = norm_systs[syst][0]
    procs = norm_systs[syst][1]
    h1 = h0.Clone()
    h2 = h0.Clone()
    h1.SetName(h0.GetName()+syst)
    h2.SetName(h0.GetName()+syst)
    for p in procs:
      print p
      hup = directory.Get(p).Clone()
      hnom = directory.Get(p).Clone()
      hup.Scale(val)
      #print syst, p, hnom.Integral(), hup.Integral()
      h1 = h0.Clone()
      h2 = h0.Clone()
      h1.Add(hup)
      h2.Add(hup,-1)
    hists.append(h1.Clone())
    hists.append(h2.Clone())

  (uncert, up, down) = TotalUnc(h0, hists)
  outfile.cd(nodename)
  uncert.Write()
  up.Write()
  down.Write()

outfile.Close()


if options.do_unrolling==0: 
  print "Finished Processing"
  exit(0)
if is_2d and not options.do_unrolling:
  print "Finished Processing"
  exit(0) # add options for is_3d as well!
plot_file = ROOT.TFile(output_name, 'READ')

#if options.method in [12,16] or (options.channel != "tt" and options.method == "18"):
#    w_os = plot_file.Get(nodename+"/W.subnodes/W_os")    
#    w_ss = plot_file.Get(nodename+"/W.subnodes/W_ss")
#    w_os_error=ROOT.Double(0.)
#    w_ss_error=ROOT.Double(0.)
#    w_os_total = w_os.IntegralAndError(0,w_os.GetNbinsX()+1,w_os_error)
#    w_ss_total = w_ss.IntegralAndError(0,w_ss.GetNbinsX()+1,w_ss_error)
#    w_os_ss = w_os_total/w_ss_total
#    w_os_ss_error = math.sqrt( (w_os_error/w_os_total)**2 + (w_ss_error/w_ss_total)**2 )*w_os_ss
#
#    print "W OS/SS ratio = ", w_os_ss, "+/-", w_os_ss_error, "("+str(100*w_os_ss_error/w_os_ss)+" %)"

if options.custom_uncerts_wt_up != "" and options.custom_uncerts_wt_down != "": 
    custom_uncerts_up_name = "total_bkg_custom_uncerts_up"
    custom_uncerts_down_name = "total_bkg_custom_uncerts_down"
elif options.custom_uncerts_up_name != '':
    custom_uncerts_up_name = options.custom_uncerts_up_name
    custom_uncerts_down_name = options.custom_uncerts_down_name


if not options.no_plot:
    if options.extra_name != '': vname = options.extra_name
    else: vname = var_name

    if options.datacard != "": plot_name = options.outputfolder+'/'+vname+'_'+options.datacard+'_'+options.channel+'_'+options.year
    else: plot_name = options.outputfolder+'/'+vname+'_'+options.cat+'_'+options.channel+'_'+options.year
    if options.do_ss: plot_name += "_ss"
    if options.log_x: plot_name += "_logx" 
    if options.log_y: plot_name += "_logy"
    titles = plotting.SetAxisTitles(options.var,options.channel)
    if options.do_unrolling and is_2d: titles2d = plotting.SetAxisTitles2D(options.var,options.channel)
    if options.x_title == "": 
      x_title = titles[0]
      if options.do_unrolling and is_2d: 
        x_title = titles2d[0][0]
    else: x_title = options.x_title
    
    if options.y_title == "": 
        y_title = titles[1]
        if options.do_unrolling and is_2d: 
          if options.norm_bins: y_title = titles2d[0][1]
          else: y_title = titles2d[0][1]
          y_var_titles = titles2d[1]
    else: y_title = options.y_title
    scheme = options.channel
    if compare_w_shapes: scheme = 'w_shape'
    if compare_qcd_shapes: scheme = 'qcd_shape'
    if options.scheme != "": scheme = options.scheme
    FF = options.method in [17,18]
    if options.ml_ff:
      options.w_ff_closure = False
      options.qcd_ff_closure = False
    if "zttEmbed" in options.cat or "jetFakes" in options.cat:
        options.blind = False
        options.x_blind_min = -1e5
        options.x_blind_max = -1e5 
    if options.do_unrolling and is_2d:
        auto_blind=False
        options.norm_bins=False
        plotting.HTTPlotUnrolled(nodename, 
        plot_file, 
        options.signal_scale, 
        options.draw_signal_mass,
        FF,
        options.norm_bins,
        options.channel,
        options.blind,
        options.x_blind_min,
        options.x_blind_max,
        auto_blind,
        options.ratio,
        options.log_y,
        options.log_x,
        options.ratio_range,
        options.custom_x_range,
        options.x_axis_min,
        options.x_axis_max,
        options.custom_y_range,
        options.y_axis_max,
        options.y_axis_min,
        x_title,
        y_title,
        options.extra_pad,
        options.do_custom_uncerts,
        options.add_stat_to_syst,
        options.add_flat_uncert,
        options.uncert_title,
        options.lumi,
        plot_name,
        custom_uncerts_up_name,
        custom_uncerts_down_name,
        scheme,
        options.cat,
        x_lines,
        [y_labels,y_var_titles],
        options.embedding,
        vbf_background,
        options.signal_scheme
        )
    elif scheme != 'signal':
      auto_blind=False
      plotting.HTTPlot(nodename, 
        plot_file, 
        options.signal_scale, 
        options.draw_signal_mass,
        FF,
        options.norm_bins,
        options.channel,
        options.blind,
        options.x_blind_min,
        options.x_blind_max,
        options.ratio,
        options.threePads,
        options.ratio_log_y,
        options.log_y,
        options.log_x,
        options.ratio_range,
        options.custom_x_range,
        options.x_axis_min,
        options.x_axis_max,
        options.custom_y_range,
        options.y_axis_max,
        options.y_axis_min,
        x_title,
        y_title,
        options.extra_pad,
        options.signal_scheme,
        options.do_custom_uncerts,
        options.add_stat_to_syst,
        options.add_flat_uncert,
        options.uncert_title,
        options.lumi,
        plot_name,
        custom_uncerts_up_name,
        custom_uncerts_down_name,
        scheme,
        options.embedding,
        vbf_background,
        options.split_sm_scheme,
        options.ggh_scheme,
        options.cat,
        split_taus,
        auto_blind,
        discrete_x_axis,
        discrete_x_labels,
        options.qcd_ff_closure,
        options.w_ff_closure,
        options.bkg_comp,
        options.plot_signals.split(",")
        )
    else:    
      plotting.HTTPlotSignal(nodename, 
        plot_file, 
        options.signal_scale, 
        options.draw_signal_mass,
        options.norm_bins,
        options.channel,
        options.blind,
        options.x_blind_min,
        options.x_blind_max,
        options.ratio,
        options.log_y,
        options.log_x,
        options.ratio_range,
        options.custom_x_range,
        options.x_axis_min,
        options.x_axis_max,
        options.custom_y_range,
        options.y_axis_max,
        options.y_axis_min,
        x_title,
        y_title,
        options.extra_pad,
        options.signal_scheme,
        options.do_custom_uncerts,
        options.add_stat_to_syst,
        options.add_flat_uncert,
        options.uncert_title,
        options.lumi,
        plot_name,
        custom_uncerts_up_name,
        custom_uncerts_down_name
        )
    
   # plotting.SoverBPlot(nodename,
   #          plot_file,
   #          options.channel,
   #          options.log_y,
   #          options.log_x,
   #          options.custom_x_range,
   #          options.x_axis_max,
   #          options.x_axis_min,
   #          options.custom_y_range,
   #          options.y_axis_max,
   #          options.y_axis_min,
   #          x_title,
   #          options.extra_pad,
   #          plot_name+'_soverb')
    
   # hists = [plot_file.Get(nodename+"/bbH-LO700"), plot_file.Get(nodename+"/bbH700") ]
   # plotting.CompareHists(hists,
   #          ['Pythia','amc@NLO'],
   #          "bb#phi 700",
   #          options.ratio,
   #          options.log_y,
   #          options.log_x,
   #          options.ratio_range,
   #          options.custom_x_range,
   #          options.x_axis_max,
   #          options.x_axis_min,
   #          options.custom_y_range,
   #          options.y_axis_max,
   #          options.y_axis_min,
   #          x_title,
   #          y_title,
   #          options.extra_pad,
   #          False,
   #          plot_name,
   #          "#mu#tau_{h}")

#norm signal yields on datacards to 1pb AFTER plotting    
outfile =  ROOT.TFile(output_name, 'UPDATE')
for add_name in add_names: 
    if options.analysis in ['mssm','mssmrun2']:
        NormSignals(outfile,add_name)

# for smsummer16 need to ad WplusH and WminusH templates into one
if options.era in ["smsummer16",'cpsummer16','cpdecay16',"legacy16",'UL_16_preVFP','UL_16_postVFP','cpsummer17','UL_17','cp18','UL_18','mvadm2016'] and options.channel != 'zmm' and options.analysis != "mssmrun2":
  outfile.cd(nodename)
  directory = outfile.Get(nodename)
  hists_to_add = []
  for key in directory.GetListOfKeys():
    hist_name = key.GetName()  
    if 'wminush' in hist_name.lower():
      hist_to_add_name = hist_name.replace('minus', 'plus')  
      hist = directory.Get(hist_name).Clone()
      hist_to_add = directory.Get(hist_to_add_name).Clone()
      hist.Add(hist_to_add)
      hist.SetName(hist_name.replace('minus',''))
      hists_to_add.append(hist)
  for hist in hists_to_add: hist.Write()

if options.analysis in ['mssmrun2','vlq']:
  RenameMSSMrun2Datacards(outfile)

outfile.Close()

print "Finished Processing"
