#!/usr/bin/env python

# python scripts/mssm_control_plots.py --channel='mt,et,tt' --year='2018' --output='control_plots_mssm_new' --control_plots=True --ff_closures=True --add_systs=True --dry_run=True

import os
import argparse
import math as math

parser = argparse.ArgumentParser()
parser.add_argument('--channel',help= 'Name of input channels', default='mttt,ettt,tttt,eett,mmtt,emtt')
parser.add_argument('--year', help= 'Name of input years', default='2016,2017,2018')
parser.add_argument('--output', help= 'Name of output folder to create', default='4tau_plots')
args = parser.parse_args()
      
# Things to loop over
channels = args.channel.split(',')
years = args.year.split(',')
output = args.output

def myround(x, base=5):
    return base * round(x/base)

def GetBinning(var,min_bin,max_bin,num_bins,log=False,round=1,final_bin_start=None):
  bins = []
  if not log:
    if final_bin_start == None:
      range_bins = max_bin - min_bin
    else:
      range_bins = final_bin_start - min_bin
    bin_interval = range_bins/num_bins
    for i in range(0,num_bins+1):
      bins.append(myround(min_bin+(i*bin_interval),base=round))
  else:
    if min_bin == 0: min_bin = 0.001
    min_bin_log = math.log(min_bin,10)
    if final_bin_start == None:
      max_bin_log = math.log(max_bin,10)
      range_bins_log = max_bin_log - min_bin_log
    else:
      final_bin_start_log = math.log(final_bin_start,10)
      range_bins_log = final_bin_start_log - min_bin_log
    bin_interval_log = range_bins_log/num_bins
    for i in range(0,num_bins+1):
      bins.append(myround(10**(min_bin_log+(i*bin_interval_log)),base=round))

  if final_bin_start != None: bins.append(max_bin)
  return var+str(bins).strip()

def CreateBatchJob(name,cmssw_base,cmd_list):
  if os.path.exists(job_file): os.system('rm %(name)s' % vars())
  os.system('echo "#!/bin/bash" >> %(name)s' % vars())
  os.system('echo "cd %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/4tau" >> %(name)s' % vars())
  os.system('echo "source /vols/grid/cms/setup.sh" >> %(name)s' % vars())
  os.system('echo "export SCRAM_ARCH=slc6_amd64_gcc481" >> %(name)s' % vars())
  os.system('echo "eval \'scramv1 runtime -sh\'" >> %(name)s' % vars())
  os.system('echo "source %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/4tau/scripts/setup_libs.sh" >> %(name)s' % vars())
  os.system('echo "ulimit -c 0" >> %(name)s' % vars())
  for cmd in cmd_list:
    os.system('echo "%(cmd)s" >> %(name)s' % vars())
  os.system('chmod +x %(name)s' % vars())
  print "Created job:",name

def SubmitBatchJob(name,time=180,memory=24,cores=1):
  error_log = name.replace('.sh','_error.log')
  output_log = name.replace('.sh','_output.log')
  if os.path.exists(error_log): os.system('rm %(error_log)s' % vars())
  if os.path.exists(output_log): os.system('rm %(output_log)s' % vars())
  if cores>1: os.system('qsub -e %(error_log)s -o %(output_log)s -V -q hep.q -pe hep.pe %(cores)s -l h_rt=0:%(time)s:0 -l h_vmem=%(memory)sG -cwd %(name)s' % vars())
  else: os.system('qsub -e %(error_log)s -o %(output_log)s -V -q hep.q -l h_rt=0:%(time)s:0 -l h_vmem=%(memory)sG -cwd %(name)s' % vars())

all_ch_variables = [
                    #GetBinning('pt_1',0,400,100,round=1),
                    #GetBinning('pt_2',0,400,100,round=1),
                    #GetBinning('pt_3',0,400,100,round=1),
                    #GetBinning('mt_lep_12',0,400,100,round=1),
                    #GetBinning('mt_lep_13',0,400,100,round=1),
                    #GetBinning('mt_lep_23',0,400,100,round=1),
                    #GetBinning('mvis_12',0,400,100,round=1),
                    #GetBinning('mvis_13',0,400,100,round=1),
                    #GetBinning('mvis_23',0,400,100,round=1),
                    GetBinning('mvis_min_sum_dR_1',0,500,100,round=1),
                    #GetBinning('mvis_min_sum_dR_2',0,300,60,round=1),
                    #GetBinning('pt_min_dphi_1',0,300,60,round=1),
                    #GetBinning('pt_min_dphi_2',0,300,60,round=1),
                    #"mvis_min_sum_dR_1,mvis_min_sum_dR_2[0,50,100,150,200,250,300,350,400],[0,50,100,150,200,250,300,350,400]",
                    ]

ch_dep_var = {"mttt":[],
              "ettt":[],
              "mmtt":[],
              "eett":[],
              "tttt":[],
              "emtt":[],
	      "ttt":[],
              }

unb_ch_dep_var = {"mttt":[],
              "ettt":[],
              "mmtt":[],
              "eett":[],
              "tttt":[],
              "emtt":[],
	      "ttt":[],
              }


config_files = {'2016':'scripts/plot_UL_2016.cfg',
                '2017':'scripts/plot_UL_2017.cfg',
                '2018':'scripts/plot_UL_2018.cfg'
               }

categories = {
              "ttt":["inclusive","nobtag","btag"],
              "tttt":["inclusive","nobtag","btag"],
              "ettt":["inclusive","nobtag","btag"],
              "mttt":["inclusive","nobtag","btag"],
              "emtt":["inclusive","nobtag","btag"],
              "eett":["z_control","2l2t_sig","inclusive","nobtag","btag","z_control_nobtag","2l2t_sig_nobtag","z_control_btag","2l2t_sig_btag"],
              "mmtt":["z_control","2l2t_sig","inclusive","nobtag","btag","z_control_nobtag","2l2t_sig_nobtag","z_control_btag","2l2t_sig_btag"],
              }

add_options = ''
add_options = '--ratio_range=\'0,3\' --plot_signals=\'phi200A100To4Tau\'  --auto_rebinning --bin_uncert_fraction=0.15 --signal_scale=10'
#add_options = '--ratio_range=\'0,3\' --plot_signals=\'phi200A100To4Tau\'  --auto_rebinning --bin_uncert_fraction=0.15 --signal_scale=10 --vsjets=None --add_wt=\'1/(idisoweight_3*idisoweight_4)\''
 
blind_options = '--blind --x_blind_min=-999 --x_blind_max=999'

# Set up output directories

cmssw_base = os.getcwd()

if not os.path.isdir('%(cmssw_base)s/%(output)s' % vars()):
  os.system("mkdir %(cmssw_base)s/%(output)s" % vars())

if not os.path.isdir('%(cmssw_base)s/%(output)s/jobs' % vars()):
  os.system("mkdir %(cmssw_base)s/%(output)s/jobs" % vars())

for channel in channels:
  if not os.path.isdir('%(cmssw_base)s/%(output)s/%(channel)s' % vars()):
    os.system("mkdir %(cmssw_base)s/%(output)s/%(channel)s" % vars())
  for year in years:
    if not os.path.isdir('%(cmssw_base)s/%(output)s/%(channel)s/%(year)s' % vars()):
      os.system("mkdir %(cmssw_base)s/%(output)s/%(channel)s/%(year)s" % vars())


# Job loop

for year in years:
  for channel in channels:
    for cat in categories[channel]:
      variables = all_ch_variables+ch_dep_var[channel]+unb_ch_dep_var[channel]
      for var in variables:
        if '[' in var: var_string = var.split('[')[0]
        elif '(' in var: var_string = var.split('(')[0]
        if var in unb_ch_dep_var[channel] and "control" not in cat: continue
        if var not in unb_ch_dep_var[channel]:
          add_options += " " + blind_options
        output_folder = '%(cmssw_base)s/%(output)s/%(channel)s/%(year)s' % vars()
        cfg = config_files[year]
        run_cmd = 'python %(cmssw_base)s/scripts/4tauPlot.py --cfg=\'%(cfg)s\' --channel=\'%(channel)s\' --var=\'%(var)s\' %(add_options)s --outputfolder=\'%(output_folder)s\' --cat=\'%(cat)s\'' % vars()
        job_file = "%(cmssw_base)s/%(output)s/jobs/%(var_string)s_%(channel)s_%(cat)s_%(year)s.sh" % vars()
        CreateBatchJob(job_file,os.getcwd().replace('src/UserCode/ICHiggsTauTau/Analysis/4tau',''),[run_cmd])
        SubmitBatchJob(job_file,time=180,memory=24,cores=1)
           
    

