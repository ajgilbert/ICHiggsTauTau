#!/usr/bin/env python

# importing packages
import sys
from optparse import OptionParser
import os
import numpy as np

config_files = {'2016_preVFP':'scripts/plot_UL_2016_preVFP.cfg',
		'2016_postVFP':'scripts/plot_UL_2016_postVFP.cfg',
		'2016preUL':'scripts/plot_mssm_2016.cfg',
                '2017preUL':'scripts/plot_mssm_2017.cfg',
                '2017':'scripts/plot_UL_2017.cfg',
                '2018':'scripts/plot_UL_2018.cfg',
                '2018preUL':'scripts/plot_mssm_2018.cfg'}
               
def split_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

def CreateBatchJob(name,cmssw_base,cmd_list):
  if os.path.exists(job_file): os.system('rm %(name)s' % vars())
  os.system('echo "#!/bin/bash" >> %(name)s' % vars())
  os.system('echo "cd %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2" >> %(name)s' % vars())
  os.system('echo "source /vols/grid/cms/setup.sh" >> %(name)s' % vars())
  os.system('echo "export SCRAM_ARCH=slc6_amd64_gcc481" >> %(name)s' % vars())
  os.system('echo "eval \'scramv1 runtime -sh\'" >> %(name)s' % vars())
  os.system('echo "source %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/setup_libs.sh" >> %(name)s' % vars())
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

def BINS(start,end,interval,offset):
  bins = np.arange(start,end+offset,interval)
  bins_string = ','.join(str(x) for x in bins)
  bins  = []
  bins.append(bins_string)
  return bins


parser = OptionParser()
parser.add_option('--channel',help= 'Name of input channels', default='mt')
parser.add_option("--years", dest="years", type='string', default='2016_preVFP,2016_postVFP,2017,2018',help="Year input")
parser.add_option('--output_folder', help= 'Name of output folder to create', default='tau_SF_datacards')
parser.add_option("--batch", dest="batch", action='store_true', default=False, help="Submit on batch.")
parser.add_option("--systs", dest="systs", action='store_true', default=False, help="Add systematic variations.")
parser.add_option('--wp', help= 'Name of VsJet WP to use', default='medium')
(options, args) = parser.parse_args()

# initialising variables
output_folder = options.output_folder
channels = options.channel.split(',')
years = options.years.split(',')
wp=options.wp

print 'Processing channels:      %(channels)s' % vars()
print 'Processing years:         %(years)s' % vars()

# cmssw directory
cmssw_base = os.getcwd().replace('src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2','')

# check whether path is in an exitisting directory
if not os.path.isdir('%(output_folder)s' % vars()):
  os.system("mkdir %(output_folder)s" % vars())
if not os.path.isdir('%(output_folder)s/jobs' % vars()):
  os.system("mkdir %(output_folder)s/jobs" % vars())
  
categories_mt = ["inclusive",'dm0','dm1','dm10','dm11']
categories_zmm = ["inclusive"]

cat_schemes = {
               'mt' : categories_mt,
	       'zmm': categories_zmm,}
var_mt = [
         #['m_vis',BINS(0,300,10,1)], 
         #['pt_2', BINS(0,200,5,1)],
         ['pt_2', '(180,20,200)'],
         ]        
var_zmm = [
         ['m_vis',BINS(0,300,5,1)], 
          ]

var_schemes = {
               'mt' : var_mt,
	       'zmm': var_zmm}

for year in years:

  mt_sep_shape_systematics = []
  common_shape_systematics = []

  if options.systs:
    mt_sep_shape_systematics = [
        ' --syst_tau_scale_grouped=\'CMS_scale_t_*group_%(year)s\'' % vars(), # Tau energy scale
        #' --syst_res_j=\'CMS_res_j_%(year)s\'' % vars(), # Jet energy resolution
        #' --syst_scale_met_unclustered=\'CMS_scale_met_unclustered_%(year)s\'' % vars(), # MET unclustered energy uncertainty
        #' --syst_scale_j=\'CMS_scale_j_%(year)s\'' % vars(), # Jet energy scale (grouped)
        #' --syst_mufake_0pi_scale="CMS_ZLShape_mt_1prong_%(year)s"' % vars(), # l to tau h fake energy scale
        #' --syst_mufake_1pi_scale="CMS_ZLShape_mt_1prong1pizero_%(year)s"' % vars(), # l to tau h fake energy scale
        #' --syst_qcd_bkg="CMS_QCD_BackgroundSubtraction"',
        ' --syst_jfake_scale="CMS_scale_jfake"',
        #' --syst_scale_j_regrouped="CMS_scale_j_*group"', # Jet energy scale (grouped)
      ]
    common_shape_systematics = [
        ' --syst_tquark=\'CMS_htt_ttbarShape\'', # Top pT re-weighting
        ' --syst_zwt=\'CMS_htt_dyShape\''
      ]
  

  if not os.path.isdir('%(output_folder)s/%(year)s' % vars()):
    os.system("mkdir %(output_folder)s/%(year)s" % vars())
  CFG = config_files[year]
  for ch in channels:
    if not os.path.isdir('%(output_folder)s/%(year)s/%(ch)s' % vars()):
      os.system("mkdir %(output_folder)s/%(year)s/%(ch)s" % vars())

    method='12' # this will relax the iso for the QCD shape to make sure it is not as strongly influenced by the background modelling
    sf='1.'

    if ch=='tt': 
      add_cond = '--add_wt=\'wt_prefire*%(sf)s/wt_tau_id_dm\'' % vars()
    elif ch in ['mt','et']: add_cond = '--add_wt=\'wt_prefire*%(sf)s/wt_tau_id_pt\'' % vars()
    else: add_cond = '--add_wt=\'wt_prefire\'' % vars()

    categories = cat_schemes[ch]
    variables = var_schemes[ch]
    for cat in categories:
      for item in variables:
           var_used = item[0]
           bin_used = item[1]
           run_cmd=''
           run_cmd_2=''
           run_cmd_3=''
           run_cmd_4=''
           run_cmd_5=''

           extra=''
           if ',' in var_used: extra = '--do_unrolling=0' 


           if options.systs:
             for s in common_shape_systematics: extra+=' %s' % s
 
           extra_mt = extra

           #if options.systs:
           #  for s in mt_sep_shape_systematics: extra_mt+=' %s' % s

           systs = ['']
           if ch == 'mt': systs+=mt_sep_shape_systematics

           for s in systs: 
             if s !='':
               extra_mt=' %s --no_default' % s 
 
             if ',' in var_used: extra_mt += ' --do_unrolling=0'       

             s_name=''
             if '--' in s and '=' in s: s_name='_'+s.split('--')[1].split('=')[0]
             if s_name!='':
               extra_mt+='  --extra_name=%(s_name)s' % vars() 
               extra+='  --extra_name=%(s_name)s' % vars() 

             if ch == 'mt':

               cuts='m_vis>40&&m_vis<80&&pt_1>25&&trg_singlemuon&&fabs(eta_2)<2.1'
               if '2016' in year: cuts=cuts.replace('trg_singlemuon','trg_singlemuon24') # make sure you use the un-prescaled trigger that matches the mu+tau monitoring path selections 
               if '2017' in year: cuts=cuts.replace('trg_singlemuon','trg_singlemuon27') # make sure you use the un-prescaled trigger that matches the mu+tau monitoring path selections 
               if year == '2017': cuts+='&&pt_1>28' # make sure muon is above 27 GeV trigger threshold as 24 GeV trigger was prescaled in 4/fb of datataking 

 
#               run_cmd = 'python %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<40)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt40" --wp=%(wp)s ' % vars()
  
               run_cmd_2 = 'python %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<30)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&iso_1<0.1)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt30" --wp=%(wp)s ' % vars()

#               run_cmd_3 = 'python %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<40)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&(trg_tt_monitoring_1||trg_tt_monitoring_2||trg_tt_monitoring_3||trg_tt_monitoring_4||trg_tt_monitoring_5||trg_tt_monitoring_6))\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt40_passTRG" --wp=%(wp)s ' % vars()

               run_cmd_4 = 'python %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<30)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&iso_1<0.1&&(trg_tt_monitoring_l1>32&&trg_tt_monitoring_1||trg_tt_monitoring_2||trg_tt_monitoring_3||trg_tt_monitoring_4||trg_tt_monitoring_5||trg_tt_monitoring_6))\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt30_passTRG" --wp=%(wp)s ' % vars()
             # remove trg_tt_monitoring_l1 after 2017 is rerun!!!
  
             elif ch == 'zmm':
              run_cmd = 'python %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --set_alias=\'%(cat)s:({%(cat)s}&&m_vis>50&&pt_1>25&&trg_singlemuon)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra)s ' % vars()
  
             commands = [run_cmd, run_cmd_2, run_cmd_3, run_cmd_4]
             if ',' in var_used: var_name=var_used.replace(',','_vs_')
             else: var_name=var_used 
 
             for i in range(0,len(commands)):
               if commands[i] is '': continue
               #s_name=''
               #if '--' in s and '=' in s: s_name='_'+s.split('--')[1].split('=')[0]
               job_file = '%(output_folder)s/jobs/%(var_name)s_%(cat)s_%(ch)s_%(year)s%(s_name)s_%(i)i.sh' % vars()
               CreateBatchJob(job_file,cmssw_base,[commands[i]])
               SubmitBatchJob(job_file,time=180,memory=24,cores=1)
  