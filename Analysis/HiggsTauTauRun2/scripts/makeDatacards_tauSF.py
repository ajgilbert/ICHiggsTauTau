#!/usr/bin/env python

#python scripts/makeDatacards_tauSF.py --systs --batch --years 2016_preVFP,2016_postVFP,2017,2018 --output_folder tau_SF_UL_220223_medium_PFMet_oldpTetacuts_v1

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
                '2018':'scripts/plot_UL_2018.cfg', #changed for v2p5:scripts/plot_UL_2018_124.cfg
                '2018preUL':'scripts/plot_mssm_2018.cfg',
                '2022_preEE':'scripts/plot_2022_preEE.cfg',
                '2022_postEE':'scripts/plot_2022_postEE.cfg'}
               
def split_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

def CreateBatchJob(name,cmssw_base,cmd_list):
  if os.path.exists(job_file): os.system('rm %(name)s' % vars())
  os.system('echo "#!/bin/bash" >> %(name)s' % vars())
  os.system('echo "cd /vols/cms/ia2318/CMSSW_10_2_19//src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2" >> %(name)s' % vars())
  #os.system('echo "source /vols/grid/cms/setup.sh" >> %(name)s' % vars())
  os.system('echo "source /cvmfs/cms.cern.ch/cmsset_default.sh" >> %(name)s' % vars()) # don't know which one to keep
  os.system('echo "export SCRAM_ARCH=slc6_amd64_gcc481" >> %(name)s' % vars())
  os.system('echo "eval \`scramv1 runtime -sh\`" >> %(name)s' % vars())
  os.system('echo "source %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/setup_libs.sh" >> %(name)s' % vars())
  os.system('echo "ulimit -c 0" >> %(name)s' % vars())
  for cmd in cmd_list:
    os.system('echo "%(cmd)s" >> %(name)s' % vars())
  os.system('chmod +x %(name)s' % vars())
  print("Created job:",name)

def CreateCondorScript(name):
  name_sh = name
  name = name.replace('.sh', '.sub')
  if os.path.exists(name): os.system('rm %(name)s' % vars())
  output = name.replace('.sub','_output.log')
  error = name.replace('.sub','_error.log')
  os.system('echo "Universe = vanilla" >> %(name)s' % vars())
  #os.system('echo "getenv = True" >> %(name)s' % vars()) # introduced perl error
  os.system('echo "executable = %(name_sh)s" >> %(name)s' % vars())
  os.system('echo "output = %(output)s" >> %(name)s' % vars())
  os.system('echo "error = %(error)s" >> %(name)s' % vars())
  os.system('echo "+MaxRuntime = 6500" >> %(name)s' % vars())
  os.system('echo "+OS = \\"EL7\\"">> %(name)s' % vars())
  os.system('echo "queue" >> %(name)s' % vars())
  print("Created condor script:",name)

def SubmitBatchJob(name,time=180,memory=24,cores=1):
  # executable = name
  # output = name.replace('.sh','_output.log')
  # error = name.replace('.sh','_error.log')
  # if os.path.exists(error): os.system('rm %(error)s' % vars())
  # if os.path.exists(output): os.system('rm %(output)s' % vars())
  os.system('condor_submit %(name)s' % vars())
  # error_log = name.replace('.sh','_error.log')
  # output_log = name.replace('.sh','_output.log')
  # if os.path.exists(error_log): os.system('rm %(error_log)s' % vars())
  # if os.path.exists(output_log): os.system('rm %(output_log)s' % vars())
  # if cores>1: os.system('qsub -e %(error_log)s -o %(output_log)s -V -q hep.q -pe hep.pe %(cores)s -l h_rt=0:%(time)s:0 -l h_vmem=%(memory)sG -cwd %(name)s' % vars())
  # else: os.system('qsub -e %(error_log)s -o %(output_log)s -V -q hep.q -l h_rt=0:%(time)s:0 -l h_vmem=%(memory)sG -cwd %(name)s' % vars())

def BINS(start,end,interval,offset):
  bins = np.arange(start,end+offset,interval)
  bins_string = ','.join(str(x) for x in bins)
  bins  = []
  bins.append(bins_string)
  return bins


parser = OptionParser()
parser.add_option('--channel',help= 'Name of input channels', default='mt,zmm')
parser.add_option("--years", dest="years", type='string', default='2016_preVFP,2016_postVFP,2017,2018',help="Year input")
parser.add_option('--output_folder', help= 'Name of output folder to create', default='tau_SF_datacards')
parser.add_option('--wp', help= 'Name of VsJet WP to use', default='medium')
parser.add_option("--batch", dest="batch", action='store_true', default=False, help="Submit on batch.")
parser.add_option("--systs", dest="systs", action='store_true', default=False, help="Add systematic variations.")
parser.add_option("--v2p5", dest="v2p5", action='store_true', default=False, help="Use DeepTau v2p5")
(options, args) = parser.parse_args()

# initialising variables
output_folder = options.output_folder
channels = options.channel.split(',')
years = options.years.split(',')
wp=options.wp

print('Processing channels:      %(channels)s' % vars())
print('Processing years:         %(years)s' % vars())

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
         ['pt_2,m_vis', '(180,20,200),(16,40,200)'],
         ]        
var_zmm = [
         ['m_vis',BINS(50,150,5,1)], 
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
        ' --syst_jfake_scale="CMS_scale_jfake"',
        ' --syst_scale_met_unclustered=\'CMS_scale_met_unclustered_%(year)s\'' % vars(), # MET unclustered energy uncertainty
        ' --syst_mufake_0pi_scale="CMS_ZLShape_mt_1prong_%(year)s"' % vars(), # l to tau h fake energy scale
        ' --syst_mufake_1pi_scale="CMS_ZLShape_mt_1prong1pizero_%(year)s"' % vars(), # l to tau h fake energy scale
        ' --syst_qcd_bkg="CMS_QCD_BackgroundSubtraction"',]
        #' --syst_scale_j=\'CMS_scale_j_%(year)s\'' % vars(), # Jet energy scale (grouped)
    if year == "2022_preEE" or year == "2022_postEE":
           mt_sep_shape_systematics.extend([
              ' --syst_scale_met_clustered=\'CMS_scale_met_clustered_%(year)s\'' % vars(),]) # MET clustered energy uncertainty
    else:
           mt_sep_shape_systematics.extend([
              ' --syst_res_j=\'CMS_res_j_%(year)s\'' % vars(), # Jet energy resolution
           ' --syst_scale_j_regrouped="CMS_scale_j_*group"', # Jet energy scale (grouped)
      ])
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

    method='9' # this will relax the iso for the QCD shape to make sure it is not as strongly influenced by the background modelling
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
           run_cmd_6=''
           run_cmd_7=''
           run_cmd_8=''

           extra=''
           if ',' in var_used: extra = '--do_unrolling=0' 


           if options.systs:
             for s in common_shape_systematics: extra+=' %s' % s
 
           extra_mt = extra

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
 
               cuts='m_vis>40&&fabs(eta_1)<2.1&&trg_singlemuon'

               if options.v2p5:
                  run_cmd_3 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<65)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt65" --wp=%(wp)s --v2p5' % vars()
                  run_cmd_4 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<65)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&deepTauVsEle_tight_2>0.5)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt65_tightVsEle" --wp=%(wp)s --v2p5' % vars()
                  run_cmd_5 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1>70)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTGt70" --wp=%(wp)s --v2p5' % vars()
                  run_cmd_6 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1>70)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&deepTauVsEle_tight_2>0.5)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTGt70_tightVsEle" --wp=%(wp)s --v2p5' % vars()
                  run_cmd_7 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<65)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s)\' --method=8 --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt65_aiso" --wp=%(wp)s --do_aiso --v2p5' % vars()
                  run_cmd_8 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<65)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&deepTauVsEle_tight_2>0.5)\' --method=8 --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt65_tightVsEle_aiso" --wp=%(wp)s --do_aiso --v2p5' % vars()

               else:
                  run_cmd_3 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<65)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt65" --wp=%(wp)s ' % vars()
                  run_cmd_4 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<65)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&deepTauVsEle_tight_2>0.5)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt65_tightVsEle" --wp=%(wp)s ' % vars()
                  run_cmd_5 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1>70)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTGt70" --wp=%(wp)s' % vars()
                  run_cmd_6 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1>70)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&deepTauVsEle_tight_2>0.5)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTGt70_tightVsEle" --wp=%(wp)s' % vars()
                  run_cmd_7 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<65)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s)\' --method=8 --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt65_aiso" --wp=%(wp)s --do_aiso' % vars()
                  run_cmd_8 = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --sel=\'(mt_1<65)\' --set_alias=\'%(cat)s:({%(cat)s}&&%(cuts)s&&deepTauVsEle_tight_2>0.5)\' --method=8 --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra_mt)s --datacard="%(cat)s_mTLt65_tightVsEle_aiso" --wp=%(wp)s --do_aiso' % vars()

               commands = [run_cmd_3, run_cmd_4, run_cmd_5, run_cmd_6, run_cmd_7,run_cmd_8]
               if ',' in var_used: var_name=var_used.replace(',','_vs_')

               for i in range(0,len(commands)):
                 if commands[i] == '': continue
                 job_file = '%(output_folder)s/jobs/%(var_name)s_%(cat)s_%(ch)s_%(year)s%(s_name)s_%(i)i.sh' % vars()
                 CreateBatchJob(job_file,cmssw_base,[commands[i]])
                 CreateCondorScript(job_file)
                 if 'syst_scale_j_regrouped' in s_name: SubmitBatchJob(job_file.replace('.sh', '.sub'),time=500,memory=24,cores=1)
                 else: SubmitBatchJob(job_file.replace('.sh', '.sub'),time=180,memory=24,cores=1)  
 
             elif ch == 'zmm':
                run_cmd = 'python2 %(cmssw_base)s/src/UserCode/ICHiggsTauTau/Analysis/HiggsTauTauRun2/scripts/HiggsTauTauPlot.py --cfg=%(CFG)s --channel=%(ch)s --set_alias=\'%(cat)s:({%(cat)s}&&m_vis>50&&pt_1>25&&fabs(eta_1)<2.1&&trg_singlemuon)\' --method=%(method)s --cat=%(cat)s --outputfolder=%(output_folder)s/%(year)s/%(ch)s --ggh_masses_powheg='' --bbh_masses_powheg='' --var=\'%(var_used)s%(bin_used)s\' %(add_cond)s --ratio_range="0.5,1.5" %(extra)s ' % vars()
                commands = [run_cmd]
                var_name = var_used
                 
                for i in range(0,len(commands)):
                  if commands[i] == '': continue
                  job_file = '%(output_folder)s/jobs/%(var_name)s_%(cat)s_%(ch)s_%(year)s%(s_name)s_%(i)i.sh' % vars()
                  CreateBatchJob(job_file,cmssw_base,[commands[i]])
                  CreateCondorScript(job_file)
                  if 'syst_scale_j_regrouped' in s_name: SubmitBatchJob(job_file.replace('.sh', '.sub'),time=500,memory=24,cores=1)
                  else: SubmitBatchJob(job_file.replace('.sh', '.sub'),time=180,memory=24,cores=1)
