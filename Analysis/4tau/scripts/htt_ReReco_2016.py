#!/usr/bin/env python

# python scripts/htt_ReReco_2016.py --bkg --data --jetmetuncerts --scales="default,scale_t_0pi,scale_t_1pi,scale_t_3prong,scale_t_3prong1pi0,scale_efake_0pi,scale_efake_1pi,scale_mufake_0pi,scale_mufake_1pi,scale_e" --submit='./scripts/submit_ic_batch_job.sh "hep.q -l h_rt=0:180:0 -l h_vmem=24G"'

import sys
from optparse import OptionParser
import os
import math
import json

os.system("voms-proxy-init --voms cms --valid 96:00 --out ~/cms.proxy")

JOBWRAPPER      = './scripts/generate_job.sh'
JOBSUBMIT       = 'true'
if "JOBWRAPPER" in os.environ:      JOBWRAPPER      = os.environ["JOBWRAPPER"]
if "JOBSUBMIT"  in os.environ:      JOBSUBMIT       = os.environ["JOBSUBMIT"]
print "Using job-wrapper:    " + JOBWRAPPER
print "Using job-submission: " + JOBSUBMIT

CONDOR_TEMPLATE = """executable = ./jobs/%(EXE)s
Proxy_path =/afs/cern.ch/user/a/adow/private/x509up
arguments = $(Proxy_path)
output                = /dev/null
error                 = /dev/null
log                   = ./jobs/%(TASK)s.$(ClusterId).log
requirements = (OpSysAndVer =?= "SLCern6")
+JobFlavour     = "longlunch"
queue
"""



def split_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))


parser = OptionParser()
parser.add_option("--wrapper", dest="wrapper",
                  help="Specify the job-wrapper script. The current wrapper is '%(JOBWRAPPER)s'."
                  " Using the --wrapper option overrides both the default and the environment variable. " % vars())
parser.add_option("--submit", dest="submit",
                  help="Specify the job-submission method. The current method is '%(JOBSUBMIT)s'"
                  " Using the --submit option overrides both the default and the environment variable. " % vars())
parser.add_option("--data", dest="proc_data", action='store_true', default=False,
                  help="Process data samples")
parser.add_option("--calc_lumi", dest="calc_lumi", action='store_true', default=False,
                  help="Run on data and only write out lumi mask jsons")
parser.add_option("--bkg", dest="proc_bkg", action='store_true', default=False,
                  help="Process background mc samples")
parser.add_option("--all", dest="proc_all", action='store_true', default=False,
                  help="Process all samples")
parser.add_option("--list_backup", dest="slbackupname", type='string', default='prevlist',
                  help="Name you want to give to the previous files_per_samples file, in case you're resubmitting a subset of jobs")
parser.add_option("--scales", dest="scales", type='string', default='default',
                  help="List of systematic shifts to process")
parser.add_option("--analysis", dest="analysis", type='string', default='sm',
                  help="Specify whether trees are produced for mssm or sm analysis")
parser.add_option("--parajobs", dest="parajobs", action='store_true', default=False,
                  help="Submit jobs parametrically")
parser.add_option("--config", dest="config", type='string', default='',
                  help="Config file")
parser.add_option("--condor", action='store_true', default=False,
                  help="Submit jobs to condor (for lxplus)")
parser.add_option("--jetmetuncerts", dest="jetmetuncerts", action='store_true', default=False,
                  help="Do JES, JER, and MET uncertainties")
parser.add_option("--batch", dest="batch", action='store_true', default=False,
                  help="Submit to imperial batch")

(options, args) = parser.parse_args()
if options.batch: options.submit='./scripts/submit_ic_batch_job.sh "hep.q -l h_rt=0:180:0 -l h_vmem=24G"'
if options.wrapper: JOBWRAPPER=options.wrapper
if options.submit:  JOBSUBMIT=options.submit
if options.condor: JOBWRAPPER = "./scripts/generate_condor_job.sh"

jetuncert_string=''
if options.jetmetuncerts:
  jetuncert_string='\\"do_jetmet_uncerts\\":true'
else:
  jetuncert_string='\\"do_jetmet_uncerts\\":false'

def getParaJobSubmit(N):
  if not options.submit: return 'true'
  sub_opts=JOBSUBMIT.split(' ',1)[1]
  if '\"' in sub_opts: sub_opts = sub_opts.replace('\"','')
  sub_opts = '%s \"%s -t 1-%i:1\"' %(JOBSUBMIT.split(' ',1)[0], sub_opts, N)
  print sub_opts
  return sub_opts

BACKUPNAME = options.slbackupname

#channels = options.channels
scales = options.scales
analysis = options.analysis
parajobs = options.parajobs


scale_list = scales.split(',')
flatjsonlist = []
flatjsonlistdysig = []
flatjsonlist.append("job:sequences:all:")
n_scales=0
for scale in scale_list: 
  n_scales+=1  
  if scale == "default":
    flatjsonlist.append("^%(scale)s"%vars())
    flatjsonlistdysig.append("^%(scale)s"%vars())
  else:
    n_scales+=1
    flatjsonlist.append("^%(scale)s_hi^%(scale)s_lo"%vars()) 
    flatjsonlistdysig.append("^%(scale)s_hi^%(scale)s_lo"%vars()) 
 
CONFIG='scripts/config_ReReco_2016.json'
if options.config != '': CONFIG = options.config

n_channels=1
output_folder = ""
svfit_folder = ""
svfit_mode = 0
with open(CONFIG,"r") as input:
  with open ("config_for_python_channels.json","w") as output:
    for line in input:
      if not '//' in line:
        output.write(line)
    output.close()
  input.close()

with open("config_for_python_channels.json") as config_file:
  cfg = json.load(config_file)
  n_channels=len(cfg["job"]["channels"])
  output_folder = cfg["sequence"]["output_folder"]
  svfit_mode = cfg["sequence"]["new_svfit_mode"]
  svfit_folder = cfg["sequence"]["svfit_folder"]

# makes sure output folder(s) (and svfit folder(s) if needed) is always created
os.system("bash scripts/make_output_folder.sh {}".format(output_folder))
if svfit_mode == 1:
    os.system("bash scripts/make_output_folder.sh {}".format(svfit_folder))
  
#scale = int(math.ceil(float(n_scales*n_channels)/50))
scale = int(math.ceil(float(n_scales*n_channels)/30))

if scale < 1: scale = 1

total = float(len(flatjsonlistdysig))
flatjsons = []
# this makes sure the JES's are submitted as seperate jobs (solves memory issues)
#for i in flatjsonlistdysig:
#  if 'scale_j' in i and 'hf' not in i and 'cent' not in i and 'full' not in i and 'relbal' not in i:
#    flatjsons.append('job:sequences:all:'+i)
#    flatjsonlistdysig.remove(i)
#    scale = int(math.ceil(float((n_scales-2)*n_channels)/50))
#    if scale < 1: scale = 1
if options.jetmetuncerts:
# when we do the jet met uncertainties we do not want to run additional systematics in the same job 
  for i in flatjsonlistdysig:
    if 'default' in i:
      flatjsons.append('job:sequences:all:'+i)
      flatjsonlistdysig.remove(i)
# split into seperate jobs if number of scales is over a value
for i in range(0,scale):
   first = i*int(math.ceil(total/scale))
   last = (i+1)*int(math.ceil(total/scale))
   temp=''.join(flatjsonlistdysig[first:last]) 
   if temp == '': continue
   temp='job:sequences:all:'+temp
   flatjsons.append(temp)
 
#FILELIST='filelists/Jan24_MC_102X_2016'
FILELIST='filelists/Sep18_2016_MC_102X'

signal_mc = [ ]
signal_vh = [ ] 
signal_mc_ww = [ ]

if os.path.isfile("./jobs/files_per_sample_2016.txt"):
  os.system("mv ./jobs/files_per_sample_2016.txt ./jobs/files_per_sample_2016-%(BACKUPNAME)s.txt"%vars())

file_persamp = open("./jobs/files_per_sample_2016.txt", "w")


if options.proc_data or options.proc_all or options.calc_lumi:
    with open(CONFIG,"r") as input:
      with open ("config_for_python.json","w") as output:
        for line in input:
          if not '//' in line:
            output.write(line)
        output.close()
      input.close()

    with open("config_for_python.json") as config_file:
      cfg = json.load(config_file)

    channels=cfg["job"]["channels"]

if options.proc_data or options.proc_all or options.calc_lumi:  

  data_samples = []
  data_eras = ['B','C','D','E','F','G','H']
  for chn in channels:
    for era in data_eras:
      if 'mttt' in chn:
          data_samples+=['SingleMuon'+era]
      #if 'mmtt' in chn:
      #    data_samples+=['DoubleMuon'+era]
      if 'ettt' in chn:
          data_samples+=['SingleElectron'+era]
      #if 'eett' in chn:
      #    data_samples+=['DoubleEG'+era]
      if 'emtt' in chn:
          data_samples+=['MuonEG'+era]
      if 'tttt' in chn or 'mttt' in chn or "ettt" in chn or "emtt" in chn or "eett" in chn or "mmtt" in chn:
          data_samples+=['Tau'+era]

  DATAFILELIST="./filelists/Sep18_2016_Data_102X"

  for sa in data_samples:
      JOB='%s_2016' % (sa)
      DATAFILELIST_ = DATAFILELIST
      user='guttley'
      prefix='Sep18_Data_102X_2016'
      if 'SingleElectron' in sa: 
        DATAFILELIST_ = "./filelists/Sep18_2016_etauRedo_Data_102X"
        prefix='Sep18_Data_102X_2016_etauRedo'
        user='dwinterb'
      JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(DATAFILELIST_)s_%(sa)s.dat\",\"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/%(user)s/%(prefix)s/\",\"sequences\":{\"em\":[],\"et\":[],\"mt\":[],\"tt\":[],\"zmm\":[],\"zee\":[]}}, \"sequence\":{\"output_name\":\"%(JOB)s\",\"is_data\":true}}' "%vars());
      nfiles = sum(1 for line in open('%(DATAFILELIST_)s_%(sa)s.dat' % vars()))

      # if 'TauB' in sa:
      #   JSONPATCH= (r"'{\"job\":{\"filelist\":\"./filelists/Oct2_Data_102X_%(sa)s.dat\",\"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/mhassans/Oct2_Data_102X_2016/\",\"sequences\":{\"em\":[],\"et\":[],\"mt\":[],\"tt\":[],\"zmm\":[],\"zee\":[]}}, \"sequence\":{\"output_name\":\"%(JOB)s\",\"is_data\":true}}' "%vars());


      nperjob = 40
      
      for i in range (0,int(math.ceil(float(nfiles)/float(nperjob)))) :  
        os.system('%(JOBWRAPPER)s "./bin/HTT --cfg=%(CONFIG)s --json=%(JSONPATCH)s --offset=%(i)d --nlines=%(nperjob)d &> jobs/%(JOB)s-%(i)d.log" jobs/%(JOB)s-%(i)s.sh' %vars())
        if not parajobs and not options.condor:
            os.system('%(JOBSUBMIT)s jobs/%(JOB)s-%(i)d.sh' % vars())
        elif not parajobs and options.condor:
            outscriptname = '{}-{}.sh'.format(JOB, i)
            subfilename = '{}_{}.sub'.format(JOB, i)
            subfile = open("jobs/{}".format(subfilename), "w")
            condor_settings = CONDOR_TEMPLATE % {
              'EXE': outscriptname,
              'TASK': "{}-{}".format(JOB, i)
            }
            subfile.write(condor_settings)
            subfile.close()
            os.system('condor_submit jobs/{}'.format(subfilename))

      if parajobs: 
        os.system('%(JOBWRAPPER)s ./jobs/%(JOB)s-\$\(\(SGE_TASK_ID-1\)\).sh  jobs/parajob_%(JOB)s.sh' %vars())
        PARAJOBSUBMIT = getParaJobSubmit(int(math.ceil(float(nfiles)/float(nperjob))))
        os.system('%(PARAJOBSUBMIT)s jobs/parajob_%(JOB)s.sh' % vars())  
      file_persamp.write("%s %d\n" %(JOB, int(math.ceil(float(nfiles)/float(nperjob)))))
        

if options.proc_bkg or options.proc_all:
  central_samples = [
    'WJetsToLNu-LO',
    'WJetsToLNu-LO-ext',
    'W1JetsToLNu-LO',
    'W2JetsToLNu-LO',
    'W2JetsToLNu-LO-ext',
    'W3JetsToLNu-LO',
    'W3JetsToLNu-LO-ext',
    'W4JetsToLNu-LO',
    'W4JetsToLNu-LO-ext1',
    'W4JetsToLNu-LO-ext2',
    'WGToLNuG-ext1',
    'WGToLNuG-ext2',
    'WGToLNuG-ext3',
    'WGstarToLNuMuMu',
    'WGstarToLNuEE',
    #'DYJetsToLL',
    'DYJetsToLL-LO-ext1',
    'DYJetsToLL-LO-ext2',
    'DY1JetsToLL-LO',
    'DY2JetsToLL-LO',
    'DY3JetsToLL-LO',
    'DY4JetsToLL-LO',
    'DYJetsToLL_M-10-50-LO',
    'VVTo2L2Nu',
    'VVTo2L2Nu-ext1',
    'ZZTo2L2Q',
    'WWTo1L1Nu2Q', 
    #'WWToLNuQQ', # didn't find filelist for these two...same as above anyway
    #'WWToLNuQQ-ext',
    'WZTo2L2Q',
    'WZTo1L3Nu',
    'WZTo1L1Nu2Q',
    'WZJToLLLNu',
    ##'ZZTo4L', # use amcat one below
    'Tbar-tW',
    'T-tW',
    'Tbar-t',
    'T-t',
    'TT',
    "TTTo2L2Nu",
    "TTToHadronic",
    "TTToSemiLeptonic",
    'WJetsToLNu',
    'WJetsToLNu-ext',
    'EWKWMinus2Jets_WToLNu',
    'EWKWMinus2Jets_WToLNu-ext1',
    'EWKWMinus2Jets_WToLNu-ext2',
    'EWKWPlus2Jets_WToLNu',
    'EWKWPlus2Jets_WToLNu-ext1',
    'EWKWPlus2Jets_WToLNu-ext2',
    'EWKZ2Jets_ZToLL',
    'EWKZ2Jets_ZToLL-ext1',
    'EWKZ2Jets_ZToLL-ext2',
    'EWKZ2Jets_ZToNuNu',
    'EWKZ2Jets_ZToNuNu-ext1',
    'EWKZ2Jets_ZToNuNu-ext2',
    'ZZTo4L-amcat',
    ]

  for sa in central_samples:
      JOB='%s_2016' % (sa)
      JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(FILELIST)s_%(sa)s.dat\"}, \"sequence\":{\"output_name\":\"%(JOB)s\",%(jetuncert_string)s}}' "%vars());
      job_num=0
      for FLATJSONPATCH in flatjsons: 
        nperjob = 20
        if 'TT' in sa or 'VVTo2L2Nu' in sa: nperjob = 10
  
        if 'DY' not in sa and 'EWKZ' not in sa:
          FLATJSONPATCH = FLATJSONPATCH.replace('^scale_efake_0pi_hi^scale_efake_0pi_lo','').replace('^scale_efake_1pi_hi^scale_efake_1pi_lo','').replace('^scale_mufake_0pi_hi^scale_mufake_0pi_lo','').replace('^scale_mufake_1pi_hi^scale_mufake_1pi_lo','')
        if 'DY' not in sa and 'JetsToLNu' not in sa and 'WG' not in sa and 'EWKZ' not in sa and 'EWKW' not in sa:
          FLATJSONPATCH = FLATJSONPATCH.replace('^scale_met_hi^scale_met_lo','').replace('^res_met_hi^res_met_lo','').replace('^scale_met_njets0_hi^scale_met_njets0_lo','').replace('^res_met_njets0_hi^res_met_njets0_lo','').replace('^scale_met_njets1_hi^scale_met_njets1_lo','').replace('^res_met_njets1_hi^res_met_njets1_lo','').replace('^scale_met_njets2_hi^scale_met_njets2_lo','').replace('^res_met_njets2_hi^res_met_njets2_lo','')
        n_scales = FLATJSONPATCH.count('_lo') + FLATJSONPATCH.count('default')
        if n_scales*n_channels>=24: nperjob = 10
        if n_scales*n_channels>=48: nperjob=5
        if 'TT' in sa or 'VVTo2L2Nu' in sa:
          nperjob = 10 
          if n_scales*n_channels>24: nperjob = 5
          if n_scales*n_channels>48: nperjob = 2
        if options.jetmetuncerts and 'default' in FLATJSONPATCH: nperjob = int(math.ceil(float(nperjob)/2))

        #nperjob = int(math.ceil(float(nperjob)/max(1.,float(n_scales)*float(n_channels)/10.)))
        nfiles = sum(1 for line in open('%(FILELIST)s_%(sa)s.dat' % vars()))
        for i in range (0,int(math.ceil(float(nfiles)/float(nperjob)))) :
          os.system('%(JOBWRAPPER)s "./bin/HTT --cfg=%(CONFIG)s --json=%(JSONPATCH)s --flatjson=%(FLATJSONPATCH)s --offset=%(i)d --nlines=%(nperjob)d &> jobs/%(JOB)s-%(job_num)d.log" jobs/%(JOB)s-%(job_num)s.sh' %vars())
          if not parajobs and not options.condor:
              os.system('%(JOBSUBMIT)s jobs/%(JOB)s-%(job_num)d.sh' % vars())
          elif not parajobs and options.condor:
              outscriptname = '{}-{}.sh'.format(JOB, job_num)
              subfilename = '{}_{}.sub'.format(JOB, job_num)
              subfile = open("jobs/{}".format(subfilename), "w")
              condor_settings = CONDOR_TEMPLATE % {
                'EXE': outscriptname,
                'TASK': "{}-{}".format(JOB, job_num)
              }
              subfile.write(condor_settings)
              subfile.close()
              os.system('condor_submit jobs/{}'.format(subfilename))
          job_num+=1
        file_persamp.write("%s %d\n" %(JOB, int(math.ceil(float(nfiles)/float(nperjob)))))
      if parajobs: 
        os.system('%(JOBWRAPPER)s ./jobs/%(JOB)s-\$\(\(SGE_TASK_ID-1\)\).sh  jobs/parajob_%(JOB)s.sh' %vars())
        PARAJOBSUBMIT = getParaJobSubmit(job_num)
        os.system('%(PARAJOBSUBMIT)s jobs/parajob_%(JOB)s.sh' % vars()) 
#if float(n_scales*n_channels)/100 > 1: nperjob = int(math.ceil(nperjob/(float(n_scales*n_channels)/100)))  

