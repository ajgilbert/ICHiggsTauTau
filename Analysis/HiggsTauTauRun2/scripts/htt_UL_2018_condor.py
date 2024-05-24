#!/usr/bin/env python

# importing libraries
import sys
from optparse import OptionParser
import os
import math
import json

#-----------------------------------------
os.system("voms-proxy-init --voms cms --valid 192:00 --out ~/cms.proxy")

JOBWRAPPER      = './scripts/generate_job.sh'
if "JOBWRAPPER" in os.environ:      JOBWRAPPER      = os.environ["JOBWRAPPER"]
print("Using job-wrapper:    " + JOBWRAPPER)

CONDOR_TEMPLATE = """executable = ./jobs/%(EXE)s
Proxy_path =/afs/cern.ch/user/a/adow/private/x509up
arguments = $(ProcId) $(Proxy_path)
output                = ./jobs/%(TASK)s.$(ClusterId).$(ProcId).out
error                 = ./jobs/%(TASK)s.$(ClusterId).$(ProcId).err
log                   = ./jobs/%(TASK)s.$(ClusterId).log
requirements = (OpSysAndVer =?= "SLCern6")
+JobFlavour     = "longlunch"
queue
"""
#-----------------------------------------

# parsing options
#-----------------------------------------
def split_callback(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()
parser.add_option("--wrapper", dest="wrapper",
                  help="Specify the job-wrapper script. The current wrapper is '%(JOBWRAPPER)s'."
                  " Using the --wrapper option overrides both the default and the environment variable. " % vars())

parser.add_option("--submit", dest="submit",
                  help="submit jobs", action='store_true')

parser.add_option("--data", dest="proc_data", action='store_true', default=False,
                  help="Process data samples")

parser.add_option("--calc_lumi", dest="calc_lumi", action='store_true', default=False,
                  help="Run on data and only write out lumi mask jsons")

parser.add_option("--bkg", dest="proc_bkg", action='store_true', default=False,
                  help="Process background mc samples")

parser.add_option("--all", dest="proc_all", action='store_true', default=False,
                  help="Process all samples")

parser.add_option("--no_json", dest="no_json", action='store_true', default=False,
                  help="Do not read the channels to process from the json to decide which datasets to run on")

parser.add_option("--scales", dest="scales", type='string', default='default',
                  help="List of systematic shifts to process")

parser.add_option("--config", dest="config", type='string', default='scripts/config_UL_2018.json',
                  help="Config file")

parser.add_option("--list_backup", dest="slbackupname", type='string', default='prevlist',
                  help="Name you want to give to the previous files_per_samples file, in case you're resubmitting a subset of jobs")

parser.add_option("--effective_events", action='store_true', default=False,
                  help="Generating effective event numbers for processing all samples")

parser.add_option("--jetmetuncerts", dest="jetmetuncerts", action='store_true', default=False,
                  help="Do JES, JER, and MET uncertainties")
#-----------------------------------------
(options, args) = parser.parse_args()
if options.wrapper: JOBWRAPPER=options.wrapper

jetuncert_string=''
if options.jetmetuncerts:
  jetuncert_string='\\"do_jetmet_uncerts\\":true'
else:
  jetuncert_string='\\"do_jetmet_uncerts\\":false'

def CreateCondorScript(name, queue=None):
  name_sh = name
  name = name.replace('.sh', '.sub')
  if os.path.exists(name): os.system('rm %(name)s' % vars())
  output = name.replace('.sub','_output.log')
  error = name.replace('.sub','_error.log')
  os.system('echo "Universe = vanilla" >> %(name)s' % vars())
  os.system('echo "executable = %(name_sh)s" >> %(name)s' % vars())
  os.system('echo "arguments = \\"\$(CLUSTER) \$(PROCESS)\\"" >> %(name)s' % vars())
  os.system('echo "output = %(output)s" >> %(name)s' % vars())
  os.system('echo "error = %(error)s" >> %(name)s' % vars())
  os.system('echo "+MaxRuntime = 10799" >> %(name)s' % vars())
  os.system('echo "+OS = \\"EL7\\"">> %(name)s' % vars())
  if queue: os.system('echo "queue %(queue)i" >> %(name)s' % vars())
  else: os.system('echo "queue" >> %(name)s' % vars())
  print("Created condor script:",name)

BACKUPNAME = options.slbackupname
scales = options.scales
no_json = options.no_json

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

CONFIG='scripts/config_UL_2018.json'
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

n_channels=4
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


scale = int(math.ceil(float(n_scales*n_channels)/30))
#scale = int(math.ceil(float(n_scales*n_channels)/8)) # change back later!
if scale < 1: scale = 1

total = float(len(flatjsonlistdysig))
flatjsons = []
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

FILELIST='filelists/Mar31_2018_MC_124X'

signal_mc = [ ]
signal_vh = [ ]
signal_mc_ww = [ ]

if os.path.isfile("./jobs/files_per_sample_2018.txt"):
    os.system("mv ./jobs/files_per_sample_2018.txt ./jobs/files_per_sample_2018-%(BACKUPNAME)s.txt"%vars())

file_persamp = open("./jobs/files_per_sample_2018.txt", "w")


if options.proc_data or options.proc_all or options.calc_lumi:
    if not no_json:
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
    else:
        channels=['mt','et','tt','em','zmm','zee']

if options.proc_data or options.proc_all or options.calc_lumi:

    data_samples = []
    data_eras = ['A','B','C','D']
  
    for chn in channels:
        for era in data_eras:
            if 'mt' in chn or 'zmm' in chn:
                if 'SingleMuon'+era not in data_samples: data_samples+=['SingleMuon'+era]
            if 'et' in chn or 'zee' in chn:
                if 'EGamma'+era not in data_samples: data_samples+=['EGamma'+era]
            if 'em' in chn:
                data_samples+=['MuonEG'+era]
            if 'tt' in chn:
                data_samples+=['Tau'+era]

    if ('mt' in channels or 'et' in channels) and 'tt' not in channels:
      for era in data_eras:
        data_samples+=['Tau'+era]

    if options.effective_events: # this is just used for the effective events case
        for era in data_eras:
          if 'SingleMuon'+era not in data_samples: data_samples+=['SingleMuon'+era]
          if 'EGamma'+era not in data_samples: data_samples+=['EGamma'+era]
          if 'Tau'+era not in data_samples: data_samples+=['Tau'+era]
          if 'MuonEG'+era not in data_samples: data_samples+=['MuonEG'+era]
          if 'DoubleMuon'+era not in data_samples: data_samples+=['DoubleMuon'+era]

    DATAFILELIST="./filelists/Mar31_2018_Data_124X"

    for sa in data_samples:
        JOB='%s_2018' % (sa)
        user='dwinterb'
        prefix='Mar31_Data_124X_2018'
        JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(DATAFILELIST)s_%(sa)s.dat\",\"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/%(user)s/%(prefix)s/\",\"sequences\":{\"em\":[],\"et\":[],\"mt\":[],\"tt\":[],\"zmm\":[],\"zee\":[]}}, \"sequence\":{\"output_name\":\"%(JOB)s\",\"is_data\":true}}' "%vars());
        nfiles = sum(1 for line in open('%(DATAFILELIST)s_%(sa)s.dat' % vars()))
        nperjob = 30
    
        if "EGammaD" in sa: nperjob = 50
        # elif "TauD" in sa: nperjob = 45
        
        for i in range (0,int(math.ceil(float(nfiles)/float(nperjob)))) :  
            os.system('%(JOBWRAPPER)s "./bin/HTT --cfg=%(CONFIG)s --json=%(JSONPATCH)s --offset=%(i)d --nlines=%(nperjob)d &> jobs/%(JOB)s-%(i)d.log" jobs/%(JOB)s-%(i)s.sh' %vars())
        os.system('%(JOBWRAPPER)s "./jobs/%(JOB)s-\$2.sh"  jobs/parajob_%(JOB)s.sh' %vars())
        njobs=int(math.ceil(float(nfiles)/float(nperjob)))
        CreateCondorScript('jobs/parajob_%(JOB)s.sh' % vars(), queue=njobs)
        if options.submit: os.system('condor_submit jobs/parajob_%(JOB)s.sub' % vars())
        file_persamp.write("%s %d\n" %(JOB, int(math.ceil(float(nfiles)/float(nperjob)))))


if options.proc_bkg or options.proc_all:
    central_samples = [
    # Drell-Yan LO
    'DY1JetsToLL-LO',
    'DY2JetsToLL-LO',
    'DY3JetsToLL-LO',
    'DY4JetsToLL-LO',
    'DYJetsToLL-LO',
    'DYJetsToLL-LO-ext1',

    # Low mass Drell Yan LO
    'DYJetsToLL_M-10to50-LO',
 
    # Drell-Yan NLO
    'DYJetsToLL-NLO',
    'DYJetsToLL_0J-NLO',
    'DYJetsToLL_1J-NLO',
    'DYJetsToLL_2J-NLO',
   
    # Electroweak W and Z
    'EWKWMinus2Jets_WToLNu',
    'EWKWPlus2Jets_WToLNu',
    'EWKZ2Jets_ZToLL',
   
    # W + Jets L0
    'WJetsToLNu-LO',
    'W1JetsToLNu-LO',
    'W2JetsToLNu-LO',
    'W3JetsToLNu-LO',
    'W4JetsToLNu-LO',
   
#    # W + Jets NLO
#    'WJetsToLNu_0J-NLO',
#    'WJetsToLNu_1J-NLO',
#    'WJetsToLNu_2J-NLO',
#    #'WJetsToLNu-NLO',
  
   # ttbar
    'TTTo2L2Nu',
    'TTToHadronic',
    'TTToSemiLeptonic',
  
    # Split diboson (Missing Files: WZTo1L3Nu, WZTo2L2Q)
    'WZTo1L1Nu2Q',
    'WZTo3LNu',
    'WWTo1L1Nu2Q',
    'WWTo2L2Nu',
    'WZTo1L3Nu',
    'WZTo2Q2L',
    'ZZTo2L2Nu',
    'ZZTo4L',
 
#    # Triboson
#    'WWZ',
#    'WWZ-ext1',
#    'WZZ',
#    'WZZ-ext1',
#    'WWW',
#    'WWW-ext1',
#    'ZZZ',
#    'ZZZ-ext1',
 
    # Other backgrounds
    'WGToLNuG',
    'Tbar-t',
    'Tbar-tW',
    'T-t',
    'T-tW',
 
#    # SM Higgs
#    'GluGluHToTauTau_M125',
#    'VBFHToTauTau_M125',
#    'WminusHToTauTau_M125',
#    'WplusHToTauTau_M125',
#    'ZHToTauTau_M125',
#    'ZHToTauTau_M125-ext1',
#    'ttHToTauTau_M125'
    ]

    if options.effective_events: # this is just used for the effective events case
      central_samples += [
                          'VBF_HToZZTo4L_M125',
                          'GluGlu_HToZZTo4L_M125',
                          'GluGluToContinToZZTo2e2mu',
                          'GluGluToContinToZZTo2e2tau',
                          'GluGluToContinToZZTo2mu2tau',
                          'GluGluToContinToZZTo4e',
                          'GluGluToContinToZZTo4mu',
                          'GluGluToContinToZZTo4tau',
                          ]


    for sa in central_samples:
        JOB='%s_2018' % (sa)
        user='dwinterb'
        PREFIX='Mar31_MC_124X_2018'
        JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(FILELIST)s_%(sa)s.dat\", \"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/%(user)s/%(PREFIX)s/\"}, \"sequence\":{\"output_name\":\"%(JOB)s\",%(jetuncert_string)s}}' "%vars());
        
        job_num=0
        for FLATJSONPATCH in flatjsons:
            #nperjob = 40
            nperjob=10
            if 'DY' not in sa and 'EWKZ' not in sa:
                FLATJSONPATCH = FLATJSONPATCH.replace('^scale_efake_0pi_hi^scale_efake_0pi_lo','').replace('^scale_efake_1pi_hi^scale_efake_1pi_lo','').replace('^scale_mufake_0pi_hi^scale_mufake_0pi_lo','').replace('^scale_mufake_1pi_hi^scale_mufake_1pi_lo','')
            if 'DY' not in sa and 'JetsToLNu' not in sa and 'WG' not in sa and 'EWKZ' not in sa and 'EWKW' not in sa:
                FLATJSONPATCH = FLATJSONPATCH.replace('^scale_met_hi^scale_met_lo','').replace('^res_met_hi^res_met_lo','').replace('^scale_met_njets0_hi^scale_met_njets0_lo','').replace('^res_met_njets0_hi^res_met_njets0_lo','').replace('^scale_met_njets1_hi^scale_met_njets1_lo','').replace('^res_met_njets1_hi^res_met_njets1_lo','').replace('^scale_met_njets2_hi^scale_met_njets2_lo','').replace('^res_met_njets2_hi^res_met_njets2_lo','')
            #else:
            #    FLATJSONPATCH = FLATJSONPATCH.replace('^met_uncl_hi^met_uncl_lo','')
            if FLATJSONPATCH == 'job:sequences:all:^^' or FLATJSONPATCH == 'job:sequences:all:': continue
            n_scales = FLATJSONPATCH.count('_lo')*2 + FLATJSONPATCH.count('default')
            if n_scales*n_channels>=24: nperjob = 5
            if n_scales*n_channels>=48: nperjob=2

            if options.jetmetuncerts and 'default' in FLATJSONPATCH: nperjob = int(math.ceil(float(nperjob)/2))

            if 'TTTo' in sa: nperjob = int(math.ceil(float(nperjob)/2)) 
            if sa == 'TTTo2L2Nu': nperjob = int(math.ceil(float(nperjob)/2))
            nfiles = sum(1 for line in open('%(FILELIST)s_%(sa)s.dat' % vars()))
            for i in range (0,int(math.ceil(float(nfiles)/float(nperjob)))) :
                os.system('%(JOBWRAPPER)s "./bin/HTT --cfg=%(CONFIG)s --json=%(JSONPATCH)s --flatjson=%(FLATJSONPATCH)s --offset=%(i)d --nlines=%(nperjob)d &> jobs/%(JOB)s-%(job_num)d.log" jobs/%(JOB)s-%(job_num)s.sh' %vars())
                job_num+=1
            file_persamp.write("%s %d\n" %(JOB, int(math.ceil(float(nfiles)/float(nperjob)))))
        os.system('%(JOBWRAPPER)s "./jobs/%(JOB)s-\$2.sh"  jobs/parajob_%(JOB)s.sh' %vars())
        CreateCondorScript('jobs/parajob_%(JOB)s.sh' % vars(), queue=job_num)
        if options.submit: os.system('condor_submit jobs/parajob_%(JOB)s.sub' % vars())



