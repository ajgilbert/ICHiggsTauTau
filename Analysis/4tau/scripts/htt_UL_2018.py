#!/usr/bin/env python

# python scripts/htt_UL_2018.py --bkg --data --jetmetuncerts --scales="default,scale_t_0pi,scale_t_1pi,scale_t_3prong,scale_t_3prong1pi0,scale_efake_0pi,scale_efake_1pi,scale_mufake_0pi,scale_mufake_1pi,scale_e" --submit='./scripts/submit_ic_batch_job.sh "hep.q -l h_rt=0:180:0 -l h_vmem=24G"' --parajobs

# python scripts/htt_UL_2018.py --bkg --signal --batch --parajobs

# importing libraries
import sys
from optparse import OptionParser
import os
import math
import json

# setting up job stuff and condor template
#-----------------------------------------
os.system("voms-proxy-init --voms cms --valid 96:00 --out ~/cms.proxy")

JOBWRAPPER      = './scripts/generate_job.sh'
JOBSUBMIT       = 'true'
if "JOBWRAPPER" in os.environ:      JOBWRAPPER      = os.environ["JOBWRAPPER"]
if "JOBSUBMIT"  in os.environ:      JOBSUBMIT       = os.environ["JOBSUBMIT"]
print "Using job-wrapper:    " + JOBWRAPPER
print "Using job-submission: " + JOBSUBMIT

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
                  help="Specify the job-submission method. The current method is '%(JOBSUBMIT)s'"
                  " Using the --submit option overrides both the default and the environment variable. " % vars())

parser.add_option("--data", dest="proc_data", action='store_true', default=False,
                  help="Process data samples")

parser.add_option("--calc_lumi", dest="calc_lumi", action='store_true', default=False,
                  help="Run on data and only write out lumi mask jsons")

parser.add_option("--bkg", dest="proc_bkg", action='store_true', default=False,
                  help="Process background mc samples")

parser.add_option("--sig", dest="proc_sig", action='store_true', default=False,
                  help="Process signal mc samples")

parser.add_option("--all", dest="proc_all", action='store_true', default=False,
                  help="Process all samples")

parser.add_option("--no_json", dest="no_json", action='store_true', default=False,
                  help="Do not read the channels to process from the json to decide which datasets to run on")

parser.add_option("--scales", dest="scales", type='string', default='default',
                  help="List of systematic shifts to process")

parser.add_option("--parajobs", dest="parajobs", action='store_true', default=False,
                  help="Submit jobs parametrically")

parser.add_option("--config", dest="config", type='string', default='scripts/config_UL_2018.json',
                  help="Config file")

parser.add_option("--list_backup", dest="slbackupname", type='string', default='prevlist',
                  help="Name you want to give to the previous files_per_samples file, in case you're resubmitting a subset of jobs")

parser.add_option("--condor", action='store_true', default=False,
                  help="Submit jobs to condor (for lxplus)")

parser.add_option("--jetmetuncerts", dest="jetmetuncerts", action='store_true', default=False,
                  help="Do JES, JER, and MET uncertainties")

parser.add_option("--batch", dest="batch", action='store_true', default=False,
                  help="Submit to imperial batch")

#-----------------------------------------
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
    return sub_opts

BACKUPNAME = options.slbackupname
scales = options.scales
no_json = options.no_json
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


FILELIST='filelists/Jan10_2018_MC_106X'

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
        channels=["tttt","mttt","ettt","emtt","eett","mmtt"]

if options.proc_data or options.proc_all or options.calc_lumi:

    data_samples = []
    data_eras = ['A','B','C','D']
    for chn in channels:
        for era in data_eras:
            if 'mttt' in chn:
                data_samples+=['SingleMuon'+era]
            if 'mmtt' in chn:
                data_samples+=['SingleMuon'+era, 'DoubleMuon'+era]
            if 'ettt' in chn or 'eett' in chn:
                data_samples+=['EGamma'+era]
            if 'emtt' in chn:
                data_samples+=['MuonEG'+era]
            if 'tttt' in chn or 'mttt' in chn or "ettt" in chn or "emtt" in chn or "eett" in chn or "mmtt" in chn or "ttt" in chn:
                data_samples+=['Tau'+era]
    data_samples = list(set(data_samples))

    DATAFILELIST="./filelists/Apr22_2018_Data_106X"

    if options.calc_lumi:
        for sa in data_samples:
            DATAFILELIST="./filelists/Apr22_2018_Data_106X"
            JOB='%s_2018' % (sa)
            JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(DATAFILELIST)s_%(sa)s.dat\",\"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/guttley/Apr22_Data_106X_2018/\",\"sequences\":{\"em\":[],\"et\":[],\"mt\":[],\"tt\":[]}}, \"sequence\":{\"output_name\":\"%(JOB)s\",\"is_data\":true,\"lumi_mask_only\":true}}' "%vars());
            if sa == "EGammaD":
                DATAFILELIST="./filelists/Feb16_2018_Data_106X"
                JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(DATAFILELIST)s_%(sa)s.dat\",\"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/guttley/Feb16_Data_106X_2018/\",\"sequences\":{\"em\":[],\"et\":[],\"mt\":[],\"tt\":[]}}, \"sequence\":{\"output_name\":\"%(JOB)s\",\"is_data\":true,\"lumi_mask_only\":true}}' "%vars()); 
	        nfiles = sum(1 for line in open('%(DATAFILELIST)s_%(sa)s.dat' % vars()))
            nperjob = 500 
            if "TauC" in sa: nperjob = 252
            elif "TauD" in sa: nperjob = 271
            for i in range (0,int(math.ceil(float(nfiles)/float(nperjob)))):
                os.system('%(JOBWRAPPER)s "./bin/HTT --cfg=%(CONFIG)s --json=%(JSONPATCH)s --offset=%(i)d --nlines=%(nperjob)d &> jobs/%(JOB)s-%(i)d.log" jobs/%(JOB)s-%(i)s.sh' %vars())
                if not parajobs:
                    os.system('%(JOBSUBMIT)s jobs/%(JOB)s-%(i)d.sh' % vars())
            if parajobs:
                os.system('%(JOBWRAPPER)s ./jobs/%(JOB)s-\$\(\(SGE_TASK_ID-1\)\).sh  jobs/parajob_%(JOB)s.sh' %vars())
                PARAJOBSUBMIT = getParaJobSubmit(int(math.ceil(float(nfiles)/float(nperjob))))
                os.system('%(PARAJOBSUBMIT)s jobs/parajob_%(JOB)s.sh' % vars())
            file_persamp.write("%s %d\n" %(JOB, int(math.ceil(float(nfiles)/float(nperjob)))))

    else:
        for sa in data_samples:
            DATAFILELIST="./filelists/Apr22_2018_Data_106X"
            JOB='%s_2018' % (sa)
            JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(DATAFILELIST)s_%(sa)s.dat\",\"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/guttley/Apr22_Data_106X_2018/\",\"sequences\":{\"em\":[],\"et\":[],\"mt\":[],\"tt\":[],\"zmm\":[],\"zee\":[]}}, \"sequence\":{\"output_name\":\"%(JOB)s\",\"is_data\":true}}' "%vars());
            if sa == "EGammaD":
                DATAFILELIST="./filelists/Feb16_2018_Data_106X"   
                JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(DATAFILELIST)s_%(sa)s.dat\",\"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/guttley/Feb16_Data_106X_2018/\",\"sequences\":{\"em\":[],\"et\":[],\"mt\":[],\"tt\":[],\"zmm\":[],\"zee\":[]}}, \"sequence\":{\"output_name\":\"%(JOB)s\",\"is_data\":true}}' "%vars()); 
            nfiles = sum(1 for line in open('%(DATAFILELIST)s_%(sa)s.dat' % vars()))
            nperjob = 30
      
            if "EGammaD" in sa: nperjob = 50
            # elif "TauD" in sa: nperjob = 45
            
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
    # Drell-Yan LO
    'DY1JetsToLL-LO',
    'DY2JetsToLL-LO',
    'DY3JetsToLL-LO',
    'DY4JetsToLL-LO',
    'DYJetsToLL-LO',
	  # Low mass Drell Yan LO
    'DYJetsToLL_M-10to50-LO',
  	#'DY1JetsToLL_M-10to50-LO',
  	#'DY2JetsToLL_M-10to50-LO',
 	  #'DY3JetsToLL_M-10to50-LO',
 	  #'DY4JetsToLL_M-10to50-LO',
	  # Drell-Yan NLO
 	  #'DYJetsToLL_0J-NLO',
    #'DYJetsToLL_1J-NLO',
	  #'DYJetsToLL_2J-NLO',
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
	  # W + Jets NLO
	  #'WJetsToLNu_0J-NLO',
	  #'WJetsToLNu_1J-NLO',
    #'WJetsToLNu_2J-NLO',
    #'WJetsToLNu-NLO',
	  # ttbar
 	  'TTTo2L2Nu',
	  'TTToHadronic',
 	  'TTToSemiLeptonic',
	  # Split diboson (Missing Files: WZTo1L3Nu, WZTo2L2Q)
 	  'WZTo1L1Nu2Q',
 	  'WZTo3LNu',
 	  'WWTo1L1Nu2Q',
 	  'WWTo2L2Nu',
 	  'ZZTo2L2Nu',
	  'ZZTo4L',
	  # Inclusive
 	  'WW',
 	  'WZ',
 	  'ZZ',
 	  # Triboson
  	'WWZ',
 	  'WWZ-ext1',
 	  'WZZ',
 	  'WZZ-ext1',
  	'WWW',
 	  'WWW-ext1',
 	  'ZZZ',
 	  'ZZZ-ext1',
	  # Other backgrounds
 	  'WGToLNuG',
 	  'Tbar-t',
 	  'Tbar-tW',
 	  'T-t',
 	  'T-tW',
	  # SM Higgs
 	 # 'GluGluHToTauTau_M125',
 	 # 'VBFHToTauTau_M125',
 	 # 'WminusHToTauTau_M125',
 	 # 'WplusHToTauTau_M125',
 	 # 'ttHToTauTau_M125',
    ]


    for sa in central_samples:
        FILELIST='filelists/Jan10_2018_MC_106X'
        JOB='%s_2018' % (sa)
        JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(FILELIST)s_%(sa)s.dat\", \"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/guttley/Jan10_MC_106X_2018/\"}, \"sequence\":{\"output_name\":\"%(JOB)s\",%(jetuncert_string)s}}' "%vars());
        if sa == "TTToSemiLeptonic":
            FILELIST='filelists/Feb16_2018_MC_106X'
            JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(FILELIST)s_%(sa)s.dat\", \"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/guttley/Feb16_MC_106X_2018/\"}, \"sequence\":{\"output_name\":\"%(JOB)s\",%(jetuncert_string)s}}' "%vars());
        job_num=0
        for FLATJSONPATCH in flatjsons:
            #nperjob = 40
            nperjob=20
            if 'DY' not in sa and 'EWKZ' not in sa:
                FLATJSONPATCH = FLATJSONPATCH.replace('^scale_efake_0pi_hi^scale_efake_0pi_lo','').replace('^scale_efake_1pi_hi^scale_efake_1pi_lo','').replace('^scale_mufake_0pi_hi^scale_mufake_0pi_lo','').replace('^scale_mufake_1pi_hi^scale_mufake_1pi_lo','')
            if 'DY' not in sa and 'JetsToLNu' not in sa and 'WG' not in sa and 'EWKZ' not in sa and 'EWKW' not in sa:
                FLATJSONPATCH = FLATJSONPATCH.replace('^scale_met_hi^scale_met_lo','').replace('^res_met_hi^res_met_lo','').replace('^scale_met_njets0_hi^scale_met_njets0_lo','').replace('^res_met_njets0_hi^res_met_njets0_lo','').replace('^scale_met_njets1_hi^scale_met_njets1_lo','').replace('^res_met_njets1_hi^res_met_njets1_lo','').replace('^scale_met_njets2_hi^scale_met_njets2_lo','').replace('^res_met_njets2_hi^res_met_njets2_lo','')
            #else:
            #    FLATJSONPATCH = FLATJSONPATCH.replace('^met_uncl_hi^met_uncl_lo','')
            if FLATJSONPATCH == 'job:sequences:all:^^' or FLATJSONPATCH == 'job:sequences:all:': continue
            n_scales = FLATJSONPATCH.count('_lo')*2 + FLATJSONPATCH.count('default')
            if n_scales*n_channels>=24: nperjob = 10
            if n_scales*n_channels>=48: nperjob=5

            if options.jetmetuncerts and 'default' in FLATJSONPATCH: nperjob = int(math.ceil(float(nperjob)/2))

            #if 'TTTo' in sa: nperjob = int(math.ceil(float(nperjob)/2)) 
            #if 'TTTo' in sa: nperjob = int(math.ceil(float(nperjob)/2)) 
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


if options.proc_sig:

    signal_mc = [
    #"phi200A100To4Tau",
    #"phi200A200To4Tau",
    "phi100A100To4Tau",
    "phi100A150To4Tau",
    "phi100A60To4Tau",
    "phi200A100To4Tau",
    "phi200A150To4Tau",
    "phi200A20To4Tau",
    "phi200A60To4Tau",
    "phi300A100To4Tau",
    "phi300A150To4Tau",
    "phi300A20To4Tau",
    "phi300A60To4Tau"
    ]

    for sa in signal_mc:

        print sa
        #SIG_DIR = 'June03_Signal_106X_2018'
        #SIG_FILELIST = "filelists/June03_2018_MC_106X"
        SIG_DIR = 'Jun15_Signal_106X_2018'
        SIG_FILELIST = "filelists/Jun15_2018_MC_106X"
        user='guttley'

        JOB='%s_2018' % (sa)
        JSONPATCH= (r"'{\"job\":{\"filelist\":\"%(SIG_FILELIST)s_%(sa)s.dat\",\"file_prefix\":\"root://gfe02.grid.hep.ph.ic.ac.uk:1097//store/user/%(user)s/%(SIG_DIR)s/\"}, \"sequence\":{\"output_name\":\"%(JOB)s\",%(jetuncert_string)s}}' "%vars());

        job_num=0
        for FLATJSONPATCH in flatjsons:
            FLATJSONPATCH = FLATJSONPATCH.replace('^scale_efake_0pi_hi^scale_efake_0pi_lo','').replace('^scale_efake_1pi_hi^scale_efake_1pi_lo','').replace('^scale_mufake_0pi_hi^scale_mufake_0pi_lo','').replace('^scale_mufake_1pi_hi^scale_mufake_1pi_lo','')
            FLATJSONPATCH = FLATJSONPATCH.replace('^met_uncl_hi^met_uncl_lo','')
            if FLATJSONPATCH == 'job:sequences:all:^^' or FLATJSONPATCH == 'job:sequences:all:': continue
            #print '%(SIG_FILELIST)s_%(sa)s.dat' %vars(), os.path.exists('%(SIG_FILELIST)s_%(sa)s.dat' %vars())
            if os.path.exists('%(SIG_FILELIST)s_%(sa)s.dat' %vars()):
                nfiles = sum(1 for line in open('%(SIG_FILELIST)s_%(sa)s.dat' % vars()))
                nperjob = 5
                n_scales = FLATJSONPATCH.count('_lo')*2 + FLATJSONPATCH.count('default')
                if n_scales*n_channels>=24: nperjob = 2
                if n_scales*n_channels>=48: nperjob=1

                if options.jetmetuncerts and 'default' in FLATJSONPATCH: nperjob = int(math.ceil(float(nperjob)/2))
 
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
                        # print('condor_submit jobs/{}'.format(subfilename))
                    job_num+=1
                file_persamp.write("%s %d\n" %(JOB, int(math.ceil(float(nfiles)/float(nperjob)))))
        if parajobs:
            os.system('%(JOBWRAPPER)s ./jobs/%(JOB)s-\$\(\(SGE_TASK_ID-1\)\).sh  jobs/parajob_%(JOB)s.sh' %vars())
            PARAJOBSUBMIT = getParaJobSubmit(job_num)
            os.system('%(PARAJOBSUBMIT)s jobs/parajob_%(JOB)s.sh' % vars())

