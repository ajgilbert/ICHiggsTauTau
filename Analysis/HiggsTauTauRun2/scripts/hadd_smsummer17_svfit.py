#!/usr/bin/env python

import sys
import os
from optparse import OptionParser
import math
import fnmatch
import glob

parser = OptionParser()

parser.add_option("--folder", dest = "folder",
                  help="Specify folder that contains the output to be hadded")
parser.add_option("--batch", dest= "batch", default=False, action='store_true',
                  help="Submit as batch jobs")

def list_paths(path):
    directories = []
    for item in os.listdir(path):
      if os.path.isdir(os.path.join(path, item)):
        directories.append(item)
    return directories


(options,args) = parser.parse_args()

if not options.folder:
  parser.error('No folder specified')

outputf = options.folder
batch = options.batch

JOBWRAPPER      = './scripts/generate_job.sh'
JOBSUBMIT       = './scripts/submit_ic_batch_job.sh "hep.q -l h_rt=0:180:0"'

sample_list = [
    "DY1JetsToLL-LO",
    "DY1JetsToLL-LO-ext",
    "DY2JetsToLL-LO",
    "DY2JetsToLL-LO-ext",
    "DY3JetsToLL-LO",
    "DY3JetsToLL-LO-ext",
    "DY4JetsToLL-LO",
    "DYJetsToLL",
    "DYJetsToLL-LO", # buggy PU
    "DYJetsToLL-LO-ext1", # buggy PU
    "DYJetsToLL-ext",
    "DYJetsToLL_M-10-50-LO",
    "DYJetsToLL_M-10-50-LO-ext1",
    "EWKWMinus2Jets",
    "EWKWPlus2Jets",
    "EWKZ2Jets",
    "EmbeddingElElB",
    "EmbeddingElElC",
    "EmbeddingElElD",
    "EmbeddingElElE",
    "EmbeddingElElF",
    "EmbeddingElMuB",
    "EmbeddingElMuC",
    "EmbeddingElMuD",
    "EmbeddingElMuE",
    "EmbeddingElMuF",
    "EmbeddingElTauB",
    "EmbeddingElTauC",
    "EmbeddingElTauD",
    "EmbeddingElTauE",
    "EmbeddingElTauF",
    "EmbeddingMuMuB",
    "EmbeddingMuMuC",
    "EmbeddingMuMuD",
    "EmbeddingMuMuE",
    "EmbeddingMuMuF",
    "EmbeddingMuTauB",
    "EmbeddingMuTauC",
    "EmbeddingMuTauD",
    "EmbeddingMuTauE",
    "EmbeddingMuTauF",
    "EmbeddingTauTauB",
    "EmbeddingTauTauC",
    "EmbeddingTauTauD",
    "EmbeddingTauTauE",
    "EmbeddingTauTauF",
    "GluGluHToTauTauUncorrelatedDecay",
    "GluGluHToTauTauUncorrelatedDecay_Filtered",
    "GluGluToHToTauTauPlusTwoJets_M125_amcatnloFXFX",
    "GluGluToHToTauTau_M-125-nospinner",
    "GluGluToHToTauTau_M-125-nospinner-filter",
    "GluGluToHToTauTau_M125_amcatnloFXFX",
    "GluGluToMaxmixHToTauTauPlusTwoJets_M125_amcatnloFXFX",
    "GluGluToPseudoscalarHToTauTauPlusTwoJets_M125_amcatnloFXFX",
    "JJH0MToTauTauPlusOneJets",
    "JJH0MToTauTauPlusOneJets_Filtered",
    "JJH0MToTauTauPlusTwoJets",
    "JJH0MToTauTauPlusTwoJets_Filtered",
    "JJH0MToTauTauPlusZeroJets",
    "JJH0MToTauTauPlusZeroJets_Filtered",
    "JJH0Mf05ph0ToTauTauPlusOneJets",
    "JJH0Mf05ph0ToTauTauPlusOneJets_Filtered",
    "JJH0Mf05ph0ToTauTauPlusTwoJets",
    "JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered",
    "JJH0Mf05ph0ToTauTauPlusZeroJets",
    "JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered",
    "JJH0PMToTauTauPlusOneJets",
    "JJH0PMToTauTauPlusOneJets_Filtered",
    "JJH0PMToTauTauPlusTwoJets",
    "JJH0PMToTauTauPlusTwoJets_Filtered",
    "JJH0PMToTauTauPlusZeroJets",
    "JJH0PMToTauTauPlusZeroJets_Filtered",
    "JJHiggs0MToTauTau",
    "JJHiggs0Mf05ph0ToTauTau",
    "JJHiggs0PMToTauTau",
    "MuonEGB",
    "MuonEGC",
    "MuonEGD",
    "MuonEGE",
    "MuonEGF",
    "SingleElectronB",
    "SingleElectronC",
    "SingleElectronD",
    "SingleElectronE",
    "SingleElectronF",
    "SingleMuonB",
    "SingleMuonC",
    "SingleMuonD",
    "SingleMuonE",
    "SingleMuonF",
    "TauB",
    "TauC",
    "TauD",
    "TauE",
    "TauF",
    "T-t",
    "T-tW",
    "TTTo2L2Nu",
    "TTToHadronic",
    "TTToSemiLeptonic",
    "Tbar-t",
    "Tbar-tW",
    "VBFHToTauTauUncorrelatedDecay",
    "VBFHToTauTauUncorrelatedDecay_Filtered",
    "VBFHToTauTau_M-125-nospinner",
    "VBFHToTauTau_M-125-nospinner-filter",
    "VBFHiggs0L1ToTauTau",
    "VBFHiggs0L1ZgToTauTau",
    "VBFHiggs0L1Zgf05ph0ToTauTau",
    "VBFHiggs0L1f05ph0ToTauTau",
    "VBFHiggs0MToTauTau",
    "VBFHiggs0Mf05ph0ToTauTau",
    "VBFHiggs0PHToTauTau",
    "VBFHiggs0PHf05ph0ToTauTau",
    "VBFHiggs0PMToTauTau",
    "VVTo2L2Nu",
    "W1JetsToLNu-LO",
    "W2JetsToLNu-LO",
    "W3JetsToLNu-LO", # buggy PU
    "W4JetsToLNu-LO",
    "WGToLNuG",
    "WHiggs0L1ToTauTau",
    "WHiggs0L1f05ph0ToTauTau",
    "WHiggs0MToTauTau",
    "WHiggs0Mf05ph0ToTauTau",
    "WHiggs0PHToTauTau",
    "WHiggs0PHf05ph0ToTauTau",
    "WHiggs0PMToTauTau",
    "WJetsToLNu-LO",
    "WJetsToLNu-LO-ext",
    "WWTo1L1Nu2Q", # buggy PU
    "WWTo2L2Nu",
    "WWToLNuQQ",
    "WZTo1L1Nu2Q",
    "WZTo1L3Nu",
    "WZTo2L2Q",
    "WZTo3LNu",
    "WminusHToTauTauUncorrelatedDecay_Filtered",
    "WplusHToTauTauUncorrelatedDecay",
    "WplusHToTauTauUncorrelatedDecay_Filtered",
    "ZHToTauTauUncorrelatedDecay_Filtered",
    "ZHiggs0L1ToTauTau",
    "ZHiggs0L1ZgToTauTau",
    "ZHiggs0L1Zgf05ph0ToTauTau",
    "ZHiggs0L1f05ph0ToTauTau",
    "ZHiggs0MToTauTau",
    "ZHiggs0Mf05ph0ToTauTau",
    "ZHiggs0PHToTauTau",
    "ZHiggs0PHf05ph0ToTauTau",
    "ZHiggs0PMToTauTau",
    "ZZTo2L2Nu",
    "ZZTo2L2Q",
    "ZZTo4L",
    "ZZTo4L-ext",
    "ttHiggs0MToTauTau",
    "ttHiggs0Mf05ph0ToTauTau",
    "ttHiggs0PMToTauTau",
    "SUSYGluGluToBBHToTauTau_M-80",
    "SUSYGluGluToBBHToTauTau_M-90",
    "SUSYGluGluToBBHToTauTau_M-100",
    "SUSYGluGluToBBHToTauTau_M-110",
    "SUSYGluGluToBBHToTauTau_M-120",
    "SUSYGluGluToBBHToTauTau_M-125",
    "SUSYGluGluToBBHToTauTau_M-130",
    "SUSYGluGluToBBHToTauTau_M-140",
    "SUSYGluGluToBBHToTauTau_M-160",
    "SUSYGluGluToBBHToTauTau_M-180",
    "SUSYGluGluToBBHToTauTau_M-200",
    "SUSYGluGluToBBHToTauTau_M-250",
    "SUSYGluGluToBBHToTauTau_M-300",
    "SUSYGluGluToBBHToTauTau_M-350",
    "SUSYGluGluToBBHToTauTau_M-400",
    "SUSYGluGluToBBHToTauTau_M-450",
    "SUSYGluGluToBBHToTauTau_M-500",
    "SUSYGluGluToBBHToTauTau_M-600",
    "SUSYGluGluToBBHToTauTau_M-700",
    "SUSYGluGluToBBHToTauTau_M-800",
    "SUSYGluGluToBBHToTauTau_M-900",
    "SUSYGluGluToBBHToTauTau_M-1000",
    "SUSYGluGluToBBHToTauTau_M-1200",
    "SUSYGluGluToBBHToTauTau_M-1400",
    "SUSYGluGluToBBHToTauTau_M-1500",
    "SUSYGluGluToBBHToTauTau_M-1600",
    "SUSYGluGluToBBHToTauTau_M-1800",
    "SUSYGluGluToBBHToTauTau_M-2000",
    "SUSYGluGluToBBHToTauTau_M-2300",
    "SUSYGluGluToBBHToTauTau_M-2600",
    "SUSYGluGluToBBHToTauTau_M-2900",
    "SUSYGluGluToBBHToTauTau_M-3200",
    "SUSYGluGluToBBHToTauTau_M-3500",
    "SUSYGluGluToBBHToTauTau_M-80-NLO",
    "SUSYGluGluToBBHToTauTau_M-90-NLO",
    "SUSYGluGluToBBHToTauTau_M-100-NLO",
    "SUSYGluGluToBBHToTauTau_M-110-NLO",
    "SUSYGluGluToBBHToTauTau_M-120-NLO",
    "SUSYGluGluToBBHToTauTau_M-125-NLO",
    "SUSYGluGluToBBHToTauTau_M-130-NLO",
    "SUSYGluGluToBBHToTauTau_M-140-NLO",
    "SUSYGluGluToBBHToTauTau_M-160-NLO",
    "SUSYGluGluToBBHToTauTau_M-180-NLO",
    "SUSYGluGluToBBHToTauTau_M-200-NLO",
    "SUSYGluGluToBBHToTauTau_M-250-NLO",
    "SUSYGluGluToBBHToTauTau_M-300-NLO",
    "SUSYGluGluToBBHToTauTau_M-350-NLO",
    "SUSYGluGluToBBHToTauTau_M-400-NLO",
    "SUSYGluGluToBBHToTauTau_M-450-NLO",
    "SUSYGluGluToBBHToTauTau_M-500-NLO",
    "SUSYGluGluToBBHToTauTau_M-600-NLO",
    "SUSYGluGluToBBHToTauTau_M-700-NLO",
    "SUSYGluGluToBBHToTauTau_M-800-NLO",
    "SUSYGluGluToBBHToTauTau_M-900-NLO",
    "SUSYGluGluToBBHToTauTau_M-1000-NLO",
    "SUSYGluGluToBBHToTauTau_M-1200-NLO",
    "SUSYGluGluToBBHToTauTau_M-1400-NLO",
    "SUSYGluGluToBBHToTauTau_M-1600-NLO",
    "SUSYGluGluToBBHToTauTau_M-1800-NLO",
    "SUSYGluGluToBBHToTauTau_M-2000-NLO",
    "SUSYGluGluToBBHToTauTau_M-2300-NLO",
    "SUSYGluGluToBBHToTauTau_M-2600-NLO",
    "SUSYGluGluToBBHToTauTau_M-2900-NLO",
    "SUSYGluGluToBBHToTauTau_M-3200-NLO",
    "SUSYGluGluToBBHToTauTau_M-3500-NLO",
    "SUSYGluGluToHToTauTau_M-80",
    "SUSYGluGluToHToTauTau_M-90",
    "SUSYGluGluToHToTauTau_M-100",
    "SUSYGluGluToHToTauTau_M-110",
    "SUSYGluGluToHToTauTau_M-120",
    "SUSYGluGluToHToTauTau_M-130",
    "SUSYGluGluToHToTauTau_M-140",
    "SUSYGluGluToHToTauTau_M-180",
    "SUSYGluGluToHToTauTau_M-200",
    "SUSYGluGluToHToTauTau_M-250",
    "SUSYGluGluToHToTauTau_M-300",
    "SUSYGluGluToHToTauTau_M-350",
    "SUSYGluGluToHToTauTau_M-400",
    "SUSYGluGluToHToTauTau_M-450",
    "SUSYGluGluToHToTauTau_M-600",
    "SUSYGluGluToHToTauTau_M-700",
    "SUSYGluGluToHToTauTau_M-800",
    "SUSYGluGluToHToTauTau_M-900",
    "SUSYGluGluToHToTauTau_M-1200",
    "SUSYGluGluToHToTauTau_M-1400",
    "SUSYGluGluToHToTauTau_M-1500",
    "SUSYGluGluToHToTauTau_M-1600",
    "SUSYGluGluToHToTauTau_M-1800",
    "SUSYGluGluToHToTauTau_M-2000",
    "SUSYGluGluToHToTauTau_M-2300",
    "SUSYGluGluToHToTauTau_M-2600",
    "SUSYGluGluToHToTauTau_M-2900",
    "SUSYGluGluToHToTauTau_M-3200",
    'GluGluHToWWTo2L2Nu_M-125',
    'VBFHToWWTo2L2Nu_M-125',
    'GluGluZH_HToWW',
    'HZJ_HToWW',
    'HWminusJ_HToWW',
    'HWplusJ_HToWW',
    'GluGluHToTauTau_M-125',
    'GluGluHToTauTau_M-125-ext',
    'GluGluToHToTauTau_M-125',
    'GluGluToHToTauTau_M-125-ext',
    'VBFHToTauTau_M-125',
    'WplusHToTauTau_M-125',
    'WminusHToTauTau_M-125',
    'ZHToTauTau_M-125',
    'ttHToTauTau',
    'ZHiggs0L1ToTauTau',
    'ZHiggs0L1ZgToTauTau',
    'ZHiggs0L1Zgf05ph0ToTauTau',
    'ZHiggs0L1f05ph0ToTauTau',
    'ZHiggs0MToTauTau',
    'ZHiggs0Mf05ph0ToTauTau',
    'ZHiggs0PHToTauTau',
    'ZHiggs0PHf05ph0ToTauTau',
    'ZHiggs0PMToTauTau',
    'GluGluToHToTauTauPlusTwoJets_M125_minlo',
    'SUSYGluGluToBBHToTauTau_M-1000_powheg',
    'SUSYGluGluToBBHToTauTau_M-100_powheg',
    'SUSYGluGluToBBHToTauTau_M-1200_powheg',
    'SUSYGluGluToBBHToTauTau_M-120_powheg',
    'SUSYGluGluToBBHToTauTau_M-125_powheg',
    'SUSYGluGluToBBHToTauTau_M-130_powheg',
    'SUSYGluGluToBBHToTauTau_M-1400_powheg',
    'SUSYGluGluToBBHToTauTau_M-140_powheg',
    'SUSYGluGluToBBHToTauTau_M-1600_powheg',
    'SUSYGluGluToBBHToTauTau_M-160_powheg',
    'SUSYGluGluToBBHToTauTau_M-1800_powheg',
    'SUSYGluGluToBBHToTauTau_M-180_powheg',
    'SUSYGluGluToBBHToTauTau_M-2000_powheg',
    'SUSYGluGluToBBHToTauTau_M-200_powheg',
    'SUSYGluGluToBBHToTauTau_M-2300_powheg',
    'SUSYGluGluToBBHToTauTau_M-250_powheg',
    'SUSYGluGluToBBHToTauTau_M-2600_powheg',
    'SUSYGluGluToBBHToTauTau_M-2900_powheg',
    'SUSYGluGluToBBHToTauTau_M-300_powheg',
    'SUSYGluGluToBBHToTauTau_M-3200_powheg',
    'SUSYGluGluToBBHToTauTau_M-3500_powheg',
    'SUSYGluGluToBBHToTauTau_M-350_powheg',
    'SUSYGluGluToBBHToTauTau_M-400_powheg',
    'SUSYGluGluToBBHToTauTau_M-450_powheg',
    'SUSYGluGluToBBHToTauTau_M-500_powheg',
    'SUSYGluGluToBBHToTauTau_M-600_powheg',
    'SUSYGluGluToBBHToTauTau_M-60_powheg',
    'SUSYGluGluToBBHToTauTau_M-700_powheg',
    'SUSYGluGluToBBHToTauTau_M-800_powheg',
    'SUSYGluGluToBBHToTauTau_M-80_powheg',
    'SUSYGluGluToBBHToTauTau_M-900_powheg',
    'SUSYGluGluToHToTauTau_M-1000_powheg',
    'SUSYGluGluToHToTauTau_M-100_powheg',
    'SUSYGluGluToHToTauTau_M-1200_powheg',
    'SUSYGluGluToHToTauTau_M-120_powheg',
    'SUSYGluGluToHToTauTau_M-125_powheg',
    'SUSYGluGluToHToTauTau_M-130_powheg',
    'SUSYGluGluToHToTauTau_M-1400_powheg',
    'SUSYGluGluToHToTauTau_M-140_powheg',
    'SUSYGluGluToHToTauTau_M-1600_powheg',
    'SUSYGluGluToHToTauTau_M-160_powheg',
    'SUSYGluGluToHToTauTau_M-1800_powheg',
    'SUSYGluGluToHToTauTau_M-180_powheg',
    'SUSYGluGluToHToTauTau_M-2000_powheg',
    'SUSYGluGluToHToTauTau_M-200_powheg',
    'SUSYGluGluToHToTauTau_M-2300_powheg',
    'SUSYGluGluToHToTauTau_M-250_powheg',
    'SUSYGluGluToHToTauTau_M-2600_powheg',
    'SUSYGluGluToHToTauTau_M-2900_powheg',
    'SUSYGluGluToHToTauTau_M-300_powheg',
    'SUSYGluGluToHToTauTau_M-3200_powheg',
    'SUSYGluGluToHToTauTau_M-3500_powheg',
    'SUSYGluGluToHToTauTau_M-350_powheg',
    'SUSYGluGluToHToTauTau_M-400_powheg',
    'SUSYGluGluToHToTauTau_M-450_powheg',
    'SUSYGluGluToHToTauTau_M-500_powheg',
    'SUSYGluGluToHToTauTau_M-600_powheg',
    'SUSYGluGluToHToTauTau_M-60_powheg',
    'SUSYGluGluToHToTauTau_M-700_powheg',
    'SUSYGluGluToHToTauTau_M-800_powheg',
    'SUSYGluGluToHToTauTau_M-80_powheg',
    'SUSYGluGluToHToTauTau_M-900_powheg',
    'SUSYGluGluToHToTauTau_M-95_powheg',
    'VBFHToTauTau_M-95',
	]

channel = ['em','et','mt','tt']

subdirs=['']
subdirs+=list_paths(outputf)

new_subdirs=[]
for d in subdirs:
  infi=os.listdir('%(outputf)s/%(d)s' % vars())
  if infi: new_subdirs.append((d,infi))
subdirs=new_subdirs

#print subdirs

failed = []

for sa in sample_list:
  sa = 'svfit_'+sa
  for ch in channel:
    remove=True
    to_remove=[]
    hadd_dirs=[]
    command=''
    year = 2017
    if batch:
      JOB='jobs/hadd_%s_%s_%s.sh' % (sa,ch,year)
      os.system('%(JOBWRAPPER)s "" %(JOB)s' %vars())
    for jsdir in subdirs:
      sdir = jsdir[0]
      infiles=jsdir[1]
      if os.path.isfile('%(outputf)s/%(sdir)s/%(sa)s_2017_%(ch)s_0_0_input.root'%vars()):
        if not batch:  
          print "Hadding in subdir %(sdir)s"%vars()
          print "Hadding %(sa)s_%(ch)s in %(sdir)s"%vars()
          os.system('hadd -f %(outputf)s/%(sdir)s/%(sa)s_%(ch)s_2017_input.root %(outputf)s/%(sdir)s/%(sa)s_2017_%(ch)s_*input.root &> ./haddout.txt'% vars()) 
          os.system("sed -i '/Warning in <TInterpreter::ReadRootmapFile>/d' ./haddout.txt")
          filetext = open("./haddout.txt").read()
          if 'Warning' in filetext or 'Error' in filetext:
            print "Hadd had a problem:"
            print filetext
            remove=False 
            failed.append(sa) 
          else :
            to_remove.append('rm %(outputf)s/%(sdir)s/%(sa)s_2017_%(ch)s_*input.root' %vars())
            os.system('rm %(outputf)s/%(sdir)s/%(sa)s_2017_%(ch)s_*input.root' %vars())
        else:
          haddout='haddout_%s_%s_%s.txt' % (sa,ch,sdir) 
          command+="echo \"Hadding %(sa)s_%(ch)s in %(sdir)s\"\necho \"Hadding %(sa)s_%(ch)s\"\nhadd -f %(outputf)s/%(sdir)s/%(sa)s_%(ch)s_2017_input.root %(outputf)s/%(sdir)s/%(sa)s_2017_%(ch)s_*input.root &> ./%(haddout)s\nsed -i '/Warning in <TInterpreter::ReadRootmapFile>/d' ./%(haddout)s\nif [ \"$(cat %(haddout)s | grep -e Warning -e Error)\"  == \"\" ]; then rm %(outputf)s/%(sdir)s/%(sa)s_2017_%(ch)s_*input.root; fi\n" % vars()    
          #command+="echo \"Hadding %(sa)s_%(ch)s in %(sdir)s\"\necho \"Hadding %(sa)s_%(ch)s\"\nhadd -f %(outputf)s/%(sdir)s/%(sa)s_%(ch)s_2017_input.root %(outputf)s/%(sdir)s/%(sa)s_2017_%(ch)s_*input.root &> ./%(haddout)s\nsed -i '/Warning in <TInterpreter::ReadRootmapFile>/d' ./%(haddout)s\nif [ \"$(cat %(haddout)s | grep -e Warning -e Error)\"  == \"\" ]; then echo 'Error in this file'; fi\n" % vars()

    if batch and command:
      with open(JOB, "a") as file: 
        file.write("\n%s" % command)
        file.write('\necho End of job &> jobs/hadd_svfit_%(sa)s_%(ch)s_%(year)s.log' % vars())
      os.system('%(JOBSUBMIT)s %(JOB)s' % vars())

  #if not batch and remove:
  #  # if all channels and systematics were hadded sucsessfully then remove the input files
  #  for x in to_remove:
  #    os.system(x)
