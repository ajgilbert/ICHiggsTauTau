#!/usr/bin/env python

import os

      
# Things to loop over
channels = ['et', 'mt','tt']
methods = ['17']
years = ['2016','2017','2018']
variables = ['mt_tot[20,30,40,50,60,70,80,90,100,110,125,140,160,175,200,225,250,280,320,350,400,450,500,560,630,710,800,890,1000]',
             'm_vis[20,30,40,50,60,70,80,90,100,110,125,140,160,175,200,225,250,280,320,350,400,450,500,560,630,710,800,890,1000]',
             'met[0,10,20,30,40,50,60,70,80,90,100,120,140,200,400]',
             'n_jets(4,0,4)',
             'n_prebjets(4,0,4)', 
             'n_deepbjets(4,0,4)',
             'pt_1[20,25,30,35,40,45,50,55,60,65,70,80,90,100,120,140,200,400]',
             'pt_2[20,25,30,35,40,45,50,55,60,65,70,80,90,100,120,140,200,400]',
             'mt_1(20,0,200)',
             'mt_2(20,0,200)',
             'mt_lep(20,0,200)',
             'iso_1[0.0,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.135,0.14,0.145,0.15]']


config_files = {'2016':'scripts/plot_mssm_2016.cfg',
                '2017':'scripts/plot_mssm_2017.cfg',
                '2018':'scripts/plot_mssm_2018.cfg'
               }


add_options = '--embedding --add_wt=\'wt_tau_trg_mssm*wt_tau_id_mssm*wt_prefire\' --ggh_masses=\'\' --bbh_nlo_masses=\'\' --ratio --norm_bins'
 
cmssw_base = os.getcwd()

output_folder = cmssw_base + '/control_plots'


for year in years:
  for channel in channels:
    for method in methods:
      for var in variables:
        if '[' in var: var_string = var.split('[')[0]
        elif '(' in var: var_string = var.split('(')[0]
   
        cfg = config_files[year]
        run_cmd = 'python %(cmssw_base)s/scripts/HiggsTauTauPlot.py --cfg=\'%(cfg)s\' --channel=%(channel)s --method=%(method)s --outputfolder=\'%(output_folder)s/\' --var=\'%(var)s\' %(add_options)s' % vars()
        run_cmd_w_closure = 'python %(cmssw_base)s/scripts/HiggsTauTauPlotW_2.py --cfg=\'%(cfg)s\' --channel=%(channel)s --method=%(method)s --outputfolder=\'%(output_folder)s/\' --var=\'%(var)s\' %(add_options)s --sel=\'mt_1>70 && n_deepbjets==0\'' % vars()
        run_cmd_qcd_closure = 'python %(cmssw_base)s/scripts/HiggsTauTauPlotQCD_2.py --cfg=\'%(cfg)s\' --channel=%(channel)s --method=%(method)s --outputfolder=\'%(output_folder)s/\' --var=\'%(var)s\' %(add_options)s --do_ss --sel=\'mt_1<50 && iso_1>0.05\'' % vars()

        if channel in ["et","mt"]: run_list = [run_cmd,run_cmd_w_closure,run_cmd_qcd_closure]
        else: run_list = [run_cmd]
        
        for cmd in run_list:
          add_name = ''
          if "W_2" in cmd: add_name='_w_dr'
          elif "QCD_2" in cmd: add_name='_qcd_dr' 
          extra_name = 'method_%(method)s_%(var_string)s%(add_name)s' % vars()
          cmd += ' --extra_name=\'%(extra_name)s\'' % vars()
          if var_string in ["mt_tot","m_vis"]: cmd += ' --log_x'
          job_file = 'jobs/mssm_control_plot_%(year)s_%(channel)s_%(extra_name)s.sh' % vars()
          if os.path.exists(job_file): os.system('rm %(job_file)s' % vars())
          os.system('echo "#!/bin/bash" >> %(job_file)s' % vars())
          os.system('echo "cd %(cmssw_base)s" >> %(job_file)s' % vars())
          os.system('echo "source /vols/grid/cms/setup.sh" >> %(job_file)s' % vars()) 
          os.system('echo "export SCRAM_ARCH=slc6_amd64_gcc481" >> %(job_file)s' % vars())
          os.system('echo "eval \'scramv1 runtime -sh\'" >> %(job_file)s' % vars())
          os.system('echo "source %(cmssw_base)s/scripts/setup_libs.sh" >> %(job_file)s' % vars())
          os.system('echo "ulimit -c 0" >> %(job_file)s' % vars())
          os.system('echo "%(cmd)s" >> %(job_file)s' % vars())
          os.system('echo "rm %(output_folder)s/*.png" >> %(job_file)s' % vars())
          os.system('chmod +x %(job_file)s' % vars())

          print "Created job:",job_file
          error_file = job_file.replace('.sh','_error.log')
          output_file = job_file.replace('.sh','_output.log')

          if os.path.exists(error_file): os.system('rm %(error_file)s' % vars())
          if os.path.exists(output_file): os.system('rm %(output_file)s' % vars())

          os.system('qsub -e %(error_file)s -o %(output_file)s -V -q hep.q -l h_rt=0:180:0 -cwd %(job_file)s' % vars())
           
    

