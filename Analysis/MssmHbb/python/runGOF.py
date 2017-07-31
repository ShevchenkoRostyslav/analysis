#!/usr/bin/python

'''Macro to run combine tool with Asymptotic method.
    
'''

from runAsymptotic import checkInput, runCombineTool
import os,sys, shutil
from subprocess import call,Popen

__author__ = "Rostyslav Shevchenko"
__maintainer__ = "Rostyslav Shevchenko"
__email__ = "rostyslav.shevchenko@desy.de"

def MakeCleanDir(path_to_dir):
    """Method to create clean dir for submission.
    
    """
    if os.path.exists(path_to_dir):
        shutil.rmtree(path_to_dir)
        # create new one
    os.makedirs(path_to_dir)
    # return the path
    path_to_dir = os.path.abspath(os.path.join(os.getcwd(), path_to_dir))
    return path_to_dir

def CreateCshFile(fileName,directory):
    #prepare new csh file to run at the farm
    newCsh = directory + '/' + fileName + '.csh'
    outCsh = open(newCsh, 'w')
    outCsh.write('#!/bin/csh' + '\n')
    outCsh.close()
    call('chmod +x ' + newCsh, shell=True)
    return newCsh

def UpdateSubmissionCsh(command, cshFile):
    """Method to update csh script.
    
    """
    if not os.path.exists(cshFile):
        raise AttributeError("ERROR: File " + cshFile + " doens't exist. Please check spelling")
    call('echo $"' + command + '" >> ' + cshFile, shell=True)
    return cshFile

def submit(cshFile):
    """Method to submit jobs.
    
    """
    command_to_submit = "qsub -cwd -V -l h_rt=5:00:00 -l h_vmem=2G"
    command_for_stdout  = "-o " + os.path.splitext(cshFile[:-4])[0] + '.o'
    command_for_stderr  = "-e " + os.path.splitext(cshFile[:-4])[0] + '.e'
    command_to_run = command_to_submit + ' ' + command_for_stdout + ' ' + command_for_stderr + ' ' + cshFile
    
#     call(command_to_submit + ' ' + command_for_stdout + ' ' + command_for_stderr + ' ' + cshFile,shell=True)
#     print command_to_submit + ' ' + cshFile + ' ' + command_for_stdout + ' ' + command_for_stderr
    proc = Popen(command_to_run,shell=True)
#     call(command_to_submit +  ' ' + cshFile,shell=True)
    

if __name__ == '__main__':

    #working directory with datacards and stored output:
    datacard_folder = '/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/datacards/201707/20/unblinded/tests/'
    checkInput(datacard_folder)
    os.chdir(datacard_folder)
    #algorithm from https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideHiggsAnalysisCombinedLimit#Goodness_of_fit_tests
    algorithm = ['saturated','KS','AD']
    # number of batch jobs
    Njobs = 300
    # Total number oftoys 
    Ntoys = 3000 
    #list of mass points for those GOFwill be calculated
#     mass = ['300','350','400','500','600','700','900','1100','1300']
    mass = ['350','700','1300']
    #combine preferences to be added
    combine_add = ''
    for m in mass:
        #Data card name
        data_card = 'hbb_mbb' + m + '_mssm-13TeV.txt'
        #Create submission dir 
        dir = MakeCleanDir('naf_' + m)
        for j in range(Njobs+1):
            for algo in algorithm:
                #submission file name
                file = data_card[:-4] + '_job' +  str(j) + '_algo_' + algo +'.txt'
                # Prepare submit files
                cshFile = CreateCshFile(file,dir)
                combine_command = 'combine -M GoodnessOfFit -n ' + str(j) + '_' + algo + ' --algo=' + algo + ' -t ' + str(Ntoys/Njobs) + ' -m ' + m + ' -s 0 ' + data_card
                UpdateSubmissionCsh(combine_command,cshFile)
                submit(cshFile)
                
