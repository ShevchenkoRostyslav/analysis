#!/usr/bin/python

'''Macro to run Bias test using RooFit cpp macro
    
'''

import os,sys, shutil
from subprocess import Popen,call

#user includes
from bias_config import *

__author__ = "Rostyslav Shevchenko"
__maintainer__ = "Rostyslav Shevchenko"
__email__ = "rostyslav.shevchenko@desy.de"

def PrepareSubmissionCommand(m,signal_strength,n_observed,n_toys,pdf,alternative,bg_only,fix_turnon,parallel,parallel_job):
    """Function to prepare command for submission
    
    """
    command = macro + ' -m ' + str(m) + ' -n ' + str(n_observed) + ' -r ' + str(signal_strength) + ' -t ' + str(n_toys) + ' -f ' + pdf + ' -g ' + alternative
    if(bg_only): command += ' --bg_only'
    if(fix_turnon): command += ' --fix_turnon '
    if(parallel): 
        command += ' -o part' + str(parallel_job)
        if parallel_job == 0: command += ' --plots'
    return command

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

def submit(cshFile,submitter,path_to_stdout = ''):
    """Method to submit jobs.
    
    """
     
    if path_to_stdout == '': path_to_stdout = os.path.splitext(os.path.abspath(cshFile))[0]
    else: 
        if not os.path.exists(path_to_stdout): path_to_stdout = MakeCleanDir(path_to_stdout)
    path_to_stdout += '/'
    
    #Name without path and extension
    name_of_csh = os.path.splitext(os.path.basename(cshFile))[0]
        
    if submitter == "naf":
        command_to_submit = "qsub -cwd -V -l h_rt=1:00:00 -l h_vmem=1G"
        command_for_stdout  = "-o " + path_to_stdout + name_of_csh + '.o'
        command_for_stderr  = "-e " + path_to_stdout + name_of_csh + '.e'
        command_to_run = command_to_submit + ' ' + command_for_stdout + ' ' + command_for_stderr + ' ' + cshFile
    
    if submitter == "shell":
        command_to_submit = "csh -f"
        command_for_stdout = ' >> ' + path_to_stdout + name_of_csh + '.o'
        command_to_run = command_to_submit + ' ' + cshFile #+ ' ' + command_for_stdout
    
    print (command_to_run)
    proc = Popen(command_to_run,shell=True)


def GetSubrange(mass):
    '''Function to get subrange string accroding to the mass value.
    
    '''
    if mass in  [300,350,400]:
        return 'sr1'
    elif mass in [500,600,700,900]:
        return 'sr2'
    elif mass in [1100,1300]:
        return 'sr3'
    else:
        raise AttributeError('Wrong mass has been provided: ' + str(mass))
    
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

if __name__ == '__main__':

    #log file path:
    log_file_path = '/nfs/dust/cms/user/shevchen/output/bias/'
    #command to run BiasTest:
    macro = 'Bias_study_binned'
    
    dir = MakeCleanDir("submit_bias")
    for m in mass:
        for signal_strength in signal_strengths:
            sub_range = GetSubrange(m)
            if sanity_check:
                alternative[sub_range] = pdf[sub_range]
            #Divide into sub-jobs
            for job in range(0,n_jobs):
                ntoys_pjob = n_toys / n_jobs
                sub_jobs = False
                if n_jobs != 1: 
                    sub_jobs = True
                command = PrepareSubmissionCommand(m, signal_strength, n_observed[sub_range], ntoys_pjob, pdf[sub_range], alternative[sub_range], bg_only,fix_turnon, sub_jobs, job)

                cshFile_name  = macro + '_m' + str(m) + '_r_' + str(signal_strength) + '_F_' + pdf[sub_range] + '_G_' + alternative[sub_range]
                if fix_turnon: cshFile_name += '_fixedTurnOn'
                cshFile_name += '_part_' + str(job)
                cshFile = CreateCshFile(cshFile_name,dir)
                UpdateSubmissionCsh(command,cshFile)
                submit(cshFile,'naf',log_file_path)
#             print(command)
#             proc = Popen(command,shell=True)
