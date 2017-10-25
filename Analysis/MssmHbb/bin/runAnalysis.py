#!/usr/bin/env python
"""Script to run MSSM Hbb analysis with config file.

   This screept can bu used to run MSSM Hbb analysis
   using naf/lxplus/shell.
"""

import argparse
import os
import glob
from subprocess import call, Popen
import shutil # for rmtree function

__author__ = "Rostyslav Shevchenko"
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

def ParseOption():
    """Method to parse user input.
    
    """
    parser = argparse.ArgumentParser(description='MSSM Hbb analysis runner script')
    parser.add_argument('-c',       dest='config_file', type=str, help='config file to be used',default="config.cfg")
    parser.add_argument('-i',       dest='input_folder',type=str, help='Input folder')
    parser.add_argument('-l',       dest='log_folder',  type=str, help='Folder for batch logs', default='submit')
#    Add subparsers for naf/lxplus/shell run options
    subparsers = parser.add_subparsers(dest='submitter')
    #create parser to run jobs inciated in setinput at lxplus machines
    parser_lxplus = subparsers.add_parser('lxplus',help='Run jobs at lxplus')
    parser_lxplus.add_argument('--q', dest='queue', required=True, type=int, help='submission queue(1nh,8nh,1nd...)')
    
    #create parser to run jobs iniciated in setinput at naf machines
    parser_naf = subparsers.add_parser('naf',help='Run jobs at naf')
    
    #create parser to run jobs interactively in the shell
    parser_shell = subparsers.add_parser('shell',help='Run jobs in shell')
    
    args = parser.parse_args()
    # convert output_dir to absolute path
    args.config_file  = os.path.abspath(os.path.join(os.getcwd(), args.config_file))
    args.input_folder = os.path.abspath(os.path.join(os.getcwd(), args.input_folder))
    args.log_folder   = os.path.abspath(os.path.join(os.getcwd(), args.log_folder))
    print 'Path is: ' + args.config_file  
    
    return args

def CreateCshFile(fileName,directory):
    #prepare new csh file to run at the farm
    newCsh = directory + '/' + fileName + '.csh'
    outCsh = open(newCsh, 'w')
    outCsh.write('#!/bin/csh' + '\n')
    outCsh.close()
    call('chmod +x ' + newCsh, shell=True)
    return newCsh

def PrepareConfig(csh_name,config_file,submission_dir,file_to_submit):
    """Method to prepare config file for AllInOneTool usage.
    
    """
    config_name = os.path.splitext(config_file)[0]
    temp_config = os.path.basename(csh_name) + ".cfg"
    call("cp " + config_file + " " + submission_dir + "/" + temp_config,shell=True)
    call("sed -i 's|INPUT|" + str(file_to_submit) + "|g' " + submission_dir + "/" + temp_config,shell=True)
    return submission_dir + "/" + temp_config

def UpdateSubmissionCsh(command, cshFile):
    """Method to update csh script.
    
    """
    if not os.path.exists(cshFile):
        raise AttributeError("ERROR: File " + cshFile + " doens't exist. Please check spelling")
    call('echo $"' + command + '" >> ' + cshFile, shell=True)
    return cshFile

def submit(cshFile,submitter):
    """Method to submit jobs.
    
    """
    if submitter == "naf":
        command_to_submit = "qsub -cwd -V -l h_rt=24:00:00 -l h_vmem=20G"
        command_for_stdout  = "-o " + os.path.splitext(cshFile)[0] + '.o'
        command_for_stderr  = "-e " + os.path.splitext(cshFile)[0] + '.e'
        command_to_run = command_to_submit + ' ' + command_for_stdout + ' ' + command_for_stderr + ' ' + cshFile
    
    if submitter == "shell":
        command_to_submit = "csh -f"
        command_for_stdout = ' >> ' + os.path.splitext(cshFile)[0] + '.o'
        command_to_run = command_to_submit + ' ' + cshFile + ' ' + command_for_stdout

    print command_to_run   
#     call(command_to_submit + ' ' + command_for_stdout + ' ' + command_for_stderr + ' ' + cshFile,shell=True)
#     print command_to_submit + ' ' + cshFile + ' ' + command_for_stdout + ' ' + command_for_stderr
#    proc = Popen(command_to_run,shell=True)
#     call(command_to_submit +  ' ' + cshFile,shell=True)
    
if __name__ == '__main__':
    
    args = ParseOption()
    #Check whether submit dir exists:
    dir = MakeCleanDir(args.log_folder)
    for file in glob.glob(args.input_folder + '/*.txt'):
        print file
        # Prepare submit files
        cshFile = CreateCshFile(os.path.basename(os.path.splitext(file)[0]),dir)
        # Create config files from the template
        config_file = PrepareConfig(os.path.splitext(file)[0],args.config_file,dir,file)
        
        command = "AllInOneTool -c " + config_file
        UpdateSubmissionCsh(command,cshFile)
        submit(cshFile,args.submitter)

    
