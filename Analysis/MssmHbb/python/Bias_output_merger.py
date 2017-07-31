'''Module to merge output from the Bias test

'''

from runBiasTest import GetSubrange,MakeCleanDir
import os
import glob
import subprocess
from shutil import copy, rmtree
#user includes
from bias_config import *

__author__ = "Rostyslav Shevchenko"
__maintainer__ = "Rostyslav Shevchenko"
__email__ = "rostyslav.shevchenko@desy.de"

def RemoveSubDirs(dirs):
    '''Module to remove the directories from the list recursively.
    
    '''
    for d in dirs:
        rmtree(d)
   

def MergeFiles(base_name,input_folder,sub_name='part',first_sub_name = True, plot_name = 'Toy'):
    '''Function to merge files accordint to the base_name 
    and partial names of sub_files.
    
    '''
    os.chdir(input_folder)
    if(first_sub_name): template_name = sub_name + '*' + base_name
    else: template_name = base_name + '*' + sub_name + '*'
    #list of the root files in the folders
    root_list = glob.glob(template_name + '/' + template_name + '.root')
    #exception if no sub-folders
    if len(root_list) == 0:
        print('No sub-folders for ' + base_name) 
        return ['']
    #clean folder for an output
    if not os.path.exists(base_name): os.makedirs(base_name)
    try:
        root_list_joined = ' '.join(root_list)
        haddd_command = 'hadd -f ' + base_name + '/' + base_name + additinoal_out_name + '.root ' + root_list_joined
        subprocess.check_call(haddd_command,shell=True)
    except subprocess.CalledProcessError:
        print ('Error for hadd -f ' + 'hadd -f ' + base_name)
        return ['']
    
    [copy(f,base_name + '/' + f.split('/')[1]) for f in glob.glob(root_list[0].split('/')[0] +  '/' + plot_name + '_*')]

    return glob.glob(template_name)
    
    

if __name__ == '__main__':
    #clean partial dirs
    clean_parts = True
    bias_out = './'
    for m in mass:
        sub_range = GetSubrange(m)
        if sanity_check: alternative[sub_range] = pdf[sub_range]
        for mu in signal_strengths:
            base_name = 'BiasTest_toy' + alternative[sub_range] + '_fit' + pdf[sub_range] + '_' + str(m) + 'GeV_r' + str(mu)
            if(bg_only): base_name += '_BGonly' 
            if(fix_turnon): base_name += '_fixedTurnOn'
            sub_dirs = MergeFiles(base_name,bias_out)
            if(clean_parts and sub_dirs[0] != ''): RemoveSubDirs(sub_dirs)
    
