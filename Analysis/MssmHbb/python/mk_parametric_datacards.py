from subprocess import call
import shutil
import os
import exceptions

cmsswBase = os.environ['CMSSW_BASE']

def SetupDics(mass, model):
    """Function that returns dictionary for current mass.
    
    """
    # Number of observations
    obs_sub_r1 = 259399#358932#
    obs_sub_r2 = 105053#147630#
    obs_sub_r3 = 26760#38906#
    #model dependent nuisances:
    #pdfas that affects signal x-section
    pdfas = {300: '1.014',
             350: '1.015',
             400: '1.016',
             500: '1.018',
             600: '1.022',
             700: '1.027',
             900: '1.034',
             1100:'1.040',
             1300:'1.046'}
    #qcd scale that affects signal xsection
    scale = {300: '1.039',
             350: '1.035',
             400: '1.032',
             500: '1.031',
             600: '1.029',
             700: '1.027',
             900: '1.026',
             1100:'1.024',
             1300:'1.024'}
    
    #Bias error for SuperDiJet-2
#     bias = {300: '2.82359',
#             350: '4.3794377',
#             400: '-1.9990915',
#             500: '0.3496101',
#             600: '0.195367725',
#             700: '-0.06933306',
#             900: '-0.061285',
#             1100:'-0.0741005',
#             1300:'-0.050357485'}
    #Bias error for Bernstein-8
#     bias = {300: '0.9882565',
#             350: '0.61973175',
#             400: '-0.02324525',
#             500: '1.1320708',
#             600: '0.195367725',
#             700: '-0.06933306',
#             900: '-0.061285',
#             1100:'-0.0741005',
#             1300:'-0.050357485'}
    #FLAT bias
    bias = {300: '7.058975',
            350: '4.131545',
            400: '2.324525',
            500: '1.66481',
            600: '0.1860645',
            700: '0.1155551',
            900: '0.0721',
            1100:'0.0741005',
            1300:'0.074055125'}
    
    if 'indep' in model:
        for key in pdfas.keys():
            pdfas[key] = ''
        for key in scale.keys():
            scale[key] = ''
    
    # add_path 
    signal_path = cmsswBase + '/src/Analysis/MssmHbb/output/'
    # Specify ditionaries with TAG vs VALUE
    dictionary = {#200 : {'OBSERVATION' : obs_sub_r1, 'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-300', 'BG_SHAPE_WS' : signal_path + 'ReReco_signal_M-300', 'DATA_SHAPE_WS' : signal_path + 'signal_M-300'},
#                  250 : {'OBSERVATION' : obs_sub_r1, 'SIGNAL_SHAPE_WS' : signal_path + 'signal_M-300', 'BG_SHAPE_WS' : signal_path + 'signal_M-300', 'DATA_SHAPE_WS' : signal_path + 'signal_M-300'},
                  300 : {'OBSERVATION'     : obs_sub_r1,
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-300/workspace/signal_workspace.root',
                         'BG_SHAPE_WS'     : signal_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-300/workspace/signal_bias_workspace.root',
                         'DATA_SHAPE_WS'   : signal_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'OFFLINE_SFL' : '1.0002',
                         'ONLINE_SFB'   : '1.0167',
                         'PILEUP'       : '1.0091',
                         'MODEL_PDFAS'  : pdfas.get(300),
                         'MODEL_SCALE'  : scale.get(300),
                         'SHAPE_BG1'    : 'peak           flatParam',
                         'SHAPE_BG2'    : 'tail           flatParam',
                         'SHAPE_BG3'    : 'width          flatParam',
                         'SHAPE_BG4'    : 'par4           flatParam',
                         'SHAPE_BG5'    : 'slope_novoeff  flatParam',#'',#'slope_novoeff  flatParam',#
                         'SHAPE_BG6'    : 'turnon_novoeff flatParam',
                         'MASS'         : '300',
                         'BIAS_ERR'     :  bias.get(300)},#''},#'turnon_novoeff flatParam'},#
                  350 : {'OBSERVATION' : obs_sub_r1,
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-350/workspace/signal_workspace.root',
                         'BG_SHAPE_WS' : signal_path   + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-350/workspace/signal_bias_workspace.root',
                         'DATA_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root', 
                         'OFFLINE_SFL' : '1.00022',
                         'ONLINE_SFB'   : '1.016',
                         'PILEUP'       : '1.0045',
                         'MODEL_PDFAS'  : pdfas.get(350),
                         'MODEL_SCALE'  : scale.get(350),
                         'SHAPE_BG1'    : 'peak           flatParam',
                         'SHAPE_BG2'    : 'tail           flatParam',
                         'SHAPE_BG3'    : 'width          flatParam',
                         'SHAPE_BG4'    : 'par4           flatParam',
                         'SHAPE_BG5'    : 'slope_novoeff  flatParam',#'',#'slope_novoeff  flatParam',
                         'SHAPE_BG6'    : 'turnon_novoeff flatParam',
                         'MASS'         : '350',
                         'BIAS_ERR'     :  bias.get(350)},#''},#'turnon_novoeff flatParam'},
                  400 : {'OBSERVATION' : obs_sub_r1, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-400/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-400/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root', 
                         'OFFLINE_SFL' : '1.00042',
                         'ONLINE_SFB'   : '1.0154',
                         'PILEUP'       : '1.0046',
                         'MODEL_PDFAS'  : pdfas.get(400),
                         'MODEL_SCALE'  : scale.get(400),
                         'SHAPE_BG1'    : 'peak           flatParam',
                         'SHAPE_BG2'    : 'tail           flatParam',
                         'SHAPE_BG3'    : 'width          flatParam',
                         'SHAPE_BG4'    : 'par4           flatParam',
                         'SHAPE_BG5'    : 'slope_novoeff  flatParam',#'',#'slope_novoeff  flatParam',
                         'SHAPE_BG6'    : 'turnon_novoeff flatParam',
                         'MASS'         : '400',
                         'BIAS_ERR'     :  bias.get(400)},#''},#'turnon_novoeff flatParam'},
                  500 : {'OBSERVATION' : obs_sub_r1, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-500/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-500/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root', 
                         'OFFLINE_SFL' : '1.0006',
                         'ONLINE_SFB'   : '1.0149',
                         'PILEUP'       : '1.0045',
                         'MODEL_PDFAS'  : pdfas.get(500),
                         'MODEL_SCALE'  : scale.get(500),
                         'SHAPE_BG1'    : 'peak           flatParam',
                         'SHAPE_BG2'    : 'tail           flatParam',
                         'SHAPE_BG3'    : 'width          flatParam',
                         'SHAPE_BG4'    : 'par4           flatParam',
                         'SHAPE_BG5'    : 'slope_novoeff  flatParam',#'',#'slope_novoeff  flatParam',
                         'SHAPE_BG6'    : 'turnon_novoeff flatParam',
                         'MASS'         : '500',
                         'BIAS_ERR'     :  bias.get(500)},#''},#'turnon_novoeff flatParam'},
                  600 : {'OBSERVATION' : obs_sub_r2, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-600/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-600/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0009',
                         'ONLINE_SFB'   : '1.0153',
                         'PILEUP'       : '1.0018',
                         'MODEL_PDFAS'  : pdfas.get(600),
                         'MODEL_SCALE'  : scale.get(600),
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail1           flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '600',
                         'BIAS_ERR'     :  bias.get(600)},
                  700 : {'OBSERVATION' : obs_sub_r2, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-700/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-700/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0012',
                         'ONLINE_SFB'   : '1.0164',
                         'PILEUP'       : '1.0033',
                         'MODEL_PDFAS'  : pdfas.get(700),
                         'MODEL_SCALE'  : scale.get(700),
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail1           flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '700',
                         'BIAS_ERR'     :  bias.get(700)},
                  900 : {'OBSERVATION' : obs_sub_r2, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-900/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-900/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0016',
                         'ONLINE_SFB'   : '1.0196',
                         'PILEUP'       : '1.0017',
                         'MODEL_PDFAS'  : pdfas.get(900),
                         'MODEL_SCALE'  : scale.get(900),
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail1           flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '900',
                         'BIAS_ERR'     : bias.get(900)},
                  1100: {'OBSERVATION' : obs_sub_r3, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-1100/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr3/background_workspace_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-1100/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr3/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0024',
                         'ONLINE_SFB'   : '1.0231',
                         'PILEUP'       : '1.0026',
                         'MODEL_PDFAS'  : pdfas.get(1100),
                         'MODEL_SCALE'  : scale.get(1100),
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail1           flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '1100',
                         'BIAS_ERR'     : bias.get(1100)},
                  1300: {'OBSERVATION' : obs_sub_r3, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'ReReco_signal_M-1300/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr3/background_workspace_5000bins.root', 
                         'BIAS_SHAPE_WS'   : signal_path + 'ReReco_signal_M-1300/workspace/signal_bias_workspace.root',
                         'DATA_SHAPE_WS' : signal_path + 'ReReco_bg_fit/sr3/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0025',
                         'ONLINE_SFB'   : '1.0272',
                         'PILEUP'       : '1.0036',
                         'MODEL_PDFAS'  : pdfas.get(1300),
                         'MODEL_SCALE'  : scale.get(1300),
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail1           flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '1300',
                         'BIAS_ERR'     : bias.get(1300)}}
    
    if dictionary[mass] != None:
        return dictionary[mass]
    else:
        raise AttributeError("No rulles for mass = " + mass + " were specified in SetupDics")
    

    
def MakeDir(path):
    """Method to create a directory.
    
    """
    if os.path.exists(path):
        shutil.rmtree(path)
    # create new one
    os.makedirs(path)
    # return the path
    return path

def ConstructDataCardName(mass):
    """Function to construct the datacard name.
    
    """
    name = 'hbb_mbb' + str(mass) + '_mssm-13TeV.txt'
    return name

def CopyTemplateDatacard(path_to_template,path_to_new): 
    if not os.path.exists(path_to_template):
        raise IOError("No txt file at " + path_to_template)
    shutil.copyfile(path_to_template,path_to_new)
    
def ReplaceSpecChar(string):
    """Method to replace special characters in string.
    
    """
    temp_string = string
    for i in "/":#"!@#$%^&*()/[]\{};:,.<>?|`~-=_+":
        if i in temp_string:
            temp_string = temp_string.replace(i , '\\'+i)
    return temp_string

def ReplaceTag(path_to_datacard,tag,value):
    """Method to actually replace current TAG for value.
    
    """
    # If Value is a PATH then all / should be used with \/:
    corrected_value = ReplaceSpecChar(str(value))
    if corrected_value != '':
        command = "sed -i 's/" + str(tag) + "/" + str(corrected_value) + "/g' " + path_to_datacard
    else:
        command = "sed -i '/" + str(tag) + "/d' " + path_to_datacard
    call(command,shell=True)

def AdjustDatacard(path_to_datacard,dictionary):
    """Function to adjust datacard for current mass according to dic.
    
    """
    # check whether it exists
    if not os.path.exists(path_to_datacard):
        raise IOException("No txt file at " + path_to_datacard)
    # replace 
    for tag, value in dictionary.iteritems():
        ReplaceTag(path_to_datacard,tag,value)

if __name__ == '__main__':
    # CMSSW Base
    CMSSW_BASE = os.environ['CMSSW_BASE']
    add_path ='/src/Analysis/MssmHbb/datacards/'
    
    #Switch between mssm and independent
    model = 'independent'
    selection = '201706/14/switch_sub_ranges/independent/no_bias/500_SR2/'
    path_to_dir = ( CMSSW_BASE + add_path + selection )
    # Make dir if it doesn't exist
    MakeDir(path_to_dir)

    #To add bias: 'bias_parametric_datacard_template.txt'
    basename_of_template = 'parametric_datacard_template.txt'
    path_to_template = CMSSW_BASE + add_path + basename_of_template

#    masses = [300,350,400,500,700,900,1100,1300]
    masses = [300,350,400,500,600,700,900,1100,1300]
    for m in masses:
        #construct datacard name
        datacard_name = ConstructDataCardName(m)
        # construct path to datacard
        path_to_datacard = path_to_dir + datacard_name
        #Copy template datacard
        CopyTemplateDatacard(path_to_template,path_to_datacard)
        # Get dictionary for particular mass:
        current_dic = SetupDics(m,model)
        AdjustDatacard(path_to_datacard,current_dic)
        print(path_to_datacard + ': DONE')

