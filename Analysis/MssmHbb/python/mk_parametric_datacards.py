from subprocess import call
import shutil
import os
import exceptions

cmsswBase = os.environ['CMSSW_BASE']
blinded = False
adjoint = False

def SetupDics(mass, model):
    """Function that returns dictionary for current mass.
    
    """
    # Number of observations
    obs_sub_r1 = 259399#358932#
    obs_sub_r2 = 105053#147630#
    obs_sub_r3 = 26760#38906#

    # Model dependent uncertainties for different models
    model_pdf_alphas = {
    'mhmodp_200' : {300: '1.0155', 350: '1.0159', 400: '1.01712', 500: '1.0176', 600: '1.023', 700: '1.02637', 900: '1.03369', 1100:'1.03984', 1300:'1.04598'},
    '2hdm' : {300: '1.0155', 350: '1.0159', 400: '1.01712', 500: '1.0176', 600: '1.023', 700: '1.02637', 900: '1.03369', 1100:'1.03984', 1300:'1.04598'},
    'tau_phobic' : {300: '1.01478', 350: '1.01636', 400: '1.03722', 500: '1.03418', 600: '1.03224', 700: '1.0284', 900: '1.02639', 1100:'1.02405', 1300:'1.02403'},
    'light_stop' : {300: '1.02464', 350: '1.02686', 400: '1.03083', 500: '1.04167', 600: '1.02439', 700: '', 900: '', 1100:'', 1300:''},
    'light_stau' : {300: '1.0191', 350: '1.02131', 400: '1.02443', 500: '1.02801', 600: '1.0365', 700: '1.0379', 900: '1.03369', 1100:'1.03984', 1300:'1.04598'},
    'hMSSM' : {300: '1.01549', 350: '1.01483', 400: '1.01907', 500: '1.02115', 600: '1.02451', 700: '1.02762', 900: '1.03436', 1100:'1.03986', 1300:'1.046'},
    'independent' : {300: '', 350: '', 400: '', 500: '', 600: '', 700: '', 900: '', 1100:'', 1300:''},
    }
    model_scale = {
    'mhmodp_200' : {300: '1.04264', 350: '1.03887', 400: '1.03253', 500: '1.02933', 600: '1.02906', 700: '1.02637', 900: '1.02639', 1100:'1.02405', 1300:'1.02403'},
    '2hdm' : {300: '1.04264', 350: '1.03887', 400: '1.03253', 500: '1.02933', 600: '1.02906', 700: '1.02637', 900: '1.02639', 1100:'1.02405', 1300:'1.02403'},
    'tau_phobic' : {300: '1.04356', 350: '1.04127', 400: '1.03722', 500: '1.03418', 600: '1.03224', 700: '1.0284', 900: '1.02639', 1100:'1.02405', 1300:'1.02403'},
    'light_stop' : {300: '1.06957', 350: '1.06989', 400: '1.063', 500: '1.06644', 600: '1.03223', 700: '', 900: '', 1100:'', 1300:''},
    'light_stau' : {300: '1.05382', 350: '1.05246', 400: '1.04723', 500: '1.04622', 600: '1.04757', 700: '1.03883', 900: '1.02639', 1100:'1.02405', 1300:'1.02403'},
    'hMSSM' : {300: '1.04425', 350: '1.04025', 400: '1.03814', 500: '1.03462', 600: '1.03105', 700: '1.02901', 900: '1.02622', 1100:'1.02403', 1300:'1.02402'},
    'independent' : {300: '', 350: '', 400: '', 500: '', 600: '', 700: '', 900: '', 1100:'', 1300:''},
    }
    
    # add_path 
    signal_path = cmsswBase + '/src/Analysis/MssmHbb/output/'
    bg_path = signal_path
    if not blinded: bg_path += 'unblinded/' 
    # Specify ditionaries with TAG vs VALUE
    dictionary = {
                  300 : {'OBSERVATION'     : obs_sub_r1,
                         'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-300/workspace/signal_workspace.root',
                         'BG_SHAPE_WS'     : bg_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-300/workspace/signal_bias_workspace.root',
                         'DATA_SHAPE_WS'   : bg_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'OFFLINE_SFL' : '1.0002',
                         'ONLINE_SFB'   : '1.0167',
                         'PILEUP'       : '1.0091',
                         'MODEL_PDFAS'  : model_pdf_alphas[model][300],
                         'MODEL_SCALE'  : model_scale[model][300],
                         'SHAPE_BG1'    : 'peak           flatParam',
                         'SHAPE_BG2'    : 'tail           flatParam',
                         'SHAPE_BG3'    : 'width          flatParam',
                         'SHAPE_BG4'    : 'par4           flatParam',
                         'SHAPE_BG5'    : 'slope_novoeff  flatParam',#'',#'slope_novoeff  flatParam',#
                         'SHAPE_BG6'    : 'turnon_novoeff flatParam',
                         'MASS'         : '300',
                         'BIAS_ERR'     :  '7.05244'},#''},#'turnon_novoeff flatParam'},#
                  350 : {'OBSERVATION' : obs_sub_r1,
                         'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-350/workspace/signal_workspace.root',
                         'BG_SHAPE_WS' : bg_path   + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-350/workspace/signal_bias_workspace.root',
                         'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root', 
                         'OFFLINE_SFL' : '1.00022',
                         'ONLINE_SFB'   : '1.016',
                         'PILEUP'       : '1.0045',
                         'MODEL_PDFAS'  : model_pdf_alphas[model][350],
                         'MODEL_SCALE'  : model_scale[model][350],
                         'SHAPE_BG1'    : 'peak           flatParam',
                         'SHAPE_BG2'    : 'tail           flatParam',
                         'SHAPE_BG3'    : 'width          flatParam',
                         'SHAPE_BG4'    : 'par4           flatParam',
                         'SHAPE_BG5'    : 'slope_novoeff  flatParam',#'',#'slope_novoeff  flatParam',
                         'SHAPE_BG6'    : 'turnon_novoeff flatParam',
                         'MASS'         : '350',
                         'BIAS_ERR'     :  '4.13302'},#''},#'turnon_novoeff flatParam'},
                  400 : {'OBSERVATION' : obs_sub_r1, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-400/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-400/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root', 
                         'OFFLINE_SFL' : '1.00042',
                         'ONLINE_SFB'   : '1.0154',
                         'PILEUP'       : '1.0046',
                         'MODEL_PDFAS'  : model_pdf_alphas[model][400],
                         'MODEL_SCALE'  : model_scale[model][400],
                         'SHAPE_BG1'    : 'peak           flatParam',
                         'SHAPE_BG2'    : 'tail           flatParam',
                         'SHAPE_BG3'    : 'width          flatParam',
                         'SHAPE_BG4'    : 'par4           flatParam',
                         'SHAPE_BG5'    : 'slope_novoeff  flatParam',#'',#'slope_novoeff  flatParam',
                         'SHAPE_BG6'    : 'turnon_novoeff flatParam',
                         'MASS'         : '400',
                         'BIAS_ERR'     :  '2.30296'},#''},#'turnon_novoeff flatParam'},
                500 : {'OBSERVATION' : obs_sub_r2, 
                       'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-500/workspace/signal_workspace.root', 
                      'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-500/workspace/signal_bias_workspace.root', 
                       'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root',
                       'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root', 
                       'OFFLINE_SFL' : '1.0006',
                       'ONLINE_SFB'   : '1.0149',
                       'PILEUP'       : '1.0045',
                         'MODEL_PDFAS'  : model_pdf_alphas[model][500],
                         'MODEL_SCALE'  : model_scale[model][500],
                       'SHAPE_BG1'    : 'peak1           flatParam',
                       'SHAPE_BG2'    : 'tail1           flatParam',
                       'SHAPE_BG3'    : 'width1          flatParam',
                        'SHAPE_BG4'    : '',
                        'SHAPE_BG5'    : '',#'',#'slope_novoeff  flatParam',
                        'SHAPE_BG6'    : '',
                       'MASS'         : '500',
                       'BIAS_ERR'     :  '0.287686'},
                  600 : {'OBSERVATION' : obs_sub_r2, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-600/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-600/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0009',
                         'ONLINE_SFB'   : '1.0153',
                         'PILEUP'       : '1.0018',
                         'MODEL_PDFAS'  : model_pdf_alphas[model][600],
                         'MODEL_SCALE'  : model_scale[model][600],
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail1           flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '600',
                         'BIAS_ERR'     :  '0.1606918'},
                  700 : {'OBSERVATION' : obs_sub_r2, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-700/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-700/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0012',
                         'ONLINE_SFB'   : '1.0164',
                         'PILEUP'       : '1.0033',
                         'MODEL_PDFAS'  : model_pdf_alphas[model][700],
                         'MODEL_SCALE'  : model_scale[model][700],
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail1           flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '700',
                         'BIAS_ERR'     :  '0.1135502'},
                  900 : {'OBSERVATION' : obs_sub_r2, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-900/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root',
                         'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-900/workspace/signal_bias_workspace.root', 
                         'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0016',
                         'ONLINE_SFB'   : '1.0196',
                         'PILEUP'       : '1.0017',
                         'MODEL_PDFAS'  : model_pdf_alphas[model][900],
                         'MODEL_SCALE'  : model_scale[model][900],
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail1           flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '900',
                         'BIAS_ERR'     : '0.0715888'},
                1100: {'OBSERVATION' : obs_sub_r3, 
                       'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-1100/workspace/signal_workspace.root', 
                       'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr3/background_workspace_5000bins.root',
                       'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-1100/workspace/signal_bias_workspace.root', 
                       'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr3/background_workspace_5000bins.root', 
                       'OFFLINE_SFL' : '1.0024',
                       'ONLINE_SFB'   : '1.0231',
                       'PILEUP'       : '1.0026',
                        'MODEL_PDFAS'  : model_pdf_alphas[model][1100],
                        'MODEL_SCALE'  : model_scale[model][1100],
                       'SHAPE_BG1'    : 'peak1           flatParam',
                       'SHAPE_BG2'    : 'tail            flatParam',
                       'SHAPE_BG3'    : 'width1          flatParam',
                       'SHAPE_BG4'    : '',
                       'SHAPE_BG5'    : '',
                       'SHAPE_BG6'    : '',
                       'MASS'         : '1100',
                       'BIAS_ERR'     : '0.073350875'},
                  1300: {'OBSERVATION' : obs_sub_r3, 
                         'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-1300/workspace/signal_workspace.root', 
                         'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr3/background_workspace_5000bins.root', 
                         'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-1300/workspace/signal_bias_workspace.root',
                         'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr3/background_workspace_5000bins.root', 
                         'OFFLINE_SFL' : '1.0025',
                         'ONLINE_SFB'   : '1.0272',
                         'PILEUP'       : '1.0036',
                        'MODEL_PDFAS'  : model_pdf_alphas[model][1300],
                        'MODEL_SCALE'  : model_scale[model][1300],
                         'SHAPE_BG1'    : 'peak1           flatParam',
                         'SHAPE_BG2'    : 'tail            flatParam',
                         'SHAPE_BG3'    : 'width1          flatParam',
                         'SHAPE_BG4'    : '',
                         'SHAPE_BG5'    : '',
                         'SHAPE_BG6'    : '',
                         'MASS'         : '1300',
                         'BIAS_ERR'     : '0.07426975'}}
    
    if(adjoint == True): 
        dictionary[500] = {'OBSERVATION' : obs_sub_r1,
                                            'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-500_SR1/workspace/signal_workspace.root',
                                            'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                                            'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-500_SR1/workspace/signal_bias_workspace.root',
                                            'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root',
                                            'OFFLINE_SFL' : '1.0006',
                                            'ONLINE_SFB'   : '1.0149',
                                            'PILEUP'       : '1.0045',
                                            'MODEL_PDFAS'  : model_pdf_alphas[model][500],
                                            'MODEL_SCALE'  : model_scale[model][500],
                                            'SHAPE_BG1'    : 'peak           flatParam',
                                            'SHAPE_BG2'    : 'tail           flatParam',
                                            'SHAPE_BG3'    : 'width          flatParam',
                                            'SHAPE_BG4'    : 'par4           flatParam',
                                            'SHAPE_BG5'    : 'slope_novoeff  flatParam',#'',#'slope_novoeff  flatParam',
                                            'SHAPE_BG6'    : 'turnon_novoeff flatParam',
                                            'MASS'         : '500',
                                            'BIAS_ERR'     :  '1.66481'}
        dictionary[1100] = {'OBSERVATION' : obs_sub_r2, 
                        'SIGNAL_SHAPE_WS' : signal_path + 'test_ReReco_signal_M-1100_SR2/workspace/signal_workspace.root', 
                        'BG_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root',
                        'BIAS_SHAPE_WS'   : signal_path + 'test_ReReco_signal_M-1100_SR2/workspace/signal_bias_workspace.root', 
                        'DATA_SHAPE_WS' : bg_path + 'ReReco_bg_fit/sr2/background_workspace_5000bins.root', 
                        'OFFLINE_SFL' : '1.0024',
                        'ONLINE_SFB'   : '1.0231',
                        'PILEUP'       : '1.0026',
                        'MODEL_PDFAS'  : model_pdf_alphas[model][1100],
                        'MODEL_SCALE'  : model_scale[model][1100],
                        'SHAPE_BG1'    : 'peak1           flatParam',
                        'SHAPE_BG2'    : 'tail1           flatParam',
                        'SHAPE_BG3'    : 'width1          flatParam',
                        'SHAPE_BG4'    : '',
                        'SHAPE_BG5'    : '',
                        'SHAPE_BG6'    : '',
                        'MASS'         : '1100',
                        'BIAS_ERR'     : '0.0741005'}
    
    dictionary[mass]['SIGNAL_GROUP'] = 'signal group = lumi_13TeV CMS_eff_l_13TeV CMS_PU_13TeV CMS_eff_bonl_13TeV pdf_gg QCDscale CMS_scale_j_13TeV CMS_res_j_13TeV CMS_eff_pTonl_13TeV CMS_eff_b_13TeV pdf_Higgs_gg bias' 
    if(model == 'independent'): dictionary[mass]['SIGNAL_GROUP'] = 'signal group = lumi_13TeV CMS_eff_l_13TeV CMS_PU_13TeV CMS_eff_bonl_13TeV CMS_scale_j_13TeV CMS_res_j_13TeV CMS_eff_pTonl_13TeV CMS_eff_b_13TeV pdf_Higgs_gg bias'
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
    
    #Switch between mssm and independent:
    # mhmodp_200 2hdm tau_phobic light_stop light_stau hMSSM independent
    model = 'independent'
    selection = '201709/20/unblinded/' + model + '/'
    if(adjoint): selection += 'adjoint_'
    selection += 'bias/'
    path_to_dir = ( CMSSW_BASE + add_path + selection )
    # Make dir if it doesn't exist
    MakeDir(path_to_dir)

    #To add bias: 'bias_parametric_datacard_template.txt'
    basename_of_template = 'bias_parametric_datacard_template.txt'
    path_to_template = CMSSW_BASE + add_path + basename_of_template

#     masses = [1100]
    masses = [300,350,400,500,600,700,900,1100,1300]
#     masses = [500,1100]
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
