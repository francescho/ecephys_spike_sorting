import os, io, json, sys
import numpy as np

### FSC update 11/26/24 to make check compatibility with main_gtlist.py
### change inputs to createInputJson to reflect local walnut version (create_input_json_nidaq.py)

if sys.platform == 'linux':
    import pwd
from helpers import SpikeGLX_utils


def create_samba_directory(samba_server, samba_share):

    if sys.platform == 'linux':
        proc_owner_uid = str(pwd.getpwnam(os.environ['USER']).pw_uid)
        share_string = 'smb-share:server={},share={}'.format(samba_server, samba_share)
        data_dir = os.path.join('/', 'var', 'run', 'user', proc_owner_uid, 'gvfs', share_string)
    else:
        data_dir = r'\\' + os.path.join(samba_server, samba_share)

    return data_dir

## FSC modified 3/10/24 to take directories from .json file rather than hardcoded

def createInputJson(output_file, 
                    npx_directory=None, 
                    ecephys_directory=None,
                    kilosort_repository=None,
                    npy_matlab_repository=None,
                    catGTPath=None,
                    tPrime_path=None,
                    cWaves_path=None,
                    default_ks_tmp=None,
                    continuous_file = None,
                    spikeGLX_data=True,
                    input_meta_path=None,
                    extracted_data_directory=None,
                    kilosort_output_directory=None,
                    ks_make_copy=False,
                    # probe_type='NP1',
                    probe_type='',
                    catGT_run_name='test',
                    gate_string='0',
                    gate_list_string='0',
                    trigger_string='0,0',
                    probe_string='0',
                    depth_est_fig = 0,
                    catGT_stream_string = '-ap',
                    catGT_car_mode = 'gbldmx',
                    catGT_loccar_min_um = 40,
                    catGT_loccar_max_um = 160,
                    # catGT_cmd_string = '-prb_fld -out_prb_fld',
                    catGT_cmd_string = '',                    
                    catGT_maxZ_um = -1,
                    noise_template_use_rf = True,
                    event_ex_param_str = '',
                    # tPrime_im_ex_list = 'SY=0,384,6,500',
                    tPrime_im_ex_list = '',                    
                    # tPrime_ni_ex_list = 'XA=0,1,3,500',
                    tPrime_ni_ex_list = '',
                    sync_period = 1.0,
                    # toStream_sync_params = 'SY=0,384,6,500',
                    toStream_sync_params = '',
                    # niStream_sync_params = 'XA=0,1,3,500',
                    niStream_sync_params = '',                    
                    tPrime_3A = False,
                    toStream_path_3A = ' ',
                    fromStream_list_3A = list(),
                    ks_helper_noise_threshold = 20,
                    ks_doFilter = 0,
                    ks_remDup = 0,                   
                    ks_finalSplits = 1,
                    ks_labelGood = 1,
                    ks_saveRez = 1,
                    ks_copy_fproc = 1,
                    ks_minfr_goodchannels = 0.1,                  
                    ks_whiteningRadius_um = 163,
                    ks_Th = '[10,4]',
                    ks_CSBseed = 1,
                    ks_LTseed = 1,
                    ks_templateRadius_um = 163,
                    ks_nblocks = 5,
                    ks_CAR = 0,
                    ks_output_tag = 'ks',
                    c_Waves_snr_um = 160,
                    wm_spread_thresh = 0.12,
                    wm_site_range = 16,
                    qm_isi_thresh = 1.5/1000,
                    include_pcs = True,
                    ks_nNeighbors_sites_fix = 0,
                    snr_min = 1,
                    halfwidth_max = 0.3,
                    fr_min = 0.05,
                    isi_viol_max = 0.2,
                    n_viol_max = 1,
                    ks_trange = '[0 Inf]',
                    ks_chanMap = '/home/users/fcho/spikeGLX_tools/chanMap.mat'
                    
                    ):   


    ### CHANGE DIRECTORIES (CONSOLIDATED AT TOP) ###
    print('\n ======= create_input_json_sherlock_gtlist_v1.py, reading in from session_config_gtlist.json======== \n')
    print('kilosort dir!: ',kilosort_repository)

    print(f'\n kilosort_output_directory in create_input_json.py: {kilosort_output_directory}')


    # hard coded paths to code on your computer and system
    # ecephys_directory = r'/home/users/fcho/ecephys_spike_sorting/ecephys_spike_sorting' # github repo on sherlock

    # location of kilosort respository and kilosort version
    # kilosort_repository = r'/oak/stanford/groups/giocomo/fcho2/ecephys_localtools/Kilosort-3.0.1/Kilosort-3.0.1' # use patched 3.0.1 (released March 7, 2024 to fix spike holes) 
    # kilosort_repository = r'/oak/stanford/groups/giocomo/fcho2/ecephys_localtools/KilosortLocal' # this is version3
    

    KS2ver = '3.0'      # must equal '3.0', '2.5' or '2.0', and match the kiilosort_repository
    # KS 3.0 does not yet output pcs.
    if KS2ver == '3.0':
        include_pcs = False  # set to false for KS2ver = '3.0'

    # npy_matlab_repository = r'/home/users/fcho/npy-matlab' # repo on sherlock
    # catGTPath = r'/oak/stanford/groups/giocomo/fcho2/ecephys_localtools/spikeGLXtools/CatGT/CatGT-linux'
    # tPrime_path = r'/oak/stanford/groups/giocomo/fcho2/ecephys_localtools/spikeGLXtools/TPrime/TPrime-linux'
    # cWaves_path = r'/oak/stanford/groups/giocomo/fcho2/ecephys_localtools/spikeGLXtools/C_Waves/C_Waves-linux'    

    # ## cool stuff (for parallellzing to prevent collision of temporary ks files)
    # # for config files and kilosort working space
    # default_ks_tmp = r'/oak/stanford/groups/giocomo/fcho2/ecephys_localtools/Kilosort_scratch' 
    kilosort_output_tmp = os.environ.get("KS_OUT_TMP_OVERRIDE", default_ks_tmp) # if first doesn't exist, it will save to default (do this if running serial)

    # derived directory names
    modules_directory = os.path.join(ecephys_directory,'modules')

    if kilosort_output_directory is None \
         and extracted_data_directory is None \
         and npx_directory is None:
        raise Exception('Must specify at least one output directory')


    ### 
    #default ephys params. For spikeGLX, these get replaced by values read from metadata
    sample_rate = 30000
    num_channels = 385    
    reference_channels = [191]
    uVPerBit = 2.34375
    acq_system = 'PXI'
     
    
    if spikeGLX_data:
        # location of the raw data is the continuous file passed from script
        # metadata file should be located in same directory
        # 
        # kilosort output will be put in the same directory as the input raw data,
        # set in kilosort_output_directory passed from script
        # kilososrt postprocessing (duplicate removal) and identification of noise
        # clusters will act on phy output in the kilosort output directory
        #
        # 
        if input_meta_path is not None:
            probe_type, sample_rate, num_channels, reference_channels, \
                uVPerBit, useGeom = SpikeGLX_utils.EphysParams(input_meta_path)  
            print('SpikeGLX params read from meta')
            print('probe type: {:s}, sample_rate: {:.5f}, num_channels: {:d}, uVPerBit: {:.4f}'.format\
                  (probe_type, sample_rate, num_channels, uVPerBit))
            print('reference channels: ' + repr(reference_channels))
        
        print('kilosort output directory: ' + kilosort_output_directory )
        print('chanmap file: ' +  ks_chanMap)
        
    else:
       print('using default values for probe params')
        

            

    # geometry params by probe type. expand the dictionaries to add types
    # vertical probe pitch vs probe type
    vpitch = {'3A': 20, 'NP1': 20, 'NP21': 15, 'NP24': 15, 'NP1100': 6, 'NP1300':20, 'NP2013': 15}  
    hpitch = {'3A': 32, 'NP1': 32, 'NP21': 32, 'NP24': 32, 'NP1100': 6, 'NP1300':48, 'NP2013': 32} 
    nColumn = {'3A': 2, 'NP1': 2, 'NP21': 2, 'NP24': 2, 'NP1100': 8,'NP1300':2, 'NP2013': 2} 
    
    
    # CatGT needs the inner and outer redii for local common average referencing
    # specified in sites

    catGT_loccar_min_sites = int(round(catGT_loccar_min_um/vpitch.get(probe_type)))
    catGT_loccar_max_sites = int(round(catGT_loccar_max_um/vpitch.get(probe_type)))
    # print('loccar min: ' + repr(catGT_loccar_min_sites))
    
    # whiteningRange is the number of sites used for whitening in KIlosort
    # preprocessing. Calculate the number of sites within the user-specified
    # whitening radius for this probe geometery
    # for a Np 1.0 probe, 163 um => 32 sites
    nrows = np.sqrt((np.square(ks_whiteningRadius_um) - np.square(hpitch.get(probe_type))))/vpitch.get(probe_type)
    ks_whiteningRange = int(round(2*nrows*nColumn.get(probe_type)))
    if ks_whiteningRange > 384:
        ks_whiteningRange = 384
    
    # nNeighbors is the number of sites kilosort includes in a template.
    # Calculate the number of sites within that radisu.
    maxNeighbors = 64 # 64 for standard build of KS
    nrows = np.sqrt((np.square(ks_templateRadius_um) - np.square(hpitch.get(probe_type))))/vpitch.get(probe_type)
    ks_nNeighbors = int(round(2*nrows*nColumn.get(probe_type)))
    
    # workaround for nonstandard patterns
    print('ks_nNeighbors_sites_fix: ', ks_nNeighbors_sites_fix)
    if ks_nNeighbors_sites_fix > 0:
        ks_nNeighbors = ks_nNeighbors_sites_fix
    
    if ks_nNeighbors > maxNeighbors:
        ks_nNeighbors = maxNeighbors          
    print('ks_nNeighbors: ' + repr(ks_nNeighbors))
    
    c_waves_radius_sites = int(round(c_Waves_snr_um/vpitch.get(probe_type)))

    # Create string designating temporary output file for KS2 (gets inserted into KS2 config.m file)
    fproc = os.path.join(kilosort_output_tmp,'temp_wh.dat') # full path for temp whitened data file
    fproc_forward_slash = fproc.replace('\\','/')
    fproc_str = "'" + fproc_forward_slash + "'"
    
     
    dictionary = \
    {

        "directories": {
            "ecephys_directory":ecephys_directory,
            "npx_directory": npx_directory,
            "extracted_data_directory": extracted_data_directory,
            "kilosort_output_directory": kilosort_output_directory,
            "kilosort_output_tmp": kilosort_output_tmp
       },

        "common_files": {
            "settings_json" : npx_directory,
            "probe_json" : os.path.join(extracted_data_directory,'probe_json.json')
        },

        "waveform_metrics" : {
            "waveform_metrics_file" : os.path.join(kilosort_output_directory, 'waveform_metrics.csv')
        },
        
        "cluster_metrics" : {
            "cluster_metrics_file" : os.path.join(kilosort_output_directory, 'metrics.csv')
        },

        "ephys_params": {
            "probe_type" : probe_type,
            "sample_rate" : sample_rate,
            "lfp_sample_rate" : 2500,
            "bit_volts" : uVPerBit,
            "num_channels" : num_channels,
            "reference_channels" : reference_channels,
            "vertical_site_spacing" : 10e-6,
            "ap_band_file" : continuous_file,
            "lfp_band_file" : continuous_file.replace('.ap.bin', '.lf.bin'),
            "reorder_lfp_channels" : False,
            "cluster_group_file_name" : 'cluster_group.tsv'
        }, 

        # "extract_from_npx_params" : {
        #     "npx_directory": npx_directory,
        #     "settings_xml": npx_directory,
        #     "npx_extractor_executable": r"C:\Users\svc_neuropix\Documents\GitHub\npxextractor\Release\NpxExtractor.exe",
        #     "npx_extractor_repo": r"C:\Users\svc_neuropix\Documents\GitHub\npxextractor"
        # },
 
        "depth_estimation_params" : {
            "hi_noise_thresh" : 50.0,
            "lo_noise_thresh" : 3.0,
            "save_figure" : depth_est_fig,
            "figure_location" : os.path.join(extracted_data_directory, 'probe_depth.png'),
            "smoothing_amount" : 5,
            "power_thresh" : 2.5,
            "diff_thresh" : -0.06,
            "freq_range" : [0, 10],
            "max_freq" : 150,
            "saline_range_um" : [3700, 3800],
            "n_passes" : 10,
            "air_gap_um" : 1000,
            "time_interval" : 5,
            "skip_s_per_pass" : 10,
            "start_time" : 10
        }, 

        # "median_subtraction_params" : {
        #     "median_subtraction_executable": "C:\\Users\\svc_neuropix\\Documents\\GitHub\\spikebandmediansubtraction\\Builds\\VisualStudio2013\\Release\\SpikeBandMedianSubtraction.exe",
        #     "median_subtraction_repo": "C:\\Users\\svc_neuropix\\Documents\\GitHub\\spikebandmediansubtraction\\",
        # },

        "kilosort_helper_params" : {

            "matlab_home_directory": kilosort_output_tmp,
            "kilosort_repository" : kilosort_repository,
            "npy_matlab_repository" : npy_matlab_repository,
            "kilosort_version" : 2,
            "spikeGLX_data" : True,
            "ks_make_copy": ks_make_copy,
            "surface_channel_buffer" : 15,
            "noise_threshold" : ks_helper_noise_threshold,

            "kilosort2_params" :
            {
                "KSver" : KS2ver,
                "remDup" : ks_remDup,       #these are expressed as int rather than Bool for matlab compatability
                "finalSplits" : ks_finalSplits,
                "labelGood" : ks_labelGood,
                "saveRez" : ks_saveRez,
                "copy_fproc" : ks_copy_fproc,
                "fproc" : fproc_str,
                "chanMap" : f"'{ks_chanMap}'",
                "fshigh" : 150,
                "minfr_goodchannels" : ks_minfr_goodchannels,
                "Th" : ks_Th,
                "lam" : 10,
                "AUCsplit" : 0.9,
                "minFR" : 0.01,#1/50.,
                "momentum" : '[20 400]',
                "sigmaMask" : 30,
                "ThPre" : 8,
                "gain" : uVPerBit,
                "CSBseed" : ks_CSBseed,
                "LTseed" : ks_LTseed,
                "whiteningRange" : ks_whiteningRange,
                "nNeighbors" : ks_nNeighbors,
                "CAR" : ks_CAR,
                "nblocks" : ks_nblocks,
                "trange" : ks_trange,
            }
        },
        
        "pykilosort_helper_params" : {
            "preprocessing_function" : 'kilosort2',           
            "copy_fproc" : ks_copy_fproc,
            "fproc" : fproc_str,
            "seed" : ks_LTseed,
            "ks2_mode" : False,
            "perform_drift_registration" : True,
            "car" : ks_CAR,
            "Th" : ks_Th,
            "ThPre" : 8,
            "lam" : 10,
            "AUCsplit" : 0.9,
            "minFR" : 1/50.,
            "momentum" : '[20 400]',
            "sig_datashift" : 20,
            "sigmaMask" : 30,
            "fshigh" : 300,
            "fslow" : 10000,
            "minfr_goodchannels" : 0,
            "whiteningRange" : ks_whiteningRange,            
            "deterministic_mode" : True,            
            "nblocks" : ks_nblocks,
            "doFilter" : ks_doFilter

        },


        "ks_postprocessing_params" : {
            "align_avg_waveform" : False,              
            "remove_duplicates" : True,
            "cWaves_path" : cWaves_path,
            "within_unit_overlap_window" : 0.00017,
            "between_unit_overlap_window" : 0.00041,
            "between_unit_dist_um" : 66,
            "deletion_mode" : 'lowAmpCluster',
            "include_pcs" : include_pcs
        },

        "mean_waveform_params" : {
        
            "mean_waveforms_file" : os.path.join(kilosort_output_directory, 'mean_waveforms.npy'),
            "samples_per_spike" : 82,
            "pre_samples" : 20,
            "num_epochs" : 1,           #epochs not implemented for c_waves
            "spikes_per_epoch" : 1000,
            "spread_threshold" : wm_spread_thresh,
            "site_range" : wm_site_range,  
            "cWaves_path" : cWaves_path,
            "use_C_Waves" : True,
            "snr_radius" : c_waves_radius_sites,
            "snr_radius_um" : c_Waves_snr_um     
        },
            

        "noise_waveform_params" : {
            "classifier_path" : os.path.join(modules_directory, 'noise_templates', 'rf_classifier.pkl'),
            "multiprocessing_worker_count" : 10,
            "use_random_forest" : noise_template_use_rf
        },

        "quality_metrics_params" : {
            "isi_threshold" : qm_isi_thresh,
            "min_isi" : 0.000166,
            "tbin_sec" : 0.001,
            "max_radius_um" : 68,
            "max_spikes_for_unit" : 500,
            "max_spikes_for_nn" : 10000,
            "n_neighbors" : 4,
            'n_silhouette' : 10000,
            "drift_metrics_interval_s" : 51,
            "drift_metrics_min_spikes_per_interval" : 10,
            "include_pcs" : include_pcs
        },
        
        "catGT_helper_params" : {
            "run_name" : catGT_run_name,
            "gate_string" : gate_string,
            "probe_string" : probe_string,
            "trigger_string": trigger_string,
            "stream_string" : catGT_stream_string,
            "car_mode" : catGT_car_mode,
            "loccar_inner" : catGT_loccar_min_sites,
            "loccar_outer": catGT_loccar_max_sites,
            "loccar_inner_um" : catGT_loccar_min_um,
            "loccar_outer_um" : catGT_loccar_max_um,
            "maxZ_um" : catGT_maxZ_um,
            'useGeom' : useGeom,
            "cmdStr" : catGT_cmd_string,
            "catGTPath" : catGTPath
        },

        "tPrime_helper_params" : {
                "tPrime_path" : tPrime_path,
                "im_ex_list" : tPrime_im_ex_list,
                "ni_ex_list" : tPrime_ni_ex_list,
                "sync_period" : sync_period,
                "toStream_sync_params" : toStream_sync_params,
                "ni_sync_params" : niStream_sync_params,
                "tPrime_3A" : tPrime_3A,
                "toStream_path_3A" : toStream_path_3A,
                "fromStream_list_3A" : fromStream_list_3A,
                "psth_ex_str": event_ex_param_str,
                "sort_out_tag": ks_output_tag
        },  
        
        "prephy_filters_params" : {
                "snr_min" : snr_min,
                "halfwidth_max" : halfwidth_max,
                "fr_min" : fr_min,
                "isi_viol_max" : isi_viol_max,
                "n_viol_max" : n_viol_max,
        }, 
                
        "psth_events": {
                "event_ex_param_str": event_ex_param_str
         }
        
    }

    with io.open(output_file, 'w', encoding='utf-8') as f:
        f.write(json.dumps(dictionary, ensure_ascii=False, sort_keys=True, indent=4))

    return dictionary
