import json
import os
import sys
import subprocess
import glob
import pandas as pd
from helpers import SpikeGLX_utils, log_from_json
from helpers import run_one_probe_sherlock_gtlist
from create_input_json_sherlock_gtlist_v1 import createInputJson

### FSC Modified 11/15/24 to use "gt_list" flag in CatGT to enable concatenating non-consecutive gates.
### Extra modifications required to deal with non-zero starting gates. 
### Change log:
### - import create_input_json_sherlock_gtlist_v1.py to replace create_input_json_sherlock.py
### - modify processed_data_dir to take in desired directory for each session from all_sessions.csv (e.g. to allow separate directories for different concatenation combinations)
### - read in 'gate_list' from session_details for each session; use this to specify first_gate and thus, specify the run_str used to build new directories and filenames (e.g. _g1_imec0)
### - modify how catGT_cmd_string is built by appending gt_list (build before calling create_input_json.py)
### - modify run_specs to be constructed using the min/max of actual gt_list, rather than default gate_string [0,9]
### - add a bunch of print statements to troubleshoot path errors


class EcephysPipeline:
    def __init__(self, session_details):

        # Session details (from csv file)
        self.uuid = session_details['UUID']
        self.file = session_details['File']
        self.rec_file_stem = self.file.split('/')[1][:-3]
        self.log_string = f'{self.uuid},{self.rec_file_stem},'
        self.animal = session_details['Animal']
        self.config_path = session_details['config']
        self.imro_probe0 = session_details['IMRO_probe0']
        self.imro_probe1 = session_details['IMRO_probe1']
        self.whichprobes = str(session_details['whichprobes'])
        self.gate_list = str(session_details['gate_list'])

        # Set which modules to run per probe; catGT and TPrime are called once for each run
        self.run_CatGT = session_details['run_CatGT']
        self.run_KS = session_details['run_KS']
        self.run_TPrime = session_details['run_TPrime']
        self.module_names = ['kilosort_helper','kilosort_postprocessing','noise_templates','mean_waveforms','quality_metrics','depth_estimation','prePhy_filters']
        self.modules = []
        for module in self.module_names:
            if session_details[module]:
                self.modules.append(module)

        # Load configuration from JSON file. 
        # ** .json file does not have comments; refer back to main 2point0_fcho.py
        # Note the .json file now contains parameters for:
        # -- brain region specific params: refPerMS_dict
        # -- threshold values appropriate for version of kilosort: ksTh_dict
        # -- CatGT params: car_mode, loccar_min, loccar_max, catGT_cmd_string
        # -- NIDAQ inputs: ni_present, ni_extract_string
        # -- KS2 or KS25 parameters: ks_remDup, ks_saveRez, ks_copy_fproc, ks_templateRadius_um,
        #                            ks_whiteningRadius_um, ks_minfr_goodchannels, ks_CAR, ks_nblocks
        #                            ks_doFilter, ks_output_tag
        # -- C_Waves: c_Waves_snr_um
        # -- PSTH_events: event_ex_param_str
        # -- TPrime: sync_period, toStream_sync_params

        with open(self.config_path, 'r') as config_file:
            self.config = json.load(config_file)


        print('\n BASE_DIR: ', self.config['base_dir'])

        # Config details (from json file; not all configs, rest in self.config)
        self.base_dir = self.config['base_dir']
        self.prefix = self.config['prefix']
        self.rec_file_list = os.path.join(self.base_dir, 'Preprocessed_Data/Metadata',
                                          f'{self.prefix}_sessions_fcho_sherlocktest.csv')
        self.raw_data_dir = os.path.join(self.base_dir, 'Raw_Data/Neural')
        self.imro_dir = os.path.join(self.base_dir, 'Raw_Data/Neural/IMRO')
        
        self.processed_data_dir_tmp = session_details['Preprocessed_Directory']
        self.task = session_details['Task']
        self.processed_data_dir = os.path.join(self.base_dir, self.processed_data_dir_tmp)
        self.log_file = os.path.join(self.base_dir, 'Preprocessed_Data/Metadata', f'{self.prefix}_ecephys_log.csv')

        # Open log stream
        self.log_stream = open(self.log_file, 'a+')
        self.log_stream.write('Index,File,Sort_Error,Sort_Error_Description\n')

        print(f'\n PROCESSED_DIR: {self.processed_data_dir}')
        print(f'Task: {self.task} \n')

    def process_sessions(self):


        npx_directory = os.path.join(self.raw_data_dir, self.animal)
        # all output will be written here.
        # Output will be in the standard SpikeGLX directory structure:
        # run_folder/probe_folder/*.bin
        catGT_dest = os.path.join(self.processed_data_dir, self.animal, 'Ecephys', self.rec_file_stem)
        # catGT_dest = os.path.join(processed_data_dir, 'g1', animal, 'Ecephys', rec_file_stem)
        json_directory = catGT_dest

        print('\n \n ==== catGT_dest: ', catGT_dest)

        if not os.path.exists(catGT_dest):
            os.mkdir(catGT_dest)
            print('\n ==== made new directory for catGT_dest')

        # -----------
        # Input data
        # -----------
        # Name for log file for this pipeline run. Log file will be saved in the
        # output destination directory catGT_dest
        # If this file exists, new run data is appended to it
        print('\n ==== RUNNING ', self.rec_file_stem)
        logName = f'{self.log_string}_log.csv'

        # run_specs = name, gate, trigger and probes to process
        # Each run_spec is a list of 4 strings:
        #   undecorated run name (no g/t specifier, the run field in CatGT)
        #   gate index range, as a string (e.g. '0,9'); can be greater than actual range
        #       later code will look for all g indices in this range with the same filename stem
        #       EAJ edit: can't currently handle gate indices with >1 digit or multiple triggers
        #   triggers to process/concatenate, as a string e.g. '0,400', '0,0 for a single file
        #           can replace first limit with 'start', last with 'end'; 'start,end'
        #           will concatenate all trials in the probe folder
        #   probes to process, as a string, e.g. '0', '0,3', '0:3'
        #   brain regions, list of strings, one per probe, to set region specific params
        #           these strings must match a key in the param dictionaries above.

 
        # pull self.whichprobes from .csv so it can be specified per run
        # to run imec0&1: '0,1'
        # to run imec0 only: '0'
        # to run imec1 only: '1,1'


        # to concatenate multiple non-consecutive gates in catGT:
        # - pull self.gate_list from .csv so it can be specified per run
        # - build gt_list flag from desired gates
        gate_list_str = self.gate_list.replace(",","")
        g = [int(item) for item in gate_list_str]
        gt_flag_tmp = []
        for item in g:
            gt_flag_tmp.append(f'{{{item},0,0}}')
            print(gt_flag_tmp)
        gt_flag = ''.join(gt_flag_tmp)
        print(f'\n ~~~~~~~~~ concatenating these gates (gt_flag): {gt_flag} ~~~~~~~~~~~~ ')

        gt_minmax = f'{self.gate_list[0]},{self.gate_list[-1]}'
        print(f'\n new gate_minmax: {gt_minmax} \n')


        ## FSC 12/3/24: modify gate
        run_specs = [
            # [rec_file_stem, 'gate', 'triggers', 'probes', 'regions']
            # [self.rec_file_stem, '0,9', '0,0', self.whichprobes, ['cortex', 'cortex']]
            [self.rec_file_stem, gt_minmax, '0,0', self.whichprobes, ['cortex', 'cortex']]            
        ]
        # nshanks = 4
        # ks_trange_list = ['[0 600]', '[603 1203]', '[1206 1806]', '[1809 Inf]']
        # chanMap_list = []
        nshanks = 1
        ks_trange_list = ['[0 Inf]']
        chanMap_list = ['_HPC', '_MEC']


        # delete the existing CatGT.log
        try:
            os.remove('CatGT.log')
        except OSError:
            pass

        # delete existing Tprime.log
        try:
            os.remove('Tprime.log')
        except OSError:
            pass

        # delete existing C_waves.log
        try:
            os.remove('C_Waves.log')
        except OSError:
            pass

        # check for existence of log file, create if not there
        logFullPath = os.path.join(catGT_dest, logName)
        if not os.path.isfile(logFullPath):
            # create the log file, write header
            log_from_json.writeHeader(logFullPath)

        for spec in run_specs:

            # Make list of probes from the probe string
            prb_list = SpikeGLX_utils.ParseProbeStr(spec[3])

            # build path to the first probe folder; look into that folder
            # to determine the range of trials if the user specified t limits as
            # start and end
            # [first_gate, last_gate] = SpikeGLX_utils.ParseGateStr(spec[1])

            ## FSC 12/3/24: instead of using all existing gates, look only in the specified gate list 
            first_gate = g[0]
            last_gate = g[-1]

            run_folder_name = spec[0] + '_g' + repr(first_gate)
            prb0_fld_name = run_folder_name + '_imec' + prb_list[0]
            prb0_fld = os.path.join(npx_directory, run_folder_name, prb0_fld_name)
            first_trig, last_trig = SpikeGLX_utils.ParseTrigStr(spec[2], prb_list[0], first_gate, prb0_fld)
            trigger_str = repr(first_trig) + ',' + repr(last_trig)


            print(f'\n run_folder_name in sglx_parallel_gtlist.py: {run_folder_name}')
            print(f'\n prb0_fld in sglx_parallel_gtlist.py: {prb0_fld} \n \n ')            



            # from MS fork
            # get list of g-indices to concatenate from data directory
            g_range = '[' + spec[1][0] + '-' + spec[1][-1] + ']'
            g_tocat = sorted(glob.glob(os.path.join(npx_directory, (self.rec_file_stem + '_g' + g_range))))
            glist = ''.join((x[-1] + '-') for x in g_tocat)[:-1]  # g inds separated by dashes, minus the last dash


            print('(OLD) Concatenating g indices ' + glist)

            # loop over all probes to build json files of input parameters
            # initalize lists for input and output json files
            catGT_input_json = []
            catGT_output_json = []
            module_input_json = []
            module_output_json = []
            session_id = []
            catgt_output_dir = []
            data_directory = []

            # first loop over probes creates json files containing parameters for
            # both preprocessing (CatGt) and sorting + postprocessing
            for i, prb in enumerate(prb_list):

                print('PRB_LIST: ', prb_list)

                # create CatGT command for this probe
                print('Creating json file for CatGT on probe: ' + prb)

                # Run CatGT
                catGT_input_json.append(
                    os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_CatGT' + '-input.json'))
                catGT_output_json.append(
                    os.path.join(json_directory, spec[0] + '_g' + glist + '_prb' + prb + '_CatGT' + '-output.json'))

                # if this is the first probe proceessed, process the ni stream with it
                if i == 0 and self.config['ni_present']:
                    catGT_stream_string = '-ni -ap -lf'
                    extract_string = self.config['ni_extract_string']
                    print(f'\n ~~~~ ni_extract_string: {extract_string} \n')
                else:
                    catGT_stream_string = '-ap -lf'
                    extract_string = ''

                # build name of first trial to be concatenated/processed;
                # allows reading of the metadata
                run_str = spec[0] + '_g' + str(first_gate)
                run_folder = run_str
                prb_folder = run_str + '_imec' + prb
                input_data_directory = os.path.join(npx_directory, run_folder, prb_folder)
                fileName = run_str + '_t' + repr(first_trig) + '.imec' + prb + '.ap.bin'
                continuous_file = os.path.join(input_data_directory, fileName)
                metaName = run_str + '_t' + repr(first_trig) + '.imec' + prb + '.ap.meta'
                input_meta_fullpath = os.path.join(input_data_directory, metaName)

                print(f'\n input_meta_fullpath: {input_meta_fullpath} \n')

                ### check catGT command
                catGT_cmd_string_gtlist=self.config['catGT_cmd_string'] + f' -gtlist={gt_flag}' + ' ' + extract_string
                print(f'\n +=+=+=+=+=+=+=+= 11/26/24: new catGT command: {catGT_cmd_string_gtlist}\n ')
                tmp2 = self.config['catGTPath']
                print(f'\n catGT_path in sglx_parallel_gtlist.py: {tmp2}')


                info = createInputJson(catGT_input_json[i],
                                       npx_directory=npx_directory,
                                       ecephys_directory=self.config['ecephys_directory'],
                                       kilosort_repository=self.config['kilosort_repository'],
                                       npy_matlab_repository=self.config['npy_matlab_repository'],
                                       catGTPath=self.config['catGTPath'],
                                       tPrime_path=self.config['tPrime_path'],
                                       cWaves_path=self.config['cWaves_path'],
                                       default_ks_tmp=self.config['default_ks_tmp'],
                                       continuous_file=continuous_file,
                                       kilosort_output_directory=catGT_dest,
                                       spikeGLX_data=True,
                                       input_meta_path=input_meta_fullpath,
                                       catGT_run_name=spec[0],
                                       gate_string=spec[1],
                                       trigger_string=trigger_str,
                                       probe_string=prb,
                                       catGT_stream_string=catGT_stream_string,
                                       catGT_car_mode=self.config['car_mode'],
                                       catGT_loccar_min_um=self.config['loccar_min'],
                                       catGT_loccar_max_um=self.config['loccar_max'],
                                       #catGT_cmd_string=self.config['catGT_cmd_string'] + ' ' + extract_string,
                                       catGT_cmd_string = catGT_cmd_string_gtlist,
                                       extracted_data_directory=catGT_dest
                                       )

                # create json files for the other modules
                session_id.append(spec[0] + '_imec' + prb)
                module_input_json.append(os.path.join(json_directory, session_id[i] + '-input.json'))

                # location of the binary created by CatGT, using -out_prb_fld
                run_folder = 'catgt_' + run_str
                prb_folder = run_str + '_imec' + prb
                catgt_output_dir = os.path.join(catGT_dest, run_folder)
                data_directory.append(os.path.join(catGT_dest, run_folder, prb_folder))
                fileName = run_str + '_tcat.imec' + prb + '.ap.bin'
                continuous_file = os.path.join(data_directory[i], fileName)
                print(f'\n catgt_output_dir in sglx_parallel_gtlist.py: {catgt_output_dir}')

                # all post-processing modules alter KS output for Phy, so store the original
                ks_make_copy = True

                for shank in range(nshanks):
                    if nshanks == 1:
                        outputName = 'imec' + prb + '_ks'

                        # FSC: deal with having different IMROs across sessions for the same probe; specify in all_sessions.csv
                        if i == 0:  #
                            imroSelect = self.imro_probe0
                        elif i == 1:
                            imroSelect = self.imro_probe1
                        #chanMap = os.path.join(imro_dir, animal+chanMap_list[i]+'.mat'
                        chanMap = os.path.join(self.imro_dir +'/' + imroSelect + '.mat')
                        #chanMap = os.path.join(self.imro_dir, str('_'), imroSelect, '.mat')
                    else:
                        outputName = 'imec' + prb + 'shank' + str(shank) + '_ks'
                        chanMap = os.path.join(self.imro_dir, 'shank' + str(shank) + 'bank0.mat')
                    ks_trange = ks_trange_list[shank]

                    kilosort_output_dir = os.path.join(data_directory[i], outputName)
                    print(f'\n kilosort_output_dir in sglx_parallel_gtlist.py: {kilosort_output_dir}')

                    # get region specific parameters
                    ks_Th = self.config['ksTh_dict'].get(spec[4][i])
                    refPerMS = self.config['refPerMS_dict'].get(spec[4][i])
                    print('ks_Th: ' + repr(ks_Th) + ' ,refPerMS: ' + repr(refPerMS))

                    info = createInputJson(module_input_json[i],
                                           npx_directory=npx_directory,
                                           ecephys_directory=self.config['ecephys_directory'],
                                           kilosort_repository=self.config['kilosort_repository'],
                                           npy_matlab_repository=self.config['npy_matlab_repository'],
                                           catGTPath=self.config['catGTPath'],
                                           tPrime_path=self.config['tPrime_path'],
                                           cWaves_path=self.config['cWaves_path'],
                                           default_ks_tmp=self.config['default_ks_tmp'],
                                           continuous_file=continuous_file,
                                           spikeGLX_data=True,
                                           input_meta_path=input_meta_fullpath,
                                           kilosort_output_directory=kilosort_output_dir,
                                           ks_make_copy=ks_make_copy,
                                           noise_template_use_rf=False,
                                           catGT_run_name=session_id[i],
                                           gate_string=spec[1],
                                           probe_string=spec[3],
                                           ks_remDup=self.config['ks_remDup'],
                                           ks_finalSplits=1,
                                           ks_labelGood=1,
                                           ks_saveRez=self.config['ks_saveRez'],
                                           ks_copy_fproc=self.config['ks_copy_fproc'],
                                           ks_minfr_goodchannels=self.config['ks_minfr_goodchannels'],
                                           ks_whiteningRadius_um=self.config['ks_whiteningRadius_um'],
                                           ks_doFilter=self.config['ks_doFilter'],
                                           ks_Th=ks_Th,
                                           ks_CSBseed=1,
                                           ks_LTseed=1,
                                           ks_templateRadius_um=self.config['ks_templateRadius_um'],
                                           ks_nblocks=self.config['ks_nblocks'],
                                           ks_CAR=self.config['ks_CAR'],
                                           extracted_data_directory=data_directory[i],
                                           event_ex_param_str=self.config['event_ex_param_str'],
                                           c_Waves_snr_um=self.config['c_Waves_snr_um'],
                                           qm_isi_thresh=refPerMS / 1000,
                                           ks_trange=ks_trange,
                                           ks_chanMap=chanMap,
                                           halfwidth_max=0.35,
                                           )


                    ## FSC 12/3/24: check which gate strings are being passed onto modules
                    tmp3 = info['catGT_helper_params']['gate_string']
                    print(f'\n json created for CatGT/KS: info[gate_string]: {tmp3} \n')

                    # Run each module --- CatGT and KS are run here ---
                    testsherlock_runoneprobe.runOne(session_id[i],
                                         json_directory,
                                         data_directory[i],
                                         self.run_CatGT,
                                         catGT_input_json[i],
                                         catGT_output_json[i],
                                         self.modules,
                                         module_input_json[i],
                                         logFullPath)

                    if self.run_TPrime:
                        # after loop over probes, run TPrime to create files of
                        # event times -- edges detected in auxialliary files and spike times
                        # from each probe -- all aligned to a reference stream.

                        # create json files for calling TPrime
                        input_json = os.path.join(json_directory,
                                                  spec[0] + '_g' + glist + '_prb' + prb + '_TPrime' + '-input.json')
                        output_json = os.path.join(json_directory,
                                                   spec[0] + '_g' + glist + '_prb' + prb + '_TPrime' + '-input.json')

                        info = createInputJson(input_json,
                                               npx_directory=npx_directory,
                                               ecephys_directory=self.config['ecephys_directory'],
                                               kilosort_repository=self.config['kilosort_repository'],
                                               npy_matlab_repository=self.config['npy_matlab_repository'],
                                               catGTPath=self.config['catGTPath'],
                                               tPrime_path=self.config['tPrime_path'],
                                               cWaves_path=self.config['cWaves_path'],
                                               default_ks_tmp=self.config['default_ks_tmp'],
                                               continuous_file=continuous_file,
                                               spikeGLX_data=True,
                                               input_meta_path=input_meta_fullpath,
                                               catGT_run_name=spec[0],
                                               kilosort_output_directory=kilosort_output_dir,
                                               extracted_data_directory=catGT_dest,
                                               tPrime_ni_ex_list=self.config['ni_extract_string'],
                                               event_ex_param_str=self.config['event_ex_param_str'],
                                               sync_period=1.0,
                                               toStream_sync_params=self.config['toStream_sync_params'],
                                               tPrime_3A=False,
                                               toStream_path_3A=' ',
                                               fromStream_list_3A=list(),
                                               gate_string=spec[1],
                                               ks_output_tag=self.config['ks_output_tag']
                                               )

                        command = sys.executable + " -W ignore -m ecephys_spike_sorting.modules." + 'tPrime_helper' + " --input_json " + input_json \
                                  + " --output_json " + output_json
                        subprocess.check_call(command.split(' '))

                    #    log_string+='False,'
        

        # here exception!
        # except Exception as e:
        #    log_string+=f'True,{e}'

        self.log_stream.write(self.log_string + '\n')
        # except:
        #    print("Error running current file")

    def close(self):
        # Close the log stream when done
        self.log_stream.close()


# Usage
if __name__ == "__main__":
    print('\n \n ==== RUNNING sglx_parallel_gtlist_v2.py ====== \n \n ')
    print('sysargv: ', sys.argv[1])
    session_details = json.loads(sys.argv[1])

    pipeline = EcephysPipeline(session_details)
    pipeline.process_sessions()
    pipeline.close()


