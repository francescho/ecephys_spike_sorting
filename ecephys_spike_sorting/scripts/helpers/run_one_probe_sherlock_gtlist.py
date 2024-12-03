# test_runOne.py

import os
import shutil
import subprocess
import sys
from helpers import SpikeGLX_utils, log_from_json

# Include the modified runOne function here
def runOne(session_id,
           json_directory,
           data_directory,
           run_CatGT,
           catGT_input_json,
           catGT_output_json,
           modules,
           module_input_json,
           logFullPath):

        if run_CatGT:
            command = [
                sys.executable,
                "-W", "ignore",
                "-m", "ecephys_spike_sorting.modules.catGT_helper",
                "--input_json", catGT_input_json,
                "--output_json", catGT_output_json
            ]
            print("Command:", command)
            subprocess.check_call(command)
            #result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            #print('Subprocess STDOUT:', result.stdout)
            #print('Subprocess STDERR:', result.stderr)
            #if result.returncode != 0:
                #raise subprocess.CalledProcessError(result.returncode, command)

        # Copy module input JSON file to data directory
        print('module_input_json: ', module_input_json)
        tmp = os.path.join(data_directory, f"{session_id}-input.json")
        print('paste dir: ', tmp)


        shutil.copy(module_input_json, os.path.join(data_directory, f"{session_id}-input.json"))

        for module in modules:
            output_json = os.path.join(json_directory, f"{session_id}-{module}-output.json")
            command = [
                sys.executable,
                "-W", "ignore",
                "-m", f"ecephys_spike_sorting.modules.{module}",
                "--input_json", module_input_json,
                "--output_json", output_json
            ]
            print("Command:", command)
            subprocess.check_call(command)
            #result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            #print('Subprocess STDOUT:', result.stdout)
            #print('Subprocess STDERR:', result.stderr)
            #if result.returncode != 0:
                #raise subprocess.CalledProcessError(result.returncode, command)

        log_from_json.addEntry(modules, json_directory, session_id, logFullPath)

