#!/usr/bin/env python3.8

# RUNNING_STEP = "Process: %(ps_name)s; Running Step: %(step)s"
# ERROR_STEP_MSG = "Process: %(ps_name)s; Going To Error Step: %(step)s"
import sys
import subprocess
import logging
import traceback

RUNNING_CMD = "Running Command: %(cmd)s"
ERROR_MSG = "Command returned %(ret_code)s: %(cmd)s"

def execute(command, input_str=None, log=True):
    try:
        # info
        logging.debug(RUNNING_CMD % {'cmd': command}) if log else print("RUNNING: " +  command)
        sys.stdout.flush() # clears stdout
        if input_str:
            p = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               input=input_str, universal_newlines=True, text=True)
        else:
            p = subprocess.run(command, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
        # return output
        return p.stdout
    except subprocess.CalledProcessError as e:
        logging.error(ERROR_MSG % {'cmd': command, 'ret_code': e.returncode}) if log else print("ERROR")
        # print("ERROR:",  p.stderr)
        print(traceback.format_exc())
        exit(1)