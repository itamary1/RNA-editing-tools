# Function for wrapping subprocesses opening
import subprocess
import traceback


def execute(command, error_message):
    try:
        print "RUNNING: " +  command
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError:
        print "ERROR: " + error_message
        print traceback.format_exc()
        exit(1)
