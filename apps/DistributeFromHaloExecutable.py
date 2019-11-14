import multiprocessing
import os
import subprocess

import numpy as np

executable = '/home/gonayl/Chaste/build/projects/cellMorphogen/apps/FromHaloExecutable'

if not(os.path.isfile(executable)):
    raise Exception('Could not find executable: ' + executable)

number_of_simulations = 50

def main():
    run_simulations()

# Create a list of commands and pass them to separate processes
def run_simulations():

    # Make a list of calls to a Chaste executable
    command_list = []
    for epiepi in np.arange(4,6):
        for motility in np.arange(12,18):
            for epibnd in np.arange(8,12):
                command_list.append("nice -n 19 /home/gonayl/Chaste/build/projects/cellMorphogen/apps/FromHaloExecutable --E "+str(epiepi)+" --Mo "+str(motility)+" --Eb "+str(epibnd))

    # Use processes equal to the number of cpus available
    count = multiprocessing.cpu_count()

    print("Starting simulations with " + str(count) + " processes")

    # Generate a pool of workers
    pool = multiprocessing.Pool(processes=count)

    # Pass the list of bash commands to the pool
    pool.map_async(execute_command, command_list).get(86400)

# This is a helper function for run_simulation that runs bash commands in separate processes
def execute_command(cmd):
    return subprocess.call(cmd, shell=True)

if __name__ == "__main__":
    main()
