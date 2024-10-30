'''
Author: Anthony Badea
Date: 01/31/24
'''

import subprocess
import multiprocessing
import numpy as np
import os
import argparse

def run_executable(executable_path, options):
    command = [executable_path] + options
    subprocess.run(command)

def run_commands(commands):
    for command in commands:
        print(command)
        if "pixelav" in command[0]:
            subprocess.run(command[1:], cwd=command[0])
        else:
            subprocess.run(command)
    

if __name__ == "__main__":

    # user options
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outDir", help="Output directory", default="./MuC_Output_BfieldTest")
    parser.add_argument("-j", "--ncpu", help="Number of cores to use", default=10, type=int)
    parser.add_argument("-p", "--pixelAVdir", help="pixelAV directory", default="~/pixelav/")
    ops = parser.parse_args()

    # get absolute path for semiparametric directory
    semiparametricDir = os.path.expanduser(ops.semiparametricDir)
    pixelAVdir = os.path.expanduser(ops.pixelAVdir)

    # get absolute path and check if outdir exists
    outDir = os.path.abspath(ops.outDir)
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    else:
        # Empty folder
        files = os.listdir(outDir)
        for f in files:
            if os.path.isfile(f"{outDir}/{f}"):
                os.remove(f"{outDir}/{f}")

    commands = []

    folders = ["/home/elizahoward/produceSmartPixMuC/BIB_tracklists", "/home/elizahoward/produceSmartPixMuC/signal_tracklists"] 

    for folder in folders:

        if "BIB" in folder:
            tag0 = "bib"
        else:
            tag0 = "sig"
        
        i = 0

        for tracklist in os.listdir(folder):

            if tag0 == 'sig' and i > 76:
                break
            
            tag = f"{tag0}{i}"

            outFileName = f"{outDir}/{tag}"

            # Run pixelAV
            pixelAV = [pixelAVdir, "./bin/ppixelav2_list_trkpy_n_2f.exe", "1", f"{folder}/{tracklist}", f"{outFileName}.out", f"{outFileName}_seed"]

            # Write parquet file
            parquet = ["python3", "./processing/datagen.py", "-f", f"{tag}.out", "-t", tag, "-d", outDir]

            # commands
            commands.append([(pixelAV, parquet,),]) # weird formatting is because pool expects a tuple at input

            i += 1
        
    # List of CPU cores to use for parallel execution
    num_cores = multiprocessing.cpu_count() if ops.ncpu == -1 else ops.ncpu

    # Create a pool of processes to run in parallel
    pool = multiprocessing.Pool(num_cores)
    
    # Launch the executable N times in parallel with different options
    # pool.starmap(run_executable, [(path_to_executable, options) for options in options_list])
    print(commands) # Anthony you are here need to make the multiprocess work with delphes tied in
    pool.starmap(run_commands, commands)
    
    # Close the pool of processes
    pool.close()
    pool.join()
