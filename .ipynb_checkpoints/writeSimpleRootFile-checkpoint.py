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
    parser.add_argument("-o", "--outDir", help="Output directory", default="./")
    parser.add_argument("-j", "--ncpu", help="Number of cores to use", default=4, type=int)
    parser.add_argument("-n", "--maxEvents", help="Number of events per bin", default=1000, type=str)
    parser.add_argument("-p", "--pixelAVdir", help="pixelAV directory", default="./pixelav/")
    ops = parser.parse_args()

    # get absolute path and check if outdir exists
    outDir = os.path.abspath(ops.outDir)
    if not os.path.isdir(outDir):
        os.makedirs(outDir)

    # ./minbias.exe <outFileName> <maxEvents> <pTHatMin> <pTHatMax>
    path_to_executable = "./bin/minbias.exe"
    # maxEvents = str(1000)
    # options_list = []
    fileList = ["testTree0", "testTree1", "testTree2"]
    commands = []
    for f in fileList:
        outFileName = f"/home/elizahoward/cmspix28-mc-sim/temp/{f}"
        
        # delphes to track list for pixelAV
        # python utils/delphesRootToPixelAvTrackList.py -i outdir/cmsMatch/10/minbias_0.30_0.40_GeV.root -o test.txt
        trackList = ["python3", "utils/delphesRootToPixelAvTrackList.py", "-i", f"{outFileName}.root", "-o", f"{outFileName}.txt"]
    
        # pixelAV
        pixelAV = [ops.pixelAVdir, "./bin/ppixelav2_list_trkpy_n_2f.exe", "1", f"{outFileName}.txt", f"{outFileName}.out", f"{outFileName}_seed"]

        # Write parquet file
        parquet = ["python3", "datagen.py", "-f", "{f}.out", "-t", f]
        
        # commands
        commands.append([(trackList, pixelAV),]) # weird formatting is because pool expects a tuple at input
        
    # List of CPU cores to use for parallel execution
    num_cores = multiprocessing.cpu_count() if ops.ncpu == -1 else ops.ncpu

    # Create a pool of processes to run in parallel
    pool = multiprocessing.Pool(num_cores)
    
    # Launch the executable N times in parallel with different options
    print(commands) # Anthony you are here need to make the multiprocess work with delphes tied in
    pool.starmap(run_commands, commands)
    
    # Close the pool of processes
    pool.close()
    pool.join()
