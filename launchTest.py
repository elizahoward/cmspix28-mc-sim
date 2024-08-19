'''
Author: Anthony Badea
Date: 01/31/24
'''

import subprocess
import multiprocessing
import numpy as np
import os
import argparse
import time
#import glob

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
    parser.add_argument("-i", "--inputFile", help="Input root file", default="./temp/minbias_0.40_0.50_GeV.root")
    parser.add_argument("-o", "--outDir", help="Output directory", default="./test/")
    parser.add_argument("-j", "--ncpu", help="Number of cores to use", default=4, type=int)
    parser.add_argument("-p", "--pixelAVdir", help="pixelAV directory", default="./pixelav/")
    parser.add_argument("-s", "--semiparametricDir", help="semiparametric directory", default="~/semiparametric")
    ops = parser.parse_args()

    # get absolute path for semiparametric directory
    semiparametricDir = os.path.expanduser(ops.semiparametricDir)

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

    # ./minbias.exe <outFileName> <maxEvents> <pTHatMin> <pTHatMax>
    path_to_executable = "./bin/minbias.exe"
    indexList = [0,2,3]
        
    for i in indexList:
        printFile = f"{outDir}/Output{i}.txt"
    
        tag = f"test{i}"
        outFileName = f"{outDir}/{tag}"

        print(outFileName, printFile)
        # Write the root file
        writeRootFile = ["python3", "writeSimpleRootFile.py", "-i", str(i), "-o", f"{outFileName}.root", "-p", printFile]
        
        # delphes to track list for pixelAV
        # python utils/delphesRootToPixelAvTrackList.py -i outdir/cmsMatch/10/minbias_0.30_0.40_GeV.root -o test.txt
        trackList = ["python3", "utils/delphesRootToPixelAvTrackList.py", "-i", f"{outFileName}.root", "-o", f"{outFileName}.txt", "-t", printFile]
    
        # pixelAV
        pixelAV = [ops.pixelAVdir, "./bin/ppixelav2_list_trkpy_n_2f.exe", "1", f"{outFileName}.txt", f"{outFileName}.out", f"{outFileName}_seed", printFile]

        # Write parquet file
        parquet = ["python3", f"{semiparametricDir}/processing/datagen.py", "-f", f"{tag}.out", "-t", tag, "-d", outDir, "-p", printFile]
        
        # Make some plots
        makePlots = ["python3", f"{semiparametricDir}/plotting/makePlots.py", "-t", tag, "-d", outDir]


        # commands
        commands= [[(writeRootFile, trackList, pixelAV, parquet, makePlots),]] # weird formatting is because pool expects a tuple at input
        
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

        
