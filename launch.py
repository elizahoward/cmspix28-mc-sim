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
        if "pixelav" in command[0]:
            subprocess.run(command[1:], cwd=command[0])
        else:
            subprocess.run(command)
    

if __name__ == "__main__":

    # user options
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--outDir", help="Output directory", default="./TestDataUnflippedVertexCut")
    parser.add_argument("-j", "--ncpu", help="Number of cores to use", default=10, type=int)
    parser.add_argument("-n", "--maxEvents", help="Number of events per bin", default='200', type=str)
    parser.add_argument("-p", "--pixelAVdir", help="pixelAV directory", default="~/pixelav/")

    ops = parser.parse_args()

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

    pixelAVdir = os.path.expanduser(ops.pixelAVdir)
    semiparametricDir = os.path.expanduser(ops.semiparametricDir)

    # ./minbias.exe <outFileName> <maxEvents> <pTHatMin> <pTHatMax>
    path_to_executable = "./bin/minbias.exe"
    pt = np.linspace(0,2,21)
    # maxEvents = str(1000)
    # options_list = []
    commands = []
    for pTHatMin, pTHatMax in zip(pt[0:-1],pt[1:]):
        # fix numpy rounding
        pTHatMin = round(pTHatMin, 3)
        pTHatMax = round(pTHatMax, 3)
        # format
        tag = f"minbias_{pTHatMin:.2f}_{pTHatMax:.2f}_GeV"
        outFileName = f"{outDir}/{tag}"
        if pTHatMin != 1.9 and pTHatMin != 0:
            print(f"skipping{tag}")
            continue

        # pythia
        pythia = ["./bin/minbias.exe", outFileName, ops.maxEvents, str(pTHatMin), str(pTHatMax)]

        # delphes
        # card = os.path.join("/opt/delphes/cards", "delphes_card_CMS.tcl")
        card = "./delphes_card_CMS_ABEdit.tcl"
        delphes = ["/opt/delphes/DelphesHepMC3", card, outFileName+".root", outFileName+".hepmc"]

        # delphes to track list for pixelAV
        # python utils/delphesRootToPixelAvTrackList.py -i outdir/cmsMatch/10/minbias_0.30_0.40_GeV.root -o test.txt
        trackList = ["python3", "utils/delphesRootToPixelAvTrackList.py", "-i", f"{outFileName}.root", "-o", f"{outFileName}.txt"]

        # pixelAV
        # ../../pixelav/bin/ppixelav2_list_trkpy_n_2f.exe 1 outdir/cmsMatch/11/minbias_0.40_0.50_GeV.txt temp/minbias_0.40_0.50_GeV.out temp/seedfile
        pixelAV = [pixelAVdir, "./bin/ppixelav2_list_trkpy_n_2f.exe", "1", f"{outFileName}.txt", f"{outFileName}.out", f"{outFileName}_seed"]
        
        # Write parquet file
        parquet = ["python3", "./processing/datagen.py", "-f", f"{tag}.out", "-t", tag, "-d", outDir]

        # commands
        commands.append([(pythia, delphes, trackList, pixelAV, parquet,),]) # weird formatting is because pool expects a tuple at input 
        
    # List of CPU cores to use for parallel execution
    num_cores = multiprocessing.cpu_count() if ops.ncpu == -1 else ops.ncpu

    # Create a pool of processes to run in parallel
    pool = multiprocessing.Pool(num_cores)
    
    # Launch the executable N times in parallel with different options
    pool.starmap(run_commands, commands)
    
    # Close the pool of processes
    pool.close()
    pool.join()
