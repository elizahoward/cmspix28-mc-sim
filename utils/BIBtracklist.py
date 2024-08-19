'''
Author: Anthony Badea
Date: 02/05/24
Purpose: Extract track parameters from delphes output root files and save to track list input for PixelAV
'''

import argparse
import glob
import os
import time
import numpy as np

if __name__ == "__main__":
    
    # user options
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--inFileName", help="Input file name", default = "/home/elizahoward/cmspix28-mc-sim/output/output.txt")
    parser.add_argument("-o", "--outFileName", help="Output file name", default = "/home/elizahoward/cmspix28-mc-sim/output/BIBtracks.txt")
    parser.add_argument("-p", "--float_precision", help="Float precision to save to track_list. ", default=5, type=int)
    ops = parser.parse_args()
    
    # track list
    tracks = [] # cota cotb p flp localx localy pT

    pionPID = 211 # plus/minus
    electronPID = 11

    lines = tuple(open(ops.inFileName, 'r'))
    lines = lines[1:] # cut header line
    for line in lines:
        cota, cotb, p, flp, localx, localy, pT, hittime, pid = list(map(float, line.split(" ")))
        if abs(pid) == electronPID or abs(pid) == pionPID:
            tracks.append([cota, cotb, p, flp, localx, localy, pT])
   
    tracks = np.array(tracks)
    print("Tracks shape: ", tracks.shape)
    
    # save to file
    with open(ops.outFileName, 'w') as file:
        for track in tracks:

            # set flp to an int
            track = list(track)
            track[3] = int(track[3])

            formatted_sublist = [f"{element:.{ops.float_precision}f}" if isinstance(element, float) else element for element in track]
            line = ' '.join(map(str, formatted_sublist)) + '\n'
            file.write(line)
        
