'''
Author: Anthony Badea
Date: 02/05/24
Purpose: Extract track parameters from delphes output root files and save to track list input for PixelAV
'''

import argparse
import uproot
import glob
import awkward as ak
import numpy as np
import os
import time

# Barrel radius
h=30

def getBetaNegative(y0, R):
    yc=(h**2 + y0**2 - (h**2*y0**2)/(h**2 + y0**2) - y0**4/(h**2 + y0**2) +  \
             y0*np.sqrt(-h**2*(h**2 + y0**2)*(h**2 - 4*R**2 + y0**2))/(h**2 + y0**2))/(2*h)
    beta = np.arccos((yc-h)/R)
    return beta

def getBetaPositive(y0, R):
    yc=(h**2 + y0**2 - (h**2*y0**2)/(h**2 + y0**2) - y0**4/(h**2 + y0**2) -  \
             y0*np.sqrt(-h**2*(h**2 + y0**2)*(h**2 - 4*R**2 + y0**2))/(h**2 + y0**2))/(2*h)
    beta = np.arccos((h-yc)/R)
    return beta

def getBeta(y0, R, q):
    return np.where(q<0, getBetaNegative(y0, R), getBetaPositive(y0, R))

def getYentryNegQ(R, beta):
    yc = R*np.cos(beta)+h
    y0 = -np.sqrt(R**2 - yc**2) + np.sqrt(-h**2 + R**2 + 2*h*yc - yc**2)
    return y0
    
def getYentryPosQ(R, beta):
    yc = -R*np.cos(beta)+h
    y0 = np.sqrt(R**2 - yc**2) - np.sqrt(-h**2 + R**2 + 2*h*yc - yc**2)
    return y0
    
def getYentry(R, q, beta):
    return np.where(q<0, getYentryNegQ(R, beta), getYentryPosQ(R, beta))
    
if __name__ == "__main__":
    
    # user options
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--inFileName", help="Input file name")
    parser.add_argument("-o", "--outFileName", help="Output file name", default="./")
    parser.add_argument("-p", "--float_precision", help="Float precision to save to track_list. ", default=5, type=int)
    parser.add_argument("-t", "--printFile", help="File to print to", default="")
    ops = parser.parse_args()

    i = 0
    while not os.path.exists(ops.inFileName):
        if i > 5:
            break
        time.sleep(5)
        i += 1
    
    # track list
    tracks = [] # cota cotb p flp localx localy pT

    # load the root files
    # files = "/local/d1/badea/tracker/smartpix/simulation/outdir/cmsMatch/10/*.root"
    tree = "Delphes"
    delphes_track_pt = []
    delphes_particle_pt = []
    branches = ["Track.PID", "Track.Charge", "Track.PT", "Track.P", "Track.EtaOuter", "Track.PhiOuter", "Track.XOuter", "Track.YOuter"]
    pionPID = 211 # plus/minus

    # for array in uproot.iterate(f"{files}:{tree}", branches):
    with uproot.open(ops.inFileName) as f:
        # load the branches
        temp = {}
        for branch in branches:
            temp[branch] = np.array(ak.flatten(f[tree][branch].array()))
        
        # selection
        cut = (abs(temp["Track.PID"])==pionPID)

        # apply selection
        for branch in branches:
            temp[branch] = temp[branch][cut]
        
        # Get particle track radius
        R = temp["Track.PT"]*5.36/(np.abs(temp["Track.Charge"])*1.60217663*3.8)*1000

        # The maximum and minimum possible entry points with respect to the whole pixel sensor:
        yentrymin=-16/2+0.08125
        yentrymax=16/2-0.08125

        # Using the max and min allowed entry points, get the max and min possible beta values
        betamin = getBeta(yentrymax, R, temp["Track.Charge"])
        betamax = getBeta(yentrymin, R, temp["Track.Charge"])

        # Randomly pick a beta within the allowed range
        beta=np.empty(len(temp["Track.PID"]))
        for i in range(len(temp["Track.PID"])):
            beta[i]=np.random.uniform(betamin[i], betamax[i])

        # Get the y entry point
        yentry=getYentry(R, temp["Track.Charge"], beta)

        # Shift yentry to ylocal
        localy=yentry+0.08125

        # Get angle of hit location with respect to barrel
        gamma0 = np.arctan2(temp["Track.YOuter"],temp["Track.XOuter"])

        # Get angle of sensor center with respect to barrel
        gamma=gamma0+np.arctan(yentry/30)

        # Define unit vector of track at tracker edge with respect to barrel
        theta=2*np.arctan(np.exp(-temp["Track.EtaOuter"]))
        phi=temp["Track.PhiOuter"]
        x=np.sin(theta)*np.cos(phi)
        y=np.sin(theta)*np.sin(phi)
        z=np.cos(theta)

        # Transform into rotated coordinate system (sensor coordinate system sort of)
        xp=x*np.cos(gamma)+y*np.sin(gamma)
        yp=-x*np.sin(gamma)+y*np.cos(gamma)

        alpha=np.arctan2(np.sqrt(xp**2+yp**2),z)
        
        p = temp["Track.P"] # [GeV]

        # zero: unflipped, 1: flipped
        flp = np.zeros(p.shape)
        
        # Since we are unflipped, we must adjust alpha and beta 
        cotb = 1./np.tan(2*np.pi-beta)
        cota = 1./np.tan(alpha+np.pi)

        localx = np.random.uniform(-8,8-1.05, size=len(temp["Track.PID"])) # [mm]
        
        pT = temp["Track.PT"] # [GeV]

        tracks.append([cota, cotb, p, flp, localx, localy, pT])

    tracks = np.concatenate(tracks,-1).T
    print("Tracks shape: ", tracks.shape)
    if ops.printFile != "":
        text_file = open(ops.printFile, "a")
        text_file.write("Track list information:\n")
        text_file.write(f"cota: {cota[0]}, cotb: {cotb[0]}, p: {p[0]}, localx: {localx[0]}, localy: {localy[0]}, pT: {pT[0]}\n\n")
        text_file.close()

    # save to file
    float_precision=4  
    with open(ops.outFileName, 'w') as file:
        for track in tracks:

            # set flp to an int
            track = list(track)
            track[3] = int(track[3])

            formatted_sublist = [f"{element:.{ops.float_precision}f}" if isinstance(element, float) else element for element in track]
            line = ' '.join(map(str, formatted_sublist)) + '\n'
            file.write(line)
        
