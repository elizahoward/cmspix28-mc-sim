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

def getBetaPositive(y0, R):
    xc=(h**2 + y0**2 - (h**2*y0**2)/(h**2 + y0**2) - y0**4/(h**2 + y0**2) +  \
             y0*np.sqrt(-h**2*(h**2 + y0**2)*(h**2 - 4*R**2 + y0**2))/(h**2 + y0**2))/(2*h)
    beta = np.arccos((xc-h)/R)
    return beta

def getBetaNegative(y0, R):
    xc=(h**2 + y0**2 - (h**2*y0**2)/(h**2 + y0**2) - y0**4/(h**2 + y0**2) -  \
             y0*np.sqrt(-h**2*(h**2 + y0**2)*(h**2 - 4*R**2 + y0**2))/(h**2 + y0**2))/(2*h)
    beta = np.arccos((h-xc)/R)
    return beta

def getBeta(y0, R, q):
    return np.where(q<0, getBetaNegative(y0, R), getBetaPositive(y0, R))

def getYentryPosQ(R, beta):
    xc = R*np.cos(beta)+h
    yc = -np.sqrt(R**2 - xc**2)
    y0 = np.sqrt(R**2 - xc**2) - np.sqrt(-h**2 + R**2 + 2*h*xc - xc**2)
    return y0, xc, yc
    
def getYentryNegQ(R, beta):
    xc = -R*np.cos(beta)+h
    yc = np.sqrt(R**2 - xc**2)
    y0 = - np.sqrt(R**2 - xc**2) + np.sqrt(-h**2 + R**2 + 2*h*xc - xc**2)
    return y0, xc, yc
    
def getYentry(R, q, beta):
    return np.where(q<0, getYentryNegQ(R, beta), getYentryPosQ(R, beta))


def getCrossingPoint2(xc,yc,R):
    y = (1/(2*(xc**2 + yc**2)))*(h**2*yc - R**2*yc + yc**3 + yc*xc**2 - np.sqrt(-h**4*xc**2 + 2*h**2*R**2*xc**2 - R**4*xc**2 + 2*h**2*yc**2*xc**2 + 2*R**2*xc**2*yc**2 - yc**4*xc**2 + 2*h**2*xc**4 + 2*R**2*xc**4 - 2*yc**2*xc**4 - xc**6))
    x = -(1/(2*xc))*(-h**2 + R**2 - xc**2 - yc**2 + (h**2*yc**2)/(xc**2 + yc**2) - (R**2*yc**2)/(xc**2 + yc**2) + yc**4/(xc**2 + yc**2) + (xc**2*yc**2)/(xc**2 + yc**2) - (1/(xc**2 + yc**2))*yc*np.sqrt(-h**4*xc**2 + 2*h**2*R**2*xc**2 - R**4*xc**2 + 2*h**2*xc**2*yc**2 + 2*R**2*xc**2*yc**2 - yc**4*xc**2 + 2*h**2*xc**4 + 2*R**2*xc**4 - 2*yc**2*xc**4 - xc**6))
    return x,y


def getCrossingPoint1(xc,yc,R):
    y = (1/(2*(xc**2 + yc**2)))*(h**2*yc - R**2*yc + yc**3 + yc*xc**2 + np.sqrt(-h**4*xc**2 + 2*h**2*R**2*xc**2 - R**4*xc**2 + 2*h**2*yc**2*xc**2 + 2*R**2*xc**2*yc**2 - yc**4*xc**2 + 2*h**2*xc**4 + 2*R**2*xc**4 - 2*yc**2*xc**4 - xc**6))
    x = -(1/(2*xc))*(-h**2 + R**2 - xc**2 - yc**2 + (h**2*yc**2)/(xc**2 + yc**2) - (R**2*yc**2)/(xc**2 + yc**2) + yc**4/(xc**2 + yc**2) + (xc**2*yc**2)/(xc**2 + yc**2) + (1/(xc**2 + yc**2))*yc*np.sqrt(-h**4*xc**2 + 2*h**2*R**2*xc**2 - R**4*xc**2 + 2*h**2*xc**2*yc**2 + 2*R**2*xc**2*yc**2 - yc**4*xc**2 + 2*h**2*xc**4 + 2*R**2*xc**4 - 2*yc**2*xc**4 - xc**6))
    return x,y

def getCrossingPointNegQ(xc, yc, R):
    return np.where(xc<0, getCrossingPoint1(xc,yc,R), getCrossingPoint2(xc,yc,R))

def getCrossingPointPosQ(xc, yc, R):
    return np.where(xc>0, getCrossingPoint1(xc,yc,R), getCrossingPoint2(xc,yc,R))

def getCrossingPoint(xc, yc, R, q):
    x,y=np.where(q<0, getCrossingPointNegQ(xc,yc,R), getCrossingPointPosQ(xc,yc,R))
    return x,y

def getGamma(xOuter, yOuter, xc, yc, R, q):
    x,y=getCrossingPoint(xc,yc,R,q)
    gamma0=np.arctan2(yOuter,xOuter)
    deltaGamma=np.arctan2(y,x)
    gamma = gamma0-deltaGamma
    return gamma

def getInfo(temp, yentry, beta):
    # Get angle of hit location with respect to barrel
    gamma0 = np.arctan2(temp["Track.YOuter"],temp["Track.XOuter"])

    # Get angle of sensor center with respect to barrel
    gamma=gamma0+np.arctan(yentry/30)

    # Define unit vector of track at tracker edge with respect to barrel
    theta=2*np.arctan(np.exp(-temp["Track.EtaOuter"])) # assume eta does not change between r = 30 mm and where it hits the sensor
    phi=beta-np.pi/2+gamma # determine phi based on beta and gamma (shift beta back into barrel coordinate system)
    x=np.sin(theta)*np.cos(phi)
    y=np.sin(theta)*np.sin(phi)
    z=np.cos(theta) 

    # Transform into rotated coordinate system (sensor coordinate system sort of)
    # xp=x*np.cos(np.pi/2-gamma)-y*np.sin(np.pi/2-gamma)
    yp=x*np.sin(np.pi/2-gamma)+y*np.cos(np.pi/2-gamma)

    alpha=np.arctan2(yp,z)

    if ops.flp == 0:
        # For the unflipped geometry, we must adjust alpha and beta 
        cotb = 1./np.tan(beta)
        cota = 1./np.tan(alpha)
    else:
        cotb = 1./np.tan(beta)
        cota = 1./np.tan(alpha)

    localx = np.random.uniform(-8,8-1.05, size=len(temp["Track.PID"])) # [mm]
    
    return cota, cotb, localx


def assumeCircle(temp):
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
    yentry, xc, yc = getYentry(R, temp["Track.Charge"], beta)

    # Shift yentry to ylocal
    localy=yentry+0.08125

    # Round ylocal to the nearest pixel
    localy /= 0.0125
    localy = np.rint(localy)
    localy *= 0.0125

    gamma = getGamma(temp["Track.XOuter"],temp["Track.YOuter"], xc, yc, R, temp["Track.Charge"])

    # Define unit vector of track at tracker edge with respect to barrel
    theta=2*np.arctan(np.exp(-temp["Track.EtaOuter"])) # assume eta does not change between r = 30 mm and where it hits the sensor
    phi=beta-np.pi/2+gamma # determine phi based on beta and gamma (shift beta back into barrel coordinate system)
    x=np.sin(theta)*np.cos(phi)
    y=np.sin(theta)*np.sin(phi)
    z=np.cos(theta) 

    # Transform into rotated coordinate system (sensor coordinate system sort of)
    # xp=x*np.cos(np.pi/2-gamma)-y*np.sin(np.pi/2-gamma)
    yp=x*np.sin(np.pi/2-gamma)+y*np.cos(np.pi/2-gamma)

    alpha=np.arctan2(yp,z)

    if ops.flp == 0:
        # For the unflipped geometry, we must adjust alpha and beta 
        cotb = 1./np.tan(beta+np.pi)
        cota = 1./np.tan(2*np.pi-alpha)
    else:
        cotb = 1./np.tan(beta)
        cota = 1./np.tan(alpha)

    # Randomly generate xlocal
    localx = np.random.uniform(-8,8-1.05, size=len(temp["Track.PID"])) # [mm]

    # Round xlocal to the nearest pixel
    localx /= 0.05
    localx = np.rint(localx)
    localx *= 0.05

    return cota, cotb, localx, localy


def assumeStraight(temp):
    # Randomly generate beta between minimum and maximum allowed values
    betamin = 1.31272
    betamax = 1.82887
    beta = np.random.uniform(betamin, betamax,size=len(temp["Track.PID"]))

    # Determine yentry and ylocal based on beta
    yentry = 30/np.tan(beta)
    localy=yentry+0.08125

    cota, cotb, localx = getInfo(temp, yentry, beta)

    return cota, cotb, localx, localy


def originalMethod(temp):
    # based on the image here https://github.com/kdp-lab/pixelav/blob/ppixelav2v2/ppixelav2_operating_inst.pdf
    # phi = beta - pi -> cot(beta) = cot(phi+pi) = cot(phi) = 1/tan(phi)
    cotb = 1./np.tan(temp["Track.PhiOuter"]) 
    # theta = alpha - pi -> cot(alpha) = cot(theta+pi) = cot(theta) = 1/tan(theta) where theta = arccos(tanh(eta))
    cota = 1./np.tan(np.arccos(np.tanh(temp["Track.EtaOuter"]))) 
    
    localx = temp["Track.XOuter"] # [mm]
    localy = temp["Track.YOuter"] # [mm]

    return cota, cotb, localx, localy
    

if __name__ == "__main__":
    
    # user options
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--inFileName", help="Input file name")
    parser.add_argument("-o", "--outFileName", help="Output file name", default="./")
    parser.add_argument("-p", "--float_precision", help="Float precision to save to track_list. ", default=5, type=int)
    parser.add_argument("-a", "--pathAssumption", help="Assumption about particle path (cirlce: 0, straight: 1, original: 2)", default=0, type=int)
    parser.add_argument("-f", "--flp", help="Flipped (1) vs Unflipped (0)", default=1, type=int)
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
    electronPID = 11 # plus/minus

    # for array in uproot.iterate(f"{files}:{tree}", branches):
    with uproot.open(ops.inFileName) as f:
        # load the branches
        temp = {}
        for branch in branches:
            temp[branch] = np.array(ak.flatten(f[tree][branch].array()))
        
        # selection
        cut = abs(temp["Track.PID"])==pionPID 

        # apply selection
        for branch in branches:
            temp[branch] = temp[branch][cut]
        
        p = temp["Track.P"] # [GeV]

        # zero: unflipped, 1: flipped
        flp = np.zeros(p.shape)

        if ops.flp == 1:
            flp += 1
        
        pT = temp["Track.PT"] # [GeV]

        if ops.pathAssumption == 0:
            cota, cotb, localx, localy = assumeCircle(temp)
        elif ops.pathAssumption == 1:
            cota, cotb, localx, localy = assumeStraight(temp)
        else:
            cota, cotb, localx, localy = originalMethod(temp)

        pid = temp["Track.PID"]

        tracks.append([cota, cotb, p, flp, localx, localy, pT, pid])

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
        
