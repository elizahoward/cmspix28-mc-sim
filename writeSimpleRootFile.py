import argparse, os, uproot, argparse
import awkward as ak
import numpy as np

tree = "Delphes"
branches = ["Track.PID", "Track.PT", "Track.P", "Track.EtaOuter", "Track.PhiOuter", "Track.XOuter", "Track.YOuter"]
pionPID = 211 # plus/minus

def reatRootFile(filename, i):
    # load the branches
    with uproot.open(filename) as f:
        
        temp = {}
        for branch in branches:
            temp[branch] = np.array(ak.flatten(f[tree][branch].array()))
        
        # selection
        cut = (abs(temp["Track.PID"])==pionPID)

        # apply selection (remove everything except pions)
        for branch in branches:
            temp[branch] = temp[branch][cut]
        
        # track properties
        # based on the image here https://github.com/kdp-lab/pixelav/blob/ppixelav2v2/ppixelav2_operating_inst.pdf
        phi = ak.Array([[temp["Track.PhiOuter"][i]]]) 
        eta = ak.Array([[temp["Track.EtaOuter"][i]]]) 
        p = ak.Array([[temp["Track.P"][i]]]) # [GeV]
        localx = ak.Array([[temp["Track.XOuter"][i]]]) # [mm]
        localy = ak.Array([[temp["Track.YOuter"][i]]]) # [mm]
        pT = ak.Array([[temp["Track.PT"][i]]]) # [GeV]
        pid = ak.Array([[temp["Track.PID"][i]]])
        
    return phi, eta, p, localx, localy, pT, pid

def writeSimpleRootFile(inputFile, outputFile, i, printFile):
    phi, eta, p, localx, localy, pT, pid = reatRootFile(inputFile, i)
    
    # Write root file with just the one entry

    with uproot.recreate(outputFile) as f:
        f["Delphes"] = {"Track.PID":pid, "Track.PT":pT, "Track.P":p, "Track.EtaOuter":eta, "Track.PhiOuter":phi, "Track.XOuter":localx, "Track.YOuter":localy}
    
    phi, eta, p, localx, localy, pT, pid = reatRootFile(inputFile, 0)

    if ops.printFile != "":
        text_file = open(ops.printFile, "a")
        text_file.write("Delphes Output:")
        text_file.write(f"\nPID: {pid[0][0]}, phi: {phi[0][0]}, eta: {eta[0][0]}, p: {p[0][0]}, XOuter: {localx[0][0]}, YOuter: {localy[0][0]}, pT: {pT[0][0]}\n")
        text_file.write(f"\u221A(XOuter\u00B2+YOuter\u00B2) = {np.sqrt(localx[0][0]**2+localy[0][0]**2)} mm\n")
        text_file.write(f"cot(\u03B1) = cot(\u03D5+\u03C0) = {1/np.tan(phi[0][0]+np.pi)}\n")
        text_file.write(f"cot(\u03B2) = cot(cos\u207B\u00B9(tanh(\u03B7))+\u03C0) = {1/np.tan(np.arccos(np.tanh(eta[0][0]))+np.pi)}\n\n")
        text_file.close()


parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--index", help="Index for root file", default=0, type=int)
parser.add_argument("-f", "--inputFile", help="Input file path", default="./temp/minbias_0.40_0.50_GeV.root")
parser.add_argument("-o", "--outputFile", help="Input file path", default="")
parser.add_argument("-p", "--printFile", help="Text file path", default="")
ops = parser.parse_args()

# get absolute path for input file
inputFile = os.path.abspath(ops.inputFile)

writeSimpleRootFile(inputFile, ops.outputFile, ops.index, ops.printFile)
