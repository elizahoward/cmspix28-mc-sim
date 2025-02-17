import sys
import numpy as np
import pandas as pd
import math

def split(outdir,tag,df1,df2,df3):

        df1.columns = df1.columns.astype(str)
        df2.columns = df2.columns.astype(str)
        df3.columns = df3.columns.astype(str)

        # unflipped, all charge                                                                         
        df1[df1['z-entry']==100].to_parquet(outdir+"labels"+tag+".parquet")
        df2[df1['z-entry']==100].to_parquet(outdir+"recon2D"+tag+".parquet")
        df3[df1['z-entry']==100].to_parquet(outdir+"recon3D"+tag+".parquet")

def parseFile(filein,tag,nevents=-1):

        with open(filein) as f:
                lines = f.readlines()

        header = lines.pop(0).strip()
        pixelstats = lines.pop(0).strip()

        print("Header: ", header)
        print("Pixelstats: ", pixelstats)

        readyToGetTruth = False
        readyToGetTimeSlice = False

        clusterctr = 0
        cluster_truth =[]
        timeslice = 0
        cur_slice = []
        cur_cluster = []
        events = []
        
        for line in lines:
                ## Start of the cluster
                if "<cluster>" in line:
                        readyToGetTruth = True
                        readyToGetTimeSlice = False
                        clusterctr += 1
                        
                        # Create an empty cluster
                        cur_cluster = []
                        timeslice = 0
                        # move to next line
                        continue

                # the line after cluster is the truth
                if readyToGetTruth:
                        cluster_truth.append(line.strip().split())
                        readyToGetTruth = False

                        # move to next line
                        continue

                ## Put cluster information into np array
                if "time slice" in line:
                        readyToGetTimeSlice = True
                        cur_slice = []
                        timeslice += 1
                        # move to next line
                        continue

                if readyToGetTimeSlice:
                        cur_row = line.strip().split()
                        cur_slice += [float(item) for item in cur_row]

                        # When you have all elements of the 2D image:
                        if len(cur_slice) == 13*21:
                                cur_cluster.append(cur_slice)

                        # When you have all time slices:
                        if len(cur_cluster) == 20:
                                events.append(cur_cluster)
                                readyToGetTimeSlice = False

        print("Number of clusters = ", len(cluster_truth))
        print("Number of events = ",len(events))
        print("Number of time slices in cluster = ", len(events[0]))

        arr_truth = np.array(cluster_truth)
        arr_events = np.array( events )

        return arr_events, arr_truth

def main():
        #tag = str(sys.argv[1])
        tag = "_0.00_0.10_GeV"
        inputdir = "/home/elizahoward/cmspix28-mc-sim/temp/"
        arr_events, arr_truth = parseFile(filein=inputdir+"minbias"+tag+".out",tag=tag)

        #truth quantities - all are dumped to DF                                                                                                                           
        df = pd.DataFrame(arr_truth, columns = ['x-entry', 'y-entry','z-entry', 'n_x', 'n_y', 'n_z', 'number_eh_pairs', 'y-local', 'pt'])
        cols = df.columns
        for col in cols:
                df[col] = df[col].astype(float)

        df['cotAlpha'] = df['n_x']/df['n_z']
        df['cotBeta'] = df['n_y']/df['n_z']

        sensor_thickness = 100 #um                                                          
        df['y-midplane'] = df['y-entry'] + df['cotBeta']*(sensor_thickness/2 - df['z-entry'])
        df['x-midplane'] = df['x-entry'] + df['cotAlpha']*(sensor_thickness/2 - df['z-entry'])

        print("The shape of the event array: ", arr_events.shape)
        print("The ndim of the event array: ", arr_events.ndim)
        print("The dtype of the event array: ", arr_events.dtype)
        print("The size of the event array: ", arr_events.size)

        df2 = {}
        df2list = []

        df3 = {}
        df3list = []

        for i, e in enumerate(arr_events):

                # Only last time slice
                df2list.append(np.array(e[-1]).flatten())

                # All time slices
                df3list.append(np.array(e).flatten())

                max_val = np.amax(e)

        df2 = pd.DataFrame(df2list)
        df3 = pd.DataFrame(df3list)  

        # split into flipped/unflipped, pos/neg charge
        split(inputdir, tag,df,df2,df3)

if __name__ == "__main__":
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
