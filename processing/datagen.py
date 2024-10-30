import os
import datagensinglefile 
import argparse
import time

parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d", "--directory", help="Directory", default="/home/elizahoward/cmspix28-mc-sim/test/")
parser.add_argument("-f", "--filename", help="File name", default="testTree.out")
parser.add_argument("-t", "--tag", help="Tag", default="testTree")
ops = parser.parse_args()

i = 0
while not os.path.exists(f"{ops.directory}/{ops.filename}"):
    if i > 5:
        break
    time.sleep(5)
    i += 1

datagensinglefile.makeParquet(ops.filename, ops.tag, ops.directory)