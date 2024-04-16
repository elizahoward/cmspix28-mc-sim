# simulations for smart pixels

# setup environment
```
setupATLAS
lsetup "python 3.9.18-x86_64-centos7"
lsetup "root 6.28.08-x86_64-centos7-gcc11-opt"
```

# first get pythia 8.307
```
wget https://pythia.org/download/pythia83/pythia8307.tgz
tar -xvf pythia8307.tgz
cd pythia8307
make
````

# update this path to your pythia path in Makefile.inc
`PREFIX=/home/abadea/pythia8307`

# setup this repo and run test
```
make all
./bin/minbias.exe cards/minbias.cmnd myTree.root 100
```