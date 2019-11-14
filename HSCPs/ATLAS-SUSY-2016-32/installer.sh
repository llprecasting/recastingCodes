#!/bin/sh

homeDIR="$( pwd )"


#Get pythia tarball
pythia="pythia8243.tgz"
#URL=http://home.thep.lu.se/~torbjorn/pythia8/$pythia
#echo "[installer] getting Pythia"; wget $URL 2>/dev/null || curl -O $URL; tar -zxf $pythia;
cd pythia8;
echo "Installing Pythia in pythia8226";
./configure --with-python-include=/usr/include/python2.7/ --with-root=$ROOTSYS --prefix=$homeDIR/pythia8
make -j4; make install;

cd $homeDIR
rm $pythia;


