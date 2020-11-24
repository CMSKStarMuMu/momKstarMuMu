#!/bin/bash
export HOME=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/momentnewphi/plugins
export WORKDIR=$PWD
cd /afs/cern.ch/user/x/xuqin
source root6.12.sh
cd $WORKDIR
if [ ! -r $HOME/Moment ]; then
    echo $HOME/Moment not found
	 exit 1
fi
cp $HOME/Moment .
echo setting HOME to $HOME
echo setting WORKDIR to $WORKDIR
echo ./Moment gen
./Moment gen
for i in 0 1
    do
        echo ./Moment reco good $i
        ./Moment reco good $i
        echo ./Moment reco mis $i
        ./Moment reco mis $i
        echo ./Moment reco mix $i
        ./Moment reco mix $i
    done

rm Moment
