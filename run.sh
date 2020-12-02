#!/bin/bash
export HOME=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/momentfinalnewphi/plugins
export WORKDIR=$PWD
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE
export CMSSWDIR=/afs/cern.ch/user/x/xuqin/work/cmssw/CMSSW_10_4_0/src

cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

if [ ! -r $SAMPLEDIR/2016/lmnr/newphi/effDataset_b0_2016.root]; then
    echo $SAMPLEDIR/2016/lmnr/newphi/effDataset_b0_2016.root not found
    exit 1
fi

if [ ! -r $SAMPLEDIR/2017/lmnr/newphi/effDataset_b0_2017.root]; then
    echo $SAMPLEDIR/2017/lmnr/newphi/effDataset_b0_2017.root not found
    exit 1
fi

if [ ! -r $SAMPLEDIR/2018/lmnr/newphi/effDataset_b0_2018.root]; then
    echo $SAMPLEDIR/2018/lmnr/newphi/effDataset_b0_2018.root not found
    exit 1
fi

cd $WORKDIR
if [ ! -r $HOME/Moment ]; then
    echo $HOME/Moment not found
	 exit 1
fi

cp $HOME/Moment .
cp $SAMPLEDIR/2016/lmnr/newphi/effDataset* .
cp $SAMPLEDIR/2017/lmnr/newphi/effDataset* .
cp $SAMPLEDIR/2018/lmnr/newphi/effDataset* .

echo setting HOME to $HOME
echo setting WORKDIR to $WORKDIR
echo setting CMSSWDIR to $CMSSWDIR


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
rm effDataset*
