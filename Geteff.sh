#!/bin/bash   
export WORKDIR=$PWD
export SAMPLEDIR=/eos/cms/store/user/fiorendi/p5prime/effKDE
export CMSSWDIR=/afs/cern.ch/user/x/xuqin/work/cmssw/CMSSW_10_4_0/src
export HOME=/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/effnewphifinal                                                                                                                 
echo setting HOME to $HOME
echo setting WORKDIR to $WORKDIR
echo setting CMSSWDIR to $CMSSWDIR

cd $CMSSWDIR
source  /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scram runtime -sh`

if [ ! -r $SAMPLEDIR/2016/lmnr/newphi/effDataset_b0_2016.root ]; then
    echo $SAMPLEDIR/2016/lmnr/newphi/effDataset_b0_2016.root not found
    exit 1
fi

if [ ! -r $SAMPLEDIR/2017/lmnr/newphi/effDataset_b0_2017.root ]; then
    echo $SAMPLEDIR/2017/lmnr/newphi/effDataset_b0_2017.root not found
    exit 1
fi

if [ ! -r $SAMPLEDIR/2018/lmnr/newphi/effDataset_b0_2018.root ]; then
    echo $SAMPLEDIR/2018/lmnr/newphi/effDataset_b0_2018.root not found
    exit 1
fi
if [ ! -r $HOME/Geteff.cc ]; then
    echo $HOME/Geteff.cc not found
	 exit 1
fi

cd $WORKDIR

cp $HOME/Geteff.cc .
cp $SAMPLEDIR/2016/lmnr/newphi/effDataset* .
cp $SAMPLEDIR/2017/lmnr/newphi/effDataset* .
cp $SAMPLEDIR/2018/lmnr/newphi/effDataset* .


for i in 0 1 2 3 5 7 
	do 
		for k in 6 7 8
			do
				for p in 0 1
					do
						echo root -q -b './Geteff.cc('${k}','${i}',200, '${p}')'
						root -q -b './Geteff.cc('${k}','${i}',200, '${p}')'
					done
			done
	done
cp 1Deff* $HOME/
rm Geteff.cc
rm effDataset*
rm 1Deff*


