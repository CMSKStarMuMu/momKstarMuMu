#!/bin/bash                                                                                                                       
source /afs/cern.ch/user/x/xuqin/root6.12.sh

cd /afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/effnewphi

for i in 0 1 2 3 5 7 
	do 
		for k in 6 7 8
			do
				for p in 0 1
					do
						root -q -b './Geteff.cc('${k}','${i}',200, '${p}')'
					done
			done
	done

