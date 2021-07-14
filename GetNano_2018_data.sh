#!/bin/bash


ProcessNamesArray=("DoubleMuonRunB" "DoubleMuonRunC" "DoubleMuonRunD")
 			        

SampleStrings=("/DoubleMuon/Run2017B-Nano14Dec2018-v1/NANOAOD"
	       "/DoubleMuon/Run2017C-Nano14Dec2018-v1/NANOAOD"
	       "/DoubleMuon/Run2018D-Nano14Dec2018_ver2-v1/NANOAOD")


cd /nfs/data/eepgkkc/nanoAOD2018


for i in ${!ProcessNamesArray[@]}; do

	tmux_string="${ProcessNamesArray[i]}"

	mkdir $tmux_string
	cd $tmux_string
	
	echo ' '
	echo $tmux_string
 	echo ${SampleStrings[i]}
	pwd
	echo ' '	

	tmux new -d -s $tmux_string 'source ~/.bashrc; source /cvmfs/cms.cern.ch/cmsset_default.sh; getnano '${SampleStrings[i]}'; sleep 86400'
	cd ..

done
