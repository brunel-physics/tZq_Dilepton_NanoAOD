#!/bin/bash

echo ' '
echo -n "$bold Which region would you like, out of: SR, SBR, ZPlusJetsCR, ttbarCR or NoChi2Cut? $normal "
read Region
echo 'You have chosen the region: ' $Region
echo ' '


ProcessNamesArray=("tZq_scaleup" "tZq_scaledown") 

ChannelArray=("_ee" "_mumu")

YearArray=("2016")



if [[ $Region == 'SR' ]] 
then
	SRInt=1
	SBRInt=0
	zjcrInt=0
	ttcrInt=0

elif [[ $Region == 'SBR' ]] 
then
        SRInt=1
        SBRInt=1
        zjcrInt=0
        ttcrInt=0

elif [[ $Region == 'ZPlusJetsCR' ]] 
then
        SRInt=0
        SBRInt=0
        zjcrInt=1
        ttcrInt=0

elif [[ $Region == 'ttbarCR' ]] 
then
        SRInt=0
        SBRInt=0
        zjcrInt=0
        ttcrInt=1

elif [[ $Region == 'NoChi2Cut' ]]
then
        SRInt=0
        SBRInt=0
        zjcrInt=0
        ttcrInt=0
else 
	echo 'ERROR: Choose a region out of SR, SBR, ZPlusJetsCR, ttbarCR or NoChi2Cut. Exiting.'
	exit
	
fi

for i in ${!ChannelArray[@]}; do

	tmux_string_ScaleUp="${ProcessNamesArray[0]}${ChannelArray[i]}_${YearArray[0]}_$Region"

	echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 2 --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys 0 --channel '$i' --dcc 0'

	tmux new -d -s $tmux_string_ScaleUp './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 2 --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys 0 --channel '$i' --dcc 0; sleep 3600'

	tmux_string_ScaleDown="${ProcessNamesArray[1]}${ChannelArray[i]}_${YearArray[0]}_$Region"

	echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 3 --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys 0 --channel '$i' --dcc 0'

        tmux new -d -s $tmux_string_ScaleDown './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 3 --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys 0 --channel '$i' --dcc 0; sleep 3600'


done



echo "The tmux sessions for tZq scale up and scale down in the $Region (ee and mumu channels for 2016) are running."


