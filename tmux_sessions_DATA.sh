#!/bin/bash

echo ' '
echo -n "$bold Which channel would you like, out of: ee, mumu or emu? $normal "
read Channel
echo 'You have chosen the channel: ' $Channel
echo ' '

echo ' '
echo -n "$bold Which year would you like, out of: 2016, 2017 or 2018? $normal "
read Year
echo 'You have chosen the year: ' $Year
echo ' '

echo ' '
echo -n "$bold Which region would you like, out of: SR, SBR, ZPlusJetsCR, ttbarCR or NoChi2Cut? $normal "
read Region
echo 'You have chosen the region: ' $Region
echo ' '

SystematicNamesArray=("Nominal"                       "PU_ScaleUp"             "PU_ScaleDown"           "BTag_ScaleUp"              "BTag_ScaleDown" 
		      "JetSmearing_ScaleUp"           "JetSmearing_ScaleDown"  "JetResolution_ScaleUp"  "JetResolution_ScaleDown"   "LeptonEfficiencies_ScaleUp"
		      "LeptonEfficiencies_ScaleDown"  "PDF_ScaleUp"  	       "PDF_ScaleDown"	        "ME_Up"			    "ME_Down"
		      "MET_Up"			      "MET_Down"	       "isr_up"		        "isr_down"		    "fsr_up"
		      "fsr_down")
	


ChannelArray=("ee" "mumu" "emu")

YearArray=("2016" "2017" "2018")

if [[ $Channel == 'ee' ]]
then
	ChannelInt=1

elif [[ $Channel == 'mumu' ]]
then
        ChannelInt=2

elif [[ $Channel == 'emu' ]]
then
        ChannelInt=3
else
	echo "ERROR: Choose a channel out of: ee, mumu or emu. Exiting."
	exit
fi



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




tmux_string="Data_${SystematicNamesArray[0]}_${ChannelArray[0]}_$Year"

k=$(($ChannelInt - 1))

echo './bin/fulleventselectionMain.exe --mc 0 -y '$Year' -p 159 --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys 0 --channel '$ChannelInt' --dcc 0'

tmux new -d -s $tmux_string './bin/fulleventselectionMain.exe --mc 0 -y '$Year' -p 159 --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys 0 --channel '$ChannelInt' --dcc 0; sleep 3600'




echo "The tmux sessions for data ($Channel channel for $Year in the $Region) are running."


