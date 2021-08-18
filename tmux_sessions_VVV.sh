#!/bin/bash

bold=$(tput bold)
normal=$(tput sgr0)

echo ' '
echo -n "$bold What systematic would you like, out of: Nominal, PU_ScaleUp, PU_ScaleDown, BTag_ScaleUp, BTag_ScaleDown, JetSmearing_ScaleUp, JetSmearing_ScaleDown, JetResolution_ScaleUp, JetResolution_ScaleDown, LeptonEfficiencies_ScaleUp, LeptonEfficiencies_ScaleDown, PDF_ScaleUp, PDF_ScaleDown, ME_Up, ME_Down, MET_Up, MET_Down, isr_up, isr_down, fsr_up, fsr_down? $normal "
read Systematic
echo 'You have chosen the systematic: ' $Systematic
echo ' '

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
 			       
ProcessNamesArray=("VVV_WWWTo4F"                        "VVV_WWWTo4F_ext"                       "VVV_WWZTo4F"
                   "VVV_WWZTo4F_ext"                    "VVV_WZZ"                               "VVV_WZZ_ext"
                   "VVV_ZZZ"                            "VVV_ZZZ_ext") 


SystematicNamesArray=("Nominal"                       "PU_ScaleUp"             "PU_ScaleDown"           "BTag_ScaleUp"              "BTag_ScaleDown" 
		      "JetSmearing_ScaleUp"           "JetSmearing_ScaleDown"  "JetResolution_ScaleUp"  "JetResolution_ScaleDown"   "LeptonEfficiencies_ScaleUp"
		      "LeptonEfficiencies_ScaleDown"  "PDF_ScaleUp"  	       "PDF_ScaleDown"	        "ME_Up"			    "ME_Down"
		      "MET_Up"			      "MET_Down"	       "isr_up"		        "isr_down"		    "fsr_up"
		      "fsr_down")
	


ChannelArray=("ee" "mumu" "emu")

YearArray=("2016" "2017" "2018")

for i in $( eval echo {0..${#SystematicNamesArray[@]}} )
do
	if [[ $Systematic == ${SystematicNamesArray[i]} ]]  
	then
		SystematicInt=$i
		echo "SystematicInt = $i."
	fi
done



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


for i in ${!ProcessNamesArray[@]}; do

	j=$(($i + 96))
	
	k=$(($ChannelInt - 1))

        tmux_string="${ProcessNamesArray[i]}_${SystematicNamesArray[$SystematicInt]}_${ChannelArray[$k]}_$Year"

        echo './bin/fulleventselectionMain.exe --mc 1 -y '$Year' -p '$j' --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys '$SystematicInt' --channel '$ChannelInt' --dcc 0'

	sleep_time=$(($i*180))

        echo $sleep_time

        #tmux new -d -s $tmux_string 'sleep '$sleep_time'; source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y '$Year' -p '$j' --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys '$SystematicInt' --channel '$ChannelInt' --dcc 0; sleep 86400'

	tmux new -d -s $tmux_string './bin/fulleventselectionMain.exe --mc 1 -y '$Year' -p '$j' --npl 0 --sr '$SRInt' --sbr '$SBRInt' --zjcr '$zjcrInt' --ttcr '$ttcrInt' --sys '$SystematicInt' --channel '$ChannelInt' --dcc 0; sleep 86400'

done



echo "The tmux sessions for VVV ($Channel channel for $Year, $Systematic, $Region) are running."


