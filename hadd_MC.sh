#!/bin/bash

bold=$(tput bold)
normal=$(tput sgr0)

echo ' '
echo -n "$bold What process would you like, out of: ZPlusJets, WPlusJets, ttbar, ttbarV, SingleTop, VVV or VVV? $normal "
read process
echo 'You have chosen the process: ' $process
echo ' '

echo ' '
echo -n "$bold What systematic would you like, out of: Nominal, PU_ScaleUp, PU_ScaleDown, BTag_ScaleUp, BTag_ScaleDown, JetSmearing_ScaleUp or JetSmearing_ScaleDown, JetResolution_ScaleUp, JetResolution_ScaleDown, LeptonEfficiencies_ScaleUp, LeptonEfficiencies_ScaleDown, PDF_ScaleUp, PDF_ScaleDown, MET_ScaleUp, MET_ScaleDown, ME_ScaleUp, ME_ScaleDown, isr_ScaleUp, isr_ScaleDown, fsr_ScaleUp or fsr_ScaleDown? $normal "
read systematic
echo 'You have chosen the systematic: ' $systematic
echo ' '

echo ' '
echo -n "$bold Which channel would you like, out of: ee, mumu or emu? $normal "
read channel
echo 'You have chosen the channel: ' $channel
echo ' '

echo ' '
echo -n "$bold Which year would you like, out of: 2016, 2017 or 2018? $normal "
read year
echo 'You have chosen the year: ' $year
echo ' '

echo ' '
echo -n "$bold Which region would you like, out of: SR, SBR, ZPlusJetsCR, ttbarCR or NoChi2Cut? $normal "
read region
echo 'You have chosen the region: ' $region
echo ' '

if [[ $region == 'SR' ]] 
then 
	RegionString='_SR___' 
elif [[ $region == 'SBR' ]] 
        
then 
	RegionString='_SR_SBR__'
elif [[ $region == 'NoChi2Cut' ]]
then
	RegionString='____'
else
	echo "ERROR: Need to add options for the z+jets and ttbar CRs. Exiting."
	exit
fi

hadd $process'_'$systematic'_'$channel'_'$region'_'$year'_Combined.root' 'Results_MC_'$process'_'*'_'$systematic'_'$channel'_'$RegionString'_'$year'.root'

