#!/bin/sh
pwd
ls
echo PATH is:
echo $PATH
echo LD_LIBRARY_PATH is:
echo $LD_LIBRARY_PATH
mkdir Results
condor_submit CondorSubmission.sub
