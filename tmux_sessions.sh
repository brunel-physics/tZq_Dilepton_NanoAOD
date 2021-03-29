#!/bin/bash

#2016

#nominal
#mumu
tmux new-session -d -s tZq_mumu_2016 'make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 0 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0'
