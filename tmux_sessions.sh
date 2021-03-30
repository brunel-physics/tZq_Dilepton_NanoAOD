#!/bin/bash
#tmux new-session -d -s tZq_mumu_2016 'make; make clean; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 0 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0'
#tmux new-session -d -s test_tmux 'vim test.txt'
tmux new-session -d -s tZq_mumu_2016 'vim Commands.sh; printf "#!/bin/bash \n make \n make clean \n ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 0 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0" >> Commands.sh; chmod+x Commands.sh; ./Commands.sh'
