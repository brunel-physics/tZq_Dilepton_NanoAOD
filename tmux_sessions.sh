#!/bin/bash
#tmux new-session -d -s tZq_mumu_2016 'make; make clean; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 0 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0'
#tmux new-session -d -s test_tmux 'vim test.txt'
vim Commands.sh
chmod +x Commands.sh
tmux new-session -d -s tZq_mumu_2016 'printf "#!/bin/bash\nmake\nmake clean\n./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 0 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0" >> Commands.sh'
