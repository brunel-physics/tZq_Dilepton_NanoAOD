#Useful youtube videos about makefiles: 
#https://www.youtube.com/watch?v=_r7i5X0rXJk&ab_channel=PaulProgramming
#https://www.youtube.com/watch?v=NLhy05sdFQI&ab_channel=CalebCurry
#https://www.youtube.com/watch?v=6Gw1rNyTJWA&ab_channel=CalebCurry

tZq_executable: main.o tZq_Dilepton_Analysis.o
	g++ main.o tZq_Dilepton_Analysis.o -o tZq_executable

main.o: main.cpp
	g++ -c main.cpp 

tZq_Dilepton_Analysis.o: tZq_Dilepton_Analysis.cpp tZq_Dilepton_Analysis.h
	g++ -c tZq_Dilepton_Analysis.cpp

LIBRARY_PATH = 	-L$(shell root-config --libdir) \
		-Llib \
                -L/cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/lib \

LIBRARIES = 	$(shell root-config --libs) \

INCLUDE_PATH = 	-Iinclude  \
                -isystem/cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/include \
		-isystem$(shell root-config --incdir) \

CFLAGS = ${INCLUDE_PATH} -std=c++17 -MMD -MP -march=native \
		 -mtune=native -pipe -O2 -fPIC -m64 -pthread

CFLAGS += -Wall -Wextra -Wpedantic -Wcast-align -Wcast-qual \
			-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self \
			-Wlogical-op -Wmissing-declarations -Wmissing-include-dirs \
			-Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls \
			-Wshadow -Wsign-promo -Wstrict-null-sentinel \
			-Wstrict-overflow=5 -Wswitch-default -Wundef \
			-Wuseless-cast -Wzero-as-null-pointer-constant -Wduplicated-cond \
			-Wduplicated-branches -Wrestrict -Wnull-dereference -Wswitch-enum \
			-Wswitch-bool -Wswitch-unreachable -Wno-coverage-mismatch \
			-Wimplicit-fallthrough=5 -Wsync-nand -Wunknown-pragmas \
			-Wstringop-overflow=4 -Wstringop-truncation -Wsuggest-final-types \
			-Wsuggest-final-methods -Wsuggest-override -Walloc-zero -Walloca \
			-Wtrampolines -Wfloat-equal -Wunsafe-loop-optimizations \
			-Wplacement-new=2 -Wunused-macros -Wconditionally-supported \
			-Wsubobject-linkage -Wdate-time -Wextra-semi \
			-Wno-aggressive-loop-optimizations -Wpacked -Winvalid-pch  \
			-Wvector-operation-performance -Wdisabled-optimization \
			-Wstack-protector -Whsa -Wsuggest-attribute=const \
			-Wsuggest-attribute=pure -Wsuggest-attribute=noreturn \
            -Wsuggest-attribute=format -Wsuggest-attribute=cold


LINK_LIBRARY_FLAGS = -shared -rdynamic ${LIBRARY_PATH} ${LIBRARIES}

LINK_EXECUTABLE_FLAGS = -rdynamic ${LIBRARY_PATH} ${LIBRARIES} \
			-lTQZanalysisTools \
			-Wl,-Rlib,-R../lib,-R${PWD}/lib,--enable-new-dtags

${LIBRARY}: ${LIBRARY_OBJECT_FILES}
	${tZq_executable} ${LINK_LIBRARY_FLAGS} ${LIBRARY_OBJECT_FILES} -o $@

${LIBRARY_OBJECT_FILES}: obj/%.o : src/%.cpp
	mkdir -p {bin,obj,lib}
	${tZq_executable} -c ${CFLAGS}  $< -o $@

clean: 
	rm *.o tZq_executable

