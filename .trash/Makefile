############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Realistic Data MCMC
##
## make all PAR=yes to compile in parallel mode
## make all to compile in normal mode
## make execute PAR=yes to execute in parallel mode
## make execute (or just make) to execute in normal mode
##
## TODO : Put the execution in a script
###############################################################################


#************ COMPILER OPTIONS **************#

PAR = no                                              # Parallel mode or not
EXEC_NAME = RealisticDataMCMC                         # Executable name
CCFLAGSFORTRAN =                                      # Fortran compilation flags
LDFLAGS = -lgfortran -lm -lfftw3 #-I./lib -L./lib -lwavelet2s # C compilation flags
# EXECFLAGS2 = -d ./data/ -c ./confRun
# EXECFLAGS2 = -d $(REP) -c $(CONF) # To run in mode : make PAR=yes REP=./data3/ CONF=./data3/confRun

ifeq ($(PAR),yes)
CPP      = mpic++
CF90     = mpif90
CPPFLAGS = -O3 -Wall -I../lib -DPAR
EXEC = mpirun
EXECFLAGS = -np 30
else
CPP      = g++
CF90     = g++
CPPFLAGS = -O3 -Wall -I../lib
EXEC = $(EXEC_NAME)
EXECFLAGS = 
endif

#******* DO NOT CHANGE ANYTHING BELOW *******#

OBJ_DIR_NAME=obj
SRC_DIR_NAME=src
BIN_DIR_NAME=bin

PWD :=  $(shell pwd)
OBJ_DIR = $(PWD)/$(OBJ_DIR_NAME)
SRC_DIR = $(PWD)/$(SRC_DIR_NAME)
BIN_DIR = $(PWD)/$(BIN_DIR_NAME)

# List of source files 
SOURCESC = $(foreach nameC, $(SRC_DIR), $(wildcard $(nameC)/*.cpp)) #$(wildcard $(SRC_DIR)/*.cpp)
SOURCESF90 = $(foreach nameF, $(SRC_DIR), $(wildcard $(nameF)/*.f90)) # $(wildcard $(SRC_DIR)/*.f90)
# if $(SRC) contains multiple directories (don't work for the moment)
#SRC = $(foreach name, $(SRC_DIR), $(wildcard $(name)/*.cpp))

# List of object files :
NAME = $(basename $(notdir $(SOURCESC))) $(basename $(notdir $(SOURCESF90)))
OBJ = $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(NAME)))

# Execute
ifeq ($(PAR),yes)
execute: all
	$(EXEC) $(EXECFLAGS) $(BIN_DIR)/$(EXEC_NAME) $(EXECFLAGS2) 
else
execute: all
	$(BIN_DIR)/$(EXEC_NAME)
endif

# To make the executable
all: make_directory $(OBJ)
	cd $(SRC_DIR)/wavelib; make $(OBJ_DIR)
	cd $(SRC_DIR)/rgs; make $(OBJ_DIR)
	cd $(SRC_DIR)/FTeik; make $(OBJ_DIR)
	$(CPP) -o $(BIN_DIR)/$(EXEC) $(OBJ) $(LDFLAGS)

# Rule to make directory
make_directory: dirobj bin

# Create all .o from the .cpp found
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) -o $@ -c $< $(CPPFLAGS)

# Make the FTeik.o out of FTeik.f90
$(OBJ_DIR)/FTeik.o: $(SRC_DIR)/FTeik.f90
	$(CF90) -o $@ -c $< $(CCFLAGFORTRAN)

# .PHONY is usefull for example if it exists a file named "bin" in the current directory "make_directory: dirobj bin" would not work
.PHONY: dirobj bin

# Create the directory $(OBJ_DIR) if needed :
ifeq ($(strip $( $(wildcard $(OBJ_DIR)) ) ), )
dirobj: 
	@mkdir -p $(OBJ_DIR) # The @ mute the command
else
dirobj: 
endif

# Create the directory $(BIN_DIR) if needed :
ifeq ($(strip $( $(wildcard $(BIN_DIR)) ) ), )
bin: 
	@mkdir -p $(BIN_DIR)
else
bin: 
endif 

.PHONY: clean purge 

# Clean remove all .o (if you don't want to see the messages use : @rm -rf *.o -> only the error message will appear)
clean:
	rm -rf $(BIN_DIR)/$(EXEC) $(OBJ_DIR) $(SRC_DIR)/*~ *~

# Clean + remove executable
purge: clean
	rm -rf $(BIN_DIR) $(OBJ_DIR) OUTPUT_FILES

# Command to compile in sequential mode :
#  g++ -O3 -Wall -I../lib -c RealisticDataMCMC.cpp -o RealisticDataMCMC.o
#  g++ -O3 -Wall -I../lib -c generalFunctions.cpp -o generalFunctions.o 
#  g++ -O3 -Wall -I../lib -c rngs.cpp -o rngs.o 
#  g++ -O3 -Wall -I../lib -c rvgs.cpp -o rvgs.o 
#  g++ -O3 -Wall -I../lib -c filesAndControl.cpp -o filesAndControl.o 
#  g++ -O3 -Wall -I../lib -c functions.cpp -o functions.o
#  g++ -O3 -Wall -I../lib -c initializations.cpp -o initializations.o 
#  g++ -c FTeik.f90 -o FTeik.o 
#  g++ -o RealisticDataMCMC ./FTeik.o ./RealisticDataMCMC.o ./filesAndControl.o ./functions.o ./generalFunctions.o ./initializations.o ./rngs.o ./rvgs.o -lgfortran -lm -L../lib -lwavelet2s -lfftw3




