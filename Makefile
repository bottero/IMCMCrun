############################ -*- Mode: Makefile -*- ###########################
## Makefile --- Realistic Data MCMC
##
## make all PAR=yes to compile in parallel mode
## make all to compile in normal mode
## make execute PAR=yes to execute in parallel mode
## make execute (or just make) to execute in normal mode
## !!! WARNING !! FOR GNU MPI COMPILERS CHANGE -lifcore TO -lgfortran
##
## TODO : Add debug options
###############################################################################

#************* PATH TO FFTW3 ****************#

FFTW3_INCLUDE=/usr/local/include     # Path to FFTW3 include directory (default /usr/local/include)
FFTW3_LIB=/usr/local/lib             # Path to FFTW3 lib directory (default /usr/local/lib)

#************ COMPILER OPTIONS **************#

PAR = no                                              # Parallel mode or not
EXEC_NAME = RealisticDataMCMC                         # Executable name

ifeq ($(PAR),yes)
CPP      = mpic++                    # C++ MPI compiler
CF90     = mpif90                    # Fortran MPI compiler
CPPFLAGS = -O3 -Wall -DPAR -I$(FFTW3_INCLUDE) -L$(FFTW3_LIB) -llib       # C++ compiler flags
F90FLAGS = -O3 #-check nobounds -ftz -implicitnone -warn all -nogen-interface # -gen-interfaces Fortran compilation flags
# for debugging: change -O3 -check nobounds to -O0 -check all -debug -g -fp-stack-check -traceback -ftrapuv
LDFLAGS = -lifcore -lm -lfftw3 -I$(FFTW3_INCLUDE) -L$(FFTW3_LIB)      # Link flags !!! WARNING !! FOR GNU MPI COMPILERS CHANGE -lifcore TO -lgfortran

else
CPP      = g++        #icc              # C++ compiler
CF90     = gfortran # ifort             # Fortran compiler
CPPFLAGS = -Wall -O3  -I$(FFTW3_INCLUDE) -L$(FFTW3_LIB) -llib  #-g # -ftz -traceback -ftrapuv -debug all # C++ compiler flags
F90FLAGS = -O3 #-g #-check all -debug -g -fp-stack-check -traceback -ftrapuv -implicitnone -gen-interfaces -warn all # Fortran compiler flags
# for debugging: change -O3 -check nobounds to -O0 -check all -debug -g -fp-stack-check -traceback -ftrapuv
LDFLAGS = -lifcore -lm -lfftw3  -I$(FFTW3_INCLUDE) -L$(FFTW3_LIB) #-lgfortran -lm -lfftw3      # Linker flags 
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
SOURCESC = $(wildcard $(SRC_DIR)/*.cpp)
SOURCESF90 = $(wildcard $(SRC_DIR)/*.f90)
# if $(SRC) contains multiple directories (don't work for the moment)
#SRC = $(foreach name, $(SRC_DIR), $(wildcard $(name)/*.cpp))

# List of object files :
NAME = $(basename $(notdir $(SOURCESC))) $(basename $(notdir $(SOURCESF90)))
OBJ = $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(NAME)))

# To make the executable
all: make_directory $(OBJ)
	$(CPP) -o $(BIN_DIR)/$(EXEC_NAME) $(OBJ) $(LDFLAGS)

# Rule to make directory
make_directory: dirobj bin

# Create all .o from the .cpp found
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) -o $@ -c $< $(CPPFLAGS)

# Make the fteik.o out of fteik.f90
$(OBJ_DIR)/fteik.o: $(SRC_DIR)/fteik.f90
	$(CF90) -o $@ -c $< $(F90FLAGS)

# .PHONY is usefull for example if it exists a file named "bin" in the current directory "make_directory: dirobj bin" would not work
.PHONY: dirobj bin

# Create the directory $(OBJ_DIR) if needed :
ifeq ($(strip $( $(wildcard $(OBJ_DIR)) ) ), )
dirobj: 
	@mkdir -p $(OBJ_DIR) # The @ mutes the command
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

# Clean remove bin and obj directories (if you don't want to see the messages use : @rm -rf *.o -> only the error message will appear)
clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR) *genmod*
	find . -name '*~' -print0 | xargs -0 rm

# Clean + remove executable
purge: clean
	rm -rf OUTPUT_FILES

# Command to compile in sequential mode :
#mkdir -p /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj
#mkdir -p /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/bin
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/filesAndControl.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/filesAndControl.cpp -O3 -Wall                 
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/functions.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/functions.cpp -O3 -Wall                 
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/generalFunctions.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/generalFunctions.cpp -O3 -Wall                 
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/initializations.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/initializations.cpp -O3 -Wall                 
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/RealisticDataMCMC.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/RealisticDataMCMC.cpp -O3 -Wall                 
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/rngs.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/rngs.cpp -O3 -Wall                 
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/rvgs.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/rvgs.cpp -O3 -Wall                 
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/wavelet2s.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/wavelet2s.cpp -O3 -Wall                 
#gfortran                   -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/fteik.o -c /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/src/fteik.f90 
#g++                        -o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/bin/RealisticDataMCMC                          /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/filesAndControl.o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/functions.o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/generalFunctions.o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/initializations.o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/RealisticDataMCMC.o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/rngs.o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/rvgs.o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/wavelet2s.o /home1/bottero/Desktop/DonneesRealistesLMA/IMCMCrun/obj/fteik.o -lgfortran -lm -lfftw3     

