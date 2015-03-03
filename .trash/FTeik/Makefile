############################ -*- Mode: Makefile -*- ###########################
## Makefile --- wavelib
## Called by the main makefile.
## make OBJ_DIR=obj
###############################################################################

#************ COMPILER OPTIONS **************#

CPP=g++                    # fortran Compiler 

#******* DO NOT CHANGE ANYTHING BELOW *******#

SRC_DIR = .
CPPFLAGS=-W -Wall -Wextra -ansi -pedantic
LDFLAGS=	

# List of source files 
SRC = $(wildcard $(SRC_DIR)/*.f90)

# List of object files :
NAME = $(basename $(notdir $(SRC)))
OBJ = $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(NAME)))

all: $(OBJ)

# Create all .o from the .cpp found
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) -o $@ -c $< $(CPPFLAGS)
	
