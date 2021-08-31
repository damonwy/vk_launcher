
CXX = g++
CXXFLAGS = -Wall -g -std=c++14
DIR = .

INCLUDES = -I$(DIR)/include

SRCS = main.cpp #$(DIR)/source

BUILD = $(DIR)/build

OBJS = $(SRCS:.c=.o)

MAIN = unittest

all: $(MAIN)
	@echo Compiler named unittest has been compiled

.c.o: $(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

.PHONY: clean    
clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

build: $(MAIN)
	
# main.o: ./include/Vector3d.h ./include/Matrix4x4.h ./include/Quaternion.h
$(MAIN): $(OBJS) 
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS)

# DO NOT DELETE
