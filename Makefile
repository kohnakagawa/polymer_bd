#CXX = icpc
CXX = g++

CXXFLAGS = -std=c++11

ifeq ($(CXX), icpc)
RELEASE_FLAG = -O3 -xHOST -no-prec-div -Wall
endif

ifeq ($(CXX), g++)
RELEASE_FLAG = -O3 -Wall -ffast-math -funroll-loops
endif

DEBUG_FLAG = -O0 -g

CXXFLAGS += $(RELEASE_FLAG)
# CXXFLAGS += $(DEBUG_FLAG)

OBJECTS = main.o observer.o
TARGET  = polymer_bd.out
LOADLIBES = -lglfw -lGLU -L/usr/lib/nvidia-331-updates -lGL

.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LOADLIBES) -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) core.* *~
