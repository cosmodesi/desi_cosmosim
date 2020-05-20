.SUFFIXES:

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
##################################################################################
FLAGS    = -std=c++11 -fopenmp -fPIC
OPTIMIZE = -g3 -unroll 
CFLAGS   =
#-Wall -Wextra -pedantic -Wno-deprecated -Wno-write-strings  
LIBS = -lm -lgsl -lgslcblas
#########################################################################
##################################################################################	
# make will execute this rule by default
OBJS = get_dens_field.cpp
SRCS = $(OBJS:.o=.cpp)
TARGET = get_dens_field.exe
slics: $(TARGET)


$(TARGET): $(OBJS)
	$(CXX) $(FLAGS) $(SRCS) $(CFLAGS) $(OPTIMIZE) -o $(TARGET)


#########################################################################
clean:
	@echo "Cleaning:"
	rm -f  core.* *o *exe *~ 
#########################################################################
