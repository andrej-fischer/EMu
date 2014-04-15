CC              = g++-mp-4.7
GSL		= /opt/local/lib/libgsl.a /opt/local/lib/libgslcblas.a
CC_FLAGS	= -Wall -O3 -DHAVE_INLINE -static-libgcc -static-libstdc++ -fopenmp
LD_FLAGS	= ${GSL} -fopenmp
MAIN            = EMu
MODEL		= MutSpecEM
PREPARE        	= EMu-prepare
OBJECTS         = $(MODEL).o $(MAIN).o
#
all:  $(OBJECTS) $(PREPARE).o
	$(CC) $(CC_FLAGS) $(OBJECTS) -o ../build/$(MAIN) $(LD_FLAGS)
	rm -f *.o
#
$(PREPARE).o: $(PREPARE).cpp
	$(CC) $(CC_FLAGS) -c $(PREPARE).cpp
	$(CC) $(CC_FLAGS) $(PREPARE).o -o ../build/$(PREPARE) $(LD_FLAGS)
#
$(MODEL).o: $(MODEL).cpp $(MODEL).h
	$(CC) $(CC_FLAGS) -c $(MODEL).cpp
#
$(MAIN).o: $(MAIN).cpp
	$(CC) $(CC_FLAGS) -c $(MAIN).cpp
#
clean:
	rm -f *.o
