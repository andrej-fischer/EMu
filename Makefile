CC              = g++-mp-4.7
INCL            = /Users/af7/local/include/
LIB             = /Users/af7/local/lib/
#
CC_FLAGS        = -Wall -O3 -I$(INCL) -DHAVE_INLINE -fopenmp
LD_FLAGS	= -L$(LIB) -lm -lgsl -lgslcblas -fopenmp
MAIN            = EMu
MODEL		= MutSpecEM
PREPARE        	= EMu-prepare
OBJECTS         = $(MODEL).o $(MAIN).o
#
all:  $(OBJECTS) $(PREPARE).o
	$(CC) $(CC_FLAGS) $(OBJECTS) -o $(MAIN) $(LD_FLAGS)
	rm -f *.o
#
$(PREPARE).o: $(PREPARE).cpp
	$(CC) $(CC_FLAGS) -c $(PREPARE).cpp
	$(CC) $(CC_FLAGS) $(PREPARE).o -o $(PREPARE) $(LD_FLAGS)
#
$(MODEL).o: $(MODEL).cpp $(MODEL).h
	$(CC) $(CC_FLAGS) -c $(MODEL).cpp
#
$(MAIN).o: $(MAIN).cpp
	$(CC) $(CC_FLAGS) -c $(MAIN).cpp
#
clean:
	rm -f *.o
