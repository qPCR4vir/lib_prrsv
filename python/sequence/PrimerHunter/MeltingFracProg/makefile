#CFLAGS      = -D_output_alignment -fno-weak
CFLAGS      = -fno-weak
LFLAGS      =
#CFLAGS     = -g -gall -ms -xtaso -stdnew -global_array_new -std ms -Zp1 -misalign -assume noaligned_objects -nomember_alignment -protect_headers all
#LFLAGS     = -g -gall -ms -xtaso -stdnew -global_array_new -std ms -Zp1 -misalign -assume noaligned_objects -nomember_alignment -protect_headers all
LIBS       = -lm
EXTRALIBS  =

PURIFY        = purify -g++ -best-effort
#CCC           = /zpr/gcc/2.8/bin/g++
#CCC           = cxx -O5
#CCC           = g++ -g -Wno-deprecated
CCC           = g++ -O3 -Wno-deprecated

INCLUDE       = -I.
LIBDIRS       =

COMPILE.cc = $(CCC) $(CFLAGS) $(INCLUDE) $(LIBDIRS)

HEADERS= nnparams.h thermalign.h galign.h dinkelbach.h
SOURCES= nnparams.cpp thermalign.cpp galign.cpp dinkelbach.cpp talign_test.cpp

#SOURCES= inttest.cpp

OBJS   = $(SOURCES:.cpp=.o) 

thermalign: $(OBJS)
	$(CCC) $(LFLAGS) -o thermalign $(INCLUDE) $(LIBDIRS) $(OBJS) $(LIBS) $(EXTRALIBS)

#inttest: $(OBJS)
#	$(CCC) $(LFLAGS) -o inttest $(INCLUDE) $(LIBDIRS) $(OBJS) $(LIBS) $(EXTRALIBS)

clean: 
	rm *.o thermalign

%.o:%.cpp $(HEADERS)
	$(COMPILE.cc) -c $< -o $@
