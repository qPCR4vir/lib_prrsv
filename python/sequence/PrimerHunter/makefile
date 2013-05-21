CFLAGS      = -fno-weak -c
LFLAGS      = -lm
MELTING     = MeltingFracProg
#CFLAGS     = -g -gall -ms -xtaso -stdnew -global_array_new -std ms -Zp1 -misalign -assume noaligned_objects -nomember_alignment -protect_headers all
#LFLAGS     = -g -gall -ms -xtaso -stdnew -global_array_new -std ms -Zp1 -misalign -assume noaligned_objects -nomember_alignment -protect_headers all
LIBS       = $(MELTING)/nnparams.o $(MELTING)/thermalign.o $(MELTING)/galign.o $(MELTING)/dinkelbach.o
EXTRALIBS  =

#PURIFY        = purify -g++ -best-effort
#CCC           = /zpr/gcc/2.8/bin/g++
#CCC           = cxx -O5
#CCC           = g++ -g -Wno-deprecated
CCC           = g++ -O3 -Wno-deprecated

INCLUDE       = -I./$(MELTING)
LIBDIRS       = 

COMPILE.cc = $(CCC) $(CFLAGS) $(INCLUDE) $(LIBDIRS)

SOURCES= primerhunter.cpp dpal.cpp

#SOURCES= inttest.cpp

OBJS   = $(SOURCES:.cpp=.o) 

all: $(OBJS)
	cd MeltingFracProg; make; cd ..
	$(CCC) $(LFLAGS) -o primerhunter $(LIBDIRS) $(OBJS) $(LIBS) $(EXTRALIBS)

clean: 
	rm *.o primerhunter; cd MeltingFracProg; make clean

%.o:%.cpp $(COMPILE.cc)  $< -o $@ 
	
