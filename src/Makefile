# (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

COMPILE.CPP = g++
LINK.CPP = g++
AR = ar -rs

FLAGS = -c -I/opt/local/include

.SUFFICES:
.SUFFICES: .cpp .o

HEADERS = \
	times.h\
	constants.h

CLASSES = \
	time.o

OBJECTS = \
	time_conversion.o

MAIN = \
	main.o

#DEPS := $(OBJS:.o=.d)
#
#-include $(DEPS)
#
#%.o: %.cpp
#	$(COMPILE.CPP) $(FLAGS) -MMD -c -o $@ $<

.cpp.o: ${HEADERS}
	${COMPILE.CPP} ${FLAGS} $<

MY_LIB = libAether.a

LIB: ${CLASSES} ${OBJECTS}
	${AR} ${MY_LIB} ${CLASSES} ${OBJECTS}


Aether: ${MAIN} LIB
	${LINK.CPP} -o aether.exe ${MAIN} ${MY_LIB} -L../lib #-lmsis -L/opt/local/lib -lnetcdf-cxx4

clean:
	rm -f *~ core *.o *.exe *.a *.so *.d

