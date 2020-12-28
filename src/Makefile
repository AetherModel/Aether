# (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

COMPILE.CPP = g++
LINK.CPP = g++
AR = ar -rs

FLAGS = -O3 -ffast-math -c -I/opt/local/include
# FLAGS = -c -I/opt/local/include

.SUFFICES:
.SUFFICES: .cpp .o

HEADERS = \
	times.h\
	constants.h\
	inputs.h\
	file_input.h

CLASSES = \
	time.o\
	inputs.o\
	euv.o\
	indices.o\
	planets.o\
	grid.o\
	neutrals.o\
	ions.o\
	chemistry.o

OBJECTS = \
	time_conversion.o\
	transform.o\
	report.o\
	file_input.o\
	read_f107_file.o\
	init_geo_grid.o\
	fill_grid.o\
	calc_neutral_derived.o\
	calc_euv.o\
	solver_conduction.o\
	solver_chemistry.o\
	advance.o\
	add_sources.o\
	output.o\
	bfield.o\
	dipole.o\
	calc_chemistry.o\
	calc_chemical_sources.o

MAIN = \
	main.o

TEST = \
	test.o

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
	${LINK.CPP} -o aether.exe ${MAIN} ${MY_LIB} -L../lib -L/opt/local/lib -lnetcdf-cxx4 #-lmsis

test: ${TEST} LIB
	${LINK.CPP} -o test.exe ${TEST} ${MY_LIB} -L../lib

clean:
	rm -f *~ core *.o *.exe *.a *.so *.d

