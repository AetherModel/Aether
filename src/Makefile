# (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

include Makefile.OS

.SUFFICES:
.SUFFICES: .cpp .o

CLASSES = \
	chemistry.o\
	euv.o\
	grid.o\
	indices.o\
	inputs.o\
	ions.o\
	neutrals.o\
	planets.o\
	time.o

OBJECTS = \
	add_sources.o\
	advance.o\
	bfield.o\
	calc_chemistry.o\
	calc_chemical_sources.o\
	calc_euv.o\
	calc_neutral_derived.o\
	dipole.o\
	file_input.o\
	fill_grid.o\
	init_geo_grid.o\
	output_binary.o\
	output_netcdf.o\
	read_f107_file.o\
	read_input_file.o\
	read_omni_file.o\
	report.o\
	solver_chemistry.o\
	solver_conduction.o\
	time_conversion.o\
	transform.o

MAIN = \
	main.o

TEST = \
	test.o

.cpp.o: ${HEADERS}
	${COMPILE.CPP} ${FLAGS} $<

MY_LIB = libAether.a

LIB: ${CLASSES} ${OBJECTS}
	${AR} ${MY_LIB} ${CLASSES} ${OBJECTS}

Aether: ${MAIN} LIB
	${LINK.CPP} -o aether.exe ${MAIN} ${MY_LIB} -L../lib -L${LIB_DIR} -l${NETCDF}

test: ${TEST} LIB
	${LINK.CPP} -o test.exe ${TEST} ${MY_LIB} -L../lib -L${LIB_DIR} -l${NETCDF}

clean:
	rm -f *~ core *.o *.exe *.a *.so *.d

