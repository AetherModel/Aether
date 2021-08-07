# (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

include Makefile.OS

.SUFFICES:
.SUFFICES: .cpp .o

CLASSES = \
	chemistry.o\
	electrodynamics.o\
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
	init_mag_grid.o\
	output.o\
	read_collision_file.o\
	read_f107_file.o\
	read_input_file.o\
	read_netcdf_electrodynamics.o\
	read_omni_file.o\
	report.o\
	solver_chemistry.o\
	solver_conduction.o\
	solver_interpolation.o\
	time_conversion.o\
	transform.o

MAIN = \
	main.o

TEST = \
	test.o

.cpp.o: ${HEADERS}
	${COMPILE.CPP} ${FLAGS} -l${NETCDF} -l${ARMADILLO} $<

MY_LIB = libAether.a

LIB: ${CLASSES} ${OBJECTS}
	${AR} ${MY_LIB} ${CLASSES} ${OBJECTS}


Aether: ${MAIN} LIB
	${LINK.CPP} -o aether.exe ${MAIN} ${MY_LIB} -L../lib -L${LIB_DIR} -l${NETCDF} -l${ARMADILLO}

test: ${TEST} LIB
	${LINK.CPP} -o test.exe ${TEST} ${MY_LIB} -L../lib -L${LIB_DIR} -l${NETCDF} -l${ARMADILLO}

clean:
	rm -f *~ core *.o *.exe *.a *.so *.d



TESTMGRID = \
	test_mgrid.o

test_mgrid: ${TESTMGRID} LIB
	${LINK.CPP} -o test_mgrid.exe ${TESTMGRID} ${MY_LIB} -L../lib -L${LIB_DIR} -l${NETCDF} -l${ARMADILLO}



