# (c) 2020, the Aether Development Team (see doc/dev_team.md for members)
# Full license can be found in License.md

RUNDIR = run

default : Aether

Aether:
	cd src; make Aether

test_mgrid:
	cd src; make test_mgrid

clean:
	cd src; make clean

# We need to fix some path stuff in this. We need an absolute path to make
# this better:

rundir:
	mkdir -p ${RUNDIR}/UA
	cd ${RUNDIR} ; mkdir UA/output ; ln -s UA/outputs . ; cd -
	cd ${RUNDIR}/UA ; mkdir restartOut ; cd -
	cd ${RUNDIR}/UA ; ln -s restartOut restartIn ; cd -
	mkdir -p ${RUNDIR}/UA/inputs ; cp inputs/* ${RUNDIR}/UA/inputs
	cd ${RUNDIR} ; ln -s ../src/aether.exe . ; cd -
	cd ${RUNDIR} ; cp UA/inputs/aether.in . ; cd -
