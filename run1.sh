########################################################################
########################################################################
###	   This file compile SI3D with some subroutines of GOTM

###	Check permissions of run.sh
###     if permission denied --> type 'sudo chmod 755 run.sh'
###     To compile --> './run.sh'

###     MAC

########################################################################
########################################################################

SHELL=$!/bin/sh

export FORTRAN_COMPILER=IFORT
export SI3DDIR=/mnt/c/Users/alcortes/Github/si3d
export GOTMDIR=/mnt/c/Users/alcortes/Github/gotm
export MODDIR=$GOTMDIR/modules
export INCDIR=$GOTMDIR/include
export BINDIR=$GOTMDIR/bin
export LIBDIR=$GOTMDIR/lib

#IF IFGOTM=false --> The modules 'util' and 'turbulence' are not compiled
export IFGOTM=true
#IF IFSI3D=false --> Only SI3D is compiled
export IFSI3D=true

if $IFGOTM
then
	cd $GOTMDIR/src/util
	make clean
	make
	cd $GOTMDIR/src/turbulence
	make clean
	make
fi

if $IFSI3D
then
	cd $SI3DDIR
	make all
	#rm *.o
	cd $MODDIR
	rm si3d*
else
	cd $SI3DDIR
	make si3d
	rm *.o
fi
