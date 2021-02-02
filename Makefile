all:
	ifort  -module $(MODDIR) -I/$(INCDIR) nspcg.f PlumeModels.f si3d_types.f90 si3d_boundaryconditions.f90 si3d_ecomod.f90 si3d_mixing.f90 si3d_procedures.f90 si3d.f90 -c -O2
	ifort  -o si3d *.o -module $(MODDIR) -I/$(INCDIR) -L/$(LIBDIR) -lturbulence_prod -lutil_prod -O2
si3d:	ifort  nspcg.f PlumeModels.f si3d_types.f90 si3d_ecomod.f90 si3d_boundaryconditions.f90 si3d_mixing.f90 si3d_procedures.f90 si3d.f90 -o si3d -O2
