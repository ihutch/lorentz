include ACCIS.mk
##################################################################
% : %.f $(ACCISX) $(LIBDEPS) makefile;
	$(FORTRAN)  -o $* $(COMPILE-SWITCHES) $*.f $(LIBPATH) $(LIBRARIES)

%.o : %.f makefile ;
	$(FORTRAN) -c $(COMPILE-SWITCHES) $*.f


