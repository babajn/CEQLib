# Keep the existing grid comparison targets available from test/chimes.
path_variables = GSL_DIR EQCLIB_DIR CHIMES_DIR CHIMES_DATA CHIMES_LOCAL PEGASE_DIR OUT

.PHONY: all eqclib run equilibrium timing courant cooling plot clean
all eqclib run equilibrium timing courant cooling plot clean:
	$(MAKE) -C shielded $@ $(foreach var,$(path_variables),$(if $($(var)),$(var)="$(abspath $($(var)))"))
