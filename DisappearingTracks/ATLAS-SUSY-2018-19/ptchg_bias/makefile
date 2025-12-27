include ../../make_opts

all: ptchg_bias

clean:
	$(RM) *.o $(BIASLIBDIR)$(BIASLIBRARY)

#
# Compilation of the module ptchg_bias
#
ptchg_bias.o: ptchg_bias.f ../bias.inc
	$(FC) $(FFLAGS) $(LDFLAGS) -c -o ptchg_bias.o ptchg_bias.f

ptchg_bias: ptchg_bias.o
	$(call CREATELIB, $(BIASLIBDIR)$(BIASLIBRARY), $^)

#
# List of the requirements for this module.
# 'VALID' is the keyword that *must* be returned if everything is in order.
#
requirements:	
	@echo "VALID"

