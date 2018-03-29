ROOT:=.
include make.inc

all:
	$(MAKE) libs
	$(MAKE) cfdpack

libs:
	$(MAKE) -C $(LIB_ROOT)/ParMetis-3.1.1
	$(MAKE) -C $(LIB_ROOT)/pARMS_3.3
	$(MAKE) -C $(LIB_ROOT)/utils

cfdpack:
	$(MAKE) -C $(BIN_ROOT)/pgmsh
	$(MAKE) -C $(BIN_ROOT)/psolv

clean:
	$(MAKE) -C $(LIB_ROOT)/ParMetis-3.1.1 clean
	$(MAKE) -C $(LIB_ROOT)/pARMS_3.3 cleanobj
	$(MAKE) -C $(LIB_ROOT)/utils clean
	$(MAKE) -C $(BIN_ROOT)/pgmsh clean
	$(MAKE) -C $(BIN_ROOT)/psolv clean
	@find ./ -name \*~ -exec rm {} \;

distclean:
	$(MAKE) -C $(LIB_ROOT)/ParMetis-3.1.1 realclean
	$(MAKE) -C $(LIB_ROOT)/pARMS_3.3 cleanall
	$(MAKE) -C $(LIB_ROOT)/utils distclean
	$(MAKE) -C $(BIN_ROOT)/pgmsh distclean
	$(MAKE) -C $(BIN_ROOT)/psolv distclean
	@find ./ -name \*~ -exec rm {} \;
	@find cas/ -name \*.msh\* -exec rm {} \;
