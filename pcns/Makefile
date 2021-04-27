include etc/make/Makefile.variables

src_path:=src/
lib_path:=lib/


LIB:=
LIB+=json-cpp
LIB+=json-f90
ifeq ($(strip $(mpi)),yes)
  LIB+=ParMetis-3.1.1
  LIB+=pARMS_3.3
endif

LIB:=$(addprefix $(lib_path), $(LIB))
LIB_MAKE:=$(foreach item, $(LIB), $(MAKE) -C $(item);)
LIB_CLEAN:=$(foreach item, $(LIB), $(MAKE) -C $(item) clean;)

.PHONY: all
all: libs $(target_name)

clean: libs-clean $(target_name)-clean rm-junk

libs:
	$(LIB_MAKE)

$(target_name):
	$(MAKE) -C $(src_path)

$(target_name)-clean:
	$(MAKE) -C $(src_path) clean

libs-clean:
	$(LIB_CLEAN)

rm-junk:
	@find ./ -name \*~ -exec rm {} \;
