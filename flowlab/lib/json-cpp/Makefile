include ../../etc/make/Makefile.variables

target:=jsoncpp

sources:=
sources+=json_reader.cpp
sources+=json_writer.cpp
sources+=json_value.cpp

objects:=$(addsuffix .o,$(basename $(sources)))

include ../../etc/make/Makefile.staticlib

include ../../etc/make/Makefile.compile
