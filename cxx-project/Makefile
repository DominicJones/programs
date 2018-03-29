cpp_compiler:=g++
cpp_flags:=-static

target_path:=../bin/
target:=project.exe

sources:=$(wildcard *.cpp)
objects:=$(sources:.cpp=.o)
obj_dep_rules:=$(sources:.cpp=.mk)

objects:=$(addprefix $(target_path), $(objects))
obj_dep_rules:=$(addprefix $(target_path), $(obj_dep_rules))
target:=$(addprefix $(target_path), $(target))

make_target_path:=mkdir -p $(target_path)

.PHONY: clean default

default: $(target)

$(target): $(objects)
	@$(make_target_path)
	$(cpp_compiler) $(cpp_flags) $^ -o $@

$(target_path)%.o: $(target_path)%.mk
	@$(make_target_path)
	$(cpp_compiler) $(cpp_flags) -c $(basename $*).cpp -o $@

$(target_path)%.mk: %.cpp
	@$(make_target_path)
	$(cpp_compiler) $(cpp_flags) -MM $< >> $@

clean:
	@rm -rf $(target_path)

ifneq ($(MAKECMDGOALS), clean)
-include $(obj_dep_rules)
endif
