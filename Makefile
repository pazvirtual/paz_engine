PROJNAME := PAZ_Engine
LIBNAME := $(shell echo $(PROJNAME) | sed 's/_//g' | tr '[:upper:]' '[:lower:]')
ifeq ($(OS), Windows_NT)
    LIBPATH := /mingw64/lib
    INCLPATH := /mingw64/include
    OSPRETTY := Windows
else
    ifeq ($(shell uname -s), Darwin)
        OSPRETTY := macOS
    else
        OSPRETTY := Linux
    endif
    LIBPATH := /usr/local/lib
    INCLPATH := /usr/local/include
endif
CXXVER := 17
OPTIM := 3
ZIPNAME := $(PROJNAME)-$(OSPRETTY)
CFLAGS := -O$(OPTIM) -Wall -Wextra -Wno-missing-braces
MINMACOSVER := 10.11
ifeq ($(OSPRETTY), macOS)
    CFLAGS += -mmacosx-version-min=$(MINMACOSVER) -Wunguarded-availability
else
    ifeq ($(OSPRETTY), Windows)
        CFLAGS += -Wno-cast-function-type
    endif
endif
CXXFLAGS := -std=c++$(CXXVER) $(CFLAGS) -Wold-style-cast
ifeq ($(OSPRETTY), Windows)
    CXXFLAGS += -Wno-deprecated-copy
endif
ARFLAGS := -rcs

ASSETS := $(wildcard assets/*)
CSRC := $(wildcard *.c) $(wildcard *.cpp)
ifeq ($(OSPRETTY), macOS)
    MACOSEXCL := $(patsubst %_macos.mm, %.cpp, $(wildcard *_macos.mm))
    CSRC := $(filter-out $(MACOSEXCL), $(CSRC))
endif
OBJCSRC := $(wildcard *.mm)
OBJ := assets.o $(patsubst %.c, %.o, $(patsubst %.cpp, %.c, $(CSRC)))
ifeq ($(OSPRETTY), macOS)
    OBJ += $(OBJCSRC:%.mm=%.o)
endif

REINSTALLHEADER := $(shell cmp -s $(PROJNAME) $(INCLPATH)/$(PROJNAME); echo $$?)
REINSTALLLIB := $(shell cmp -s lib$(LIBNAME).a $(LIBPATH)/lib$(LIBNAME).a; echo $$?)

print-% : ; @echo $* = $($*)

.PHONY: test
default: test

lib$(LIBNAME).a: $(OBJ)
	$(RM) lib$(LIBNAME).a
	ar $(ARFLAGS) lib$(LIBNAME).a $^

ifneq ($(REINSTALLHEADER), 0)
    ifneq ($(REINSTALLLIB), 0)
install: $(PROJNAME) lib$(LIBNAME).a
	cp $(PROJNAME) $(INCLPATH)/
	cp lib$(LIBNAME).a $(LIBPATH)/
    else
install: $(PROJNAME) lib$(LIBNAME).a
	cp $(PROJNAME) $(INCLPATH)/
    endif
else
    ifneq ($(REINSTALLLIB), 0)
install: $(PROJNAME) lib$(LIBNAME).a
	cp lib$(LIBNAME).a $(LIBPATH)/
    else
install: $(PROJNAME) lib$(LIBNAME).a
    endif
endif

test: lib$(LIBNAME).a
	$(MAKE) -C test
	test/test

analyze: $(OBJCSRC)
	$(foreach n, $(OBJCSRC), clang++ --analyze $(n) $(CXXFLAGS) && $(RM) $(n:%.mm=%.plist);)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

%.o: %.mm
	$(CC) -c -o $@ $< $(CXXFLAGS)

assets.o: assets.paz
ifeq ($(OSPRETTY), macOS)
	@printf ".section assets, \"rd\"\n.global _paz_binary_assets_paz_start\n.global _paz_binary_assets_paz_end\n_paz_binary_assets_paz_start: .incbin \"assets.paz\"\n_paz_binary_assets_paz_end:" > include-assets.s
	as -mmacos-version-min=$(MINMACOSVER) include-assets.s -o assets.o
	$(RM) include-assets.s
else
	ld -r -b binary -o assets.o assets.paz
	objcopy --prefix-symbols=_paz assets.o
endif

assets.paz: $(ASSETS)
	paz-archive assets

clean:
	$(RM) $(OBJ) lib$(LIBNAME).a assets.paz
	$(MAKE) -C test clean

zip: $(PROJNAME) lib$(LIBNAME).a
	zip -j $(ZIPNAME).zip $(PROJNAME) lib$(LIBNAME).a
