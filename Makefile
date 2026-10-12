PROJNAME := PAZ_Engine
CXXVER := 17
MINMACOSVER := 11.0

LIBNAME := $(shell echo $(PROJNAME) | sed 's/_//g' | tr '[:upper:]' '[:lower:]')
ifeq ($(OS), Windows_NT)
    OSPRETTY := Windows
    ifeq ($(MSYSTEM), UCRT64)
        CC := gcc
        CXX := g++
        LIBPATH := /ucrt64/lib
        INCLPATH := /ucrt64/include
    else ifeq ($(MSYSTEM), CLANG64)
        CC := clang
        CXX := clang++
        LIBPATH := /clang64/lib
        INCLPATH := /clang64/include
    else
        $(error Unsupported Windows environment)
    endif
else
    ifeq ($(shell uname -s), Darwin)
        OSPRETTY := macOS
        CC := clang
        CXX := clang++
    else
        OSPRETTY := Linux
        CC := gcc
        CXX := g++
    endif
    LIBPATH := /usr/local/lib
    INCLPATH := /usr/local/include
endif
OPTIM := 3
ZIPNAME := $(PROJNAME)-$(OSPRETTY)
ZIPCONTENTS := $(PROJNAME) lib$(LIBNAME).a
CFLAGS := -O$(OPTIM) -Wall -Wextra -Wno-missing-braces
ifeq ($(OSPRETTY), macOS)
    CFLAGS += -mmacosx-version-min=$(MINMACOSVER) -Wunguarded-availability
else
    ifeq ($(OSPRETTY), Windows)
        CFLAGS += -Wno-cast-function-type
    endif
endif
CXXFLAGS := -std=c++$(CXXVER) $(CFLAGS) -Wold-style-cast
ifeq ($(CC), clang)
    CXXFLAGS += -Wno-string-plus-int
endif
ARFLAGS := -rcs

ASSETS := $(wildcard assets/*)
CSRC := $(wildcard *.c) $(wildcard *.cpp)
ifeq ($(OSPRETTY), macOS)
    MACOSEXCL := $(patsubst %_macos.mm, %.cpp, $(wildcard *_macos.mm))
    CSRC := $(filter-out $(MACOSEXCL), $(CSRC))
endif
OBJCSRC := $(wildcard *.mm)
ifeq ($(OSPRETTY), macOS)
    ARMOBJ := assets_arm64.o $(patsubst %.c, %_arm64.o, $(patsubst %.cpp, %.c, \
        $(CSRC))) $(OBJCSRC:%.mm=%_arm64.o)
    INTOBJ := assets_x86_64.o $(patsubst %.c, %_x86_64.o, $(patsubst %.cpp, \
        %.c, $(CSRC))) $(OBJCSRC:%.mm=%_x86_64.o)
else
    OBJ := assets.o $(patsubst %.c, %.o, $(patsubst %.cpp, %.c, $(CSRC)))
endif

print-% : ; @echo $* = $($*)

.PHONY: test
default: test

ifeq ($(OSPRETTY), macOS)
lib$(LIBNAME).a: lib$(LIBNAME)_arm64.a lib$(LIBNAME)_x86_64.a
	lipo -create -output $@ $^

lib$(LIBNAME)_arm64.a: $(ARMOBJ)
	$(RM) $@
	ar $(ARFLAGS) $@ $^

lib$(LIBNAME)_x86_64.a: $(INTOBJ)
	$(RM) $@
	ar $(ARFLAGS) $@ $^
else
lib$(LIBNAME).a: $(OBJ)
	$(RM) $@
	ar $(ARFLAGS) $@ $^
endif

install: $(PROJNAME) lib$(LIBNAME).a
	@[[ -d $(INCLPATH) ]] || \
	    { echo "mkdir -p $(INCLPATH)"; mkdir -p $(INCLPATH); }
	@cmp -s $(PROJNAME) $(INCLPATH)/$(PROJNAME) && \
	    echo "Nothing to do for $(INCLPATH)/$(PROJNAME)" || \
	    { echo "cp $(PROJNAME) $(INCLPATH)/"; \
	    cp $(PROJNAME) $(INCLPATH)/; }
	@[[ -d $(LIBPATH) ]] || \
	    { echo "mkdir -p $(LIBPATH)"; mkdir -p $(LIBPATH); }
	@cmp -s lib$(LIBNAME).a $(LIBPATH)/lib$(LIBNAME).a && \
	    echo "Nothing to do for $(LIBPATH)/lib$(LIBNAME).a" || \
	    { echo "cp lib$(LIBNAME).a $(LIBPATH)/"; \
	    cp lib$(LIBNAME).a $(LIBPATH)/; }

test: lib$(LIBNAME).a
	$(MAKE) -C test
	test/test 2

analyze: $(OBJCSRC)
	$(foreach n, $(OBJCSRC), clang++ --analyze $(n) $(CXXFLAGS) && $(RM) \
	    $(n:%.mm=%.plist);)

%_arm64.o: %.cpp
	$(CXX) -arch arm64 -c -o $@ $< $(CXXFLAGS)

%_x86_64.o: %.cpp
	$(CXX) -arch x86_64 -c -o $@ $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS)

%_arm64.o: %.c
	$(CC) -arch arm64 -c -o $@ $< $(CFLAGS)

%_x86_64.o: %.c
	$(CC) -arch x86_64 -c -o $@ $< $(CFLAGS)

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

%_arm64.o: %.mm
	$(CC) -arch arm64 -c -o $@ $< $(CXXFLAGS)

%_x86_64.o: %.mm
	$(CC) -arch x86_64 -c -o $@ $< $(CXXFLAGS)

%.o: %.mm
	$(CC) -c -o $@ $< $(CXXFLAGS)

ifeq ($(OSPRETTY), macOS)
assets_arm64.o assets_x86_64.o: assets.pazarchive
	@printf ".section assets, \"rd\"\n.global _paz_binary_assets_pazarchive_start\n.global _paz_binary_assets_pazarchive_end\n_paz_binary_assets_pazarchive_start: .incbin \"assets.pazarchive\"\n_paz_binary_assets_pazarchive_end:" > include-assets.s
	as -arch arm64 -mmacosx-version-min=$(MINMACOSVER) include-assets.s -o assets_arm64.o
	as -arch x86_64 -mmacosx-version-min=$(MINMACOSVER) include-assets.s -o assets_x86_64.o
	$(RM) include-assets.s
assets_x86_64.o: assets_arm64.o
else ifeq ($(CC), gcc)
assets.o: assets.pazarchive
	ld -r -b binary -o assets.o assets.pazarchive
	objcopy --prefix-symbols=_paz assets.o
else
assets.o: assets.pazarchive
	@printf ".section assets, \"rd\"\n.global _paz_binary_assets_pazarchive_start\n.global _paz_binary_assets_pazarchive_end\n_paz_binary_assets_pazarchive_start: .incbin \"assets.pazarchive\"\n_paz_binary_assets_pazarchive_end:" > include-assets.s
	$(CC) -c include-assets.s -o assets.o
	$(RM) include-assets.s
endif

assets.pazarchive: $(ASSETS)
	paz-archive -c assets

clean:
	$(RM) *.o *.a *.pazarchive
	$(MAKE) -C test clean

zip: $(ZIPCONTENTS)
	zip -j $(ZIPNAME).zip $(ZIPCONTENTS)
