ifeq ($(OS), Windows_NT)
    LIBPATH := /mingw64/lib
    OSPRETTY := Windows
else
    ifeq ($(shell uname -s), Darwin)
        OSPRETTY := macOS
    else
        OSPRETTY := Linux
    endif
    LIBPATH := /usr/local/lib
endif
CXXVER := 14
CXXFLAGS := -std=c++$(CXXVER) -O3 -Wall -Wextra -Wno-missing-braces -Wold-style-cast
MINMACOSVER := 10.11
ifeq ($(OSPRETTY), macOS)
    CXXFLAGS += -mmacosx-version-min=$(MINMACOSVER) -Wunguarded-availability
else
    ifeq ($(OSPRETTY), Windows)
        CXXFLAGS += -Wno-deprecated-copy -static -mwindows
    endif
    CXXFLAGS += -pthread
endif
CXXFLAGS += -I..
LDLIBS := ../libpazengine.a $(LIBPATH)/libpazgraphics.a $(LIBPATH)/libpazaudio.a $(LIBPATH)/libpazio.a $(LIBPATH)/libpazmath.a
ifeq ($(OSPRETTY), macOS)
    LDLIBS += -framework MetalKit -framework Metal -framework Cocoa -framework IOKit -framework CoreAudio -framework AudioUnit
else
    LDLIBS += $(LIBPATH)/libglfw3.a $(LIBPATH)/libportaudio.a
    ifeq ($(OSPRETTY), Linux)
        LDLIBS += -lGL -lX11 -ldl -lasound
    else
        LDLIBS += -ld3d11 -ldxgi -ld3dcompiler -ldxguid -lwinmm -lole32 -lsetupapi
    endif
endif

ASSETS := $(wildcard assets/*)
SRC := $(wildcard *.cpp)
OBJ := build/assets.o $(SRC:%.cpp=build/%.o)

.PHONY: test
test: build $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDLIBS)

build/assets.o: build/assets.paz
ifeq ($(OSPRETTY), macOS)
	@printf ".section assets, "rd"\n.global _binary_assets_paz_start\n.global _binary_assets_paz_end\n_binary_assets_paz_start: .incbin \"build/assets.paz\"\n_binary_assets_paz_end:" > build/include-assets.s
	as -mmacosx-version-min=$(MINMACOSVER) build/include-assets.s -o build/assets.o
else
	cd build && ld -r -b binary -o assets.o assets.paz
endif

build/assets.paz: $(ASSETS)
	paz-archive assets build/assets.paz

build/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

build:
	mkdir -p $@

clean:
	$(RM) test build/*

print-% : ; @echo $* = $($*)
