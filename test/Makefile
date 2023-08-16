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
CXXFLAGS := -std=c++$(CXXVER) -Ofast -Wall -Wextra -Wno-missing-braces -Wold-style-cast
ifeq ($(OSPRETTY), macOS)
    CXXFLAGS += -mmacosx-version-min=10.10
else
    ifeq ($(OSPRETTY), Windows)
        CXXFLAGS += -Wno-deprecated-copy -static -mwindows
    endif
    CXXFLAGS += -pthread
endif
CXXFLAGS += -I..
LDLIBS := ../libpazengine.a $(LIBPATH)/libpazgraphics.a $(LIBPATH)/libpazaudio.a $(LIBPATH)/libpazio.a $(LIBPATH)/libpazmath.a
ifeq ($(OSPRETTY), macOS)
    LDLIBS += -framework MetalKit -framework Metal -framework Cocoa -framework IOKit
else
    LDLIBS += $(LIBPATH)/libglfw3.a $(LIBPATH)/libportaudio.a
    ifeq ($(OSPRETTY), Linux)
        LDLIBS += -lGL -lX11 -ldl -lasound
    else
        LDLIBS += -lgdi32 -lopengl32 -lwinmm -lole32 -lsetupapi
    endif
endif

ASSETS := $(wildcard assets/*)

.PHONY: test.cpp
test: test.cpp build/assets.o

build/assets.o: build/assets.paz
ifeq ($(OSPRETTY), macOS)
	@printf "section .rodata\nglobal _binary_assets_paz_start\nglobal _binary_assets_paz_end\n_binary_assets_paz_start: incbin 'build/assets.paz'\n_binary_assets_paz_end:" > build/include-assets.s
	nasm -f macho64 build/include-assets.s -o build/assets.o
else
	cd build && ld -r -b binary -o assets.o assets.paz
endif

build/assets.paz: $(ASSETS) build
	paz-archive assets build/assets.paz

build:
	mkdir -p $@

clean:
	$(RM) test build/*

print-% : ; @echo $* = $($*)
