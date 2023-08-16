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
LDLIBS := ../libpazengine.a $(LIBPATH)/libpazgraphics.a $(LIBPATH)/libpazio.a $(LIBPATH)/libpazmath.a
ifeq ($(OSPRETTY), macOS)
    LDLIBS += -framework MetalKit -framework Metal -framework Cocoa
else
    LDLIBS += $(LIBPATH)/libglfw3.a
    ifeq ($(OSPRETTY), Linux)
        LDLIBS += -lGL -lX11 -ldl
    else
        LDLIBS += -lgdi32 -lopengl32
    endif
endif

ASSETS := $(wildcard assets/*)

.PHONY: test.cpp
test: test.cpp build/assets.o

build/assets.o: build/assets.paz
ifeq ($(OSPRETTY), macOS)
	@printf "section .rodata\nglobal binary_assets_paz_start\nglobal binary_assets_paz_end\n_assets_start: incbin 'build/assets.paz'\n_assets_end:" > build/include-assets.s
	nasm -f $(ASMARCH) build/include-assets.s -o build/assets.o
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
