CXXVER := 17
MINMACOSVER := 10.12

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
CXXFLAGS := -std=c++$(CXXVER) -O3 -Wall -Wextra -Wno-missing-braces -Wold-style-cast
ifeq ($(OSPRETTY), macOS)
    CXXFLAGS += -mmacosx-version-min=$(MINMACOSVER) -Wunguarded-availability -Wno-string-plus-int
else
    ifeq ($(OSPRETTY), Windows)
        CXXFLAGS += -Wno-deprecated-copy -static -mwindows
    else
        CXXFLAGS += -pthread
    endif
endif
CXXFLAGS += -I..
LDLIBS := ../libpazengine.a $(LIBPATH)/libpazgraphics.a $(LIBPATH)/libpazaudio.a $(LIBPATH)/libpazio.a $(LIBPATH)/libpazmath.a
ifeq ($(OSPRETTY), macOS)
    LDLIBS += -framework MetalKit -framework Metal -framework Cocoa -framework IOKit -framework CoreAudio -framework AudioUnit
else
    ifeq ($(OSPRETTY), Linux)
        # Try a few locations for GLFW.
        ifneq (, $(wildcard $(LIBPATH)/libglfw3.a))
            LDLIBS += $(LIBPATH)/libglfw3.a
        else ifneq (, $(wildcard /usr/lib/libglfw3.a))
            LDLIBS += /usr/lib/libglfw3.a
        else ifneq (, $(wildcard /usr/lib64/libglfw3.a))
            LDLIBS += /usr/lib64/libglfw3.a
        else ifneq (, $(wildcard /usr/lib/x86_64-linux-gnu/libglfw3.a))
            LDLIBS += /usr/lib/x86_64-linux-gnu/libglfw3.a
        else
            $(error Could not find "libglfw3.a".)
        endif
        # Try a few locations for PortAudio.
        ifneq (, $(wildcard $(LIBPATH)/libportaudio.a))
            LDLIBS += $(LIBPATH)/libportaudio.a
        else ifneq (, $(wildcard /usr/lib/libportaudio.a))
            LDLIBS += /usr/lib/libportaudio.a
        else ifneq (, $(wildcard /usr/lib/x86_64-linux-gnu/libportaudio.a))
            LDLIBS += /usr/lib/x86_64-linux-gnu/libportaudio.a
        else
            $(error Could not find "libportaudio.a".)
        endif
        LDLIBS += -lGL -lX11 -ldl -lasound -ljack
    else
        LDLIBS += $(LIBPATH)/libportaudio.a -ld3d11 -ldxgi -ld3dcompiler -ldxguid -lwinmm -lole32 -lsetupapi -Wl,-Bstatic -lstdc++ -lpthread -Wl,-Bdynamic
        LDFLAGS += -static-libgcc -static-libstdc++
    endif
endif

ASSETS := $(wildcard assets/*)
SRC := $(wildcard *.cpp)
ifeq ($(OSPRETTY), macOS)
    ARMOBJ := build/assets_arm64.o $(SRC:%.cpp=build/%_arm64.o)
    INTOBJ := build/assets_x86_64.o $(SRC:%.cpp=build/%_x86_64.o)
else
    OBJ := build/assets.o $(SRC:%.cpp=build/%.o)
endif

print-% : ; @echo $* = $($*)

ifeq ($(OSPRETTY), macOS)
.PHONY: build/test_arm64 build/test_x86_64
test: build/test_arm64 build/test_x86_64
	lipo -create -output $@ $^

build/test_arm64: build $(ARMOBJ)
	$(CXX) -arch arm64 $(CXXFLAGS) -o $@ $(ARMOBJ) $(LDFLAGS) $(LDLIBS)

build/test_x86_64: build $(INTOBJ)
	$(CXX) -arch x86_64 $(CXXFLAGS) -o $@ $(INTOBJ) $(LDFLAGS) $(LDLIBS)
else
.PHONY: test
test: build $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ) $(LDFLAGS) $(LDLIBS)
endif

ifeq ($(OSPRETTY), macOS)
build/assets_arm64.o build/assets_x86_64.o: build/assets.pazarchive
	@printf ".section assets, \"rd\"\n.global _binary_assets_pazarchive_start\n.global _binary_assets_pazarchive_end\n_binary_assets_pazarchive_start: .incbin \"build/assets.pazarchive\"\n_binary_assets_pazarchive_end:" > build/include-assets.s
	as -arch arm64 -mmacosx-version-min=$(MINMACOSVER) build/include-assets.s -o build/assets_arm64.o
	as -arch x86_64 -mmacosx-version-min=$(MINMACOSVER) build/include-assets.s -o build/assets_x86_64.o
build/assets_x86_64.o: build/assets_arm64.o
else
build/assets.o: build/assets.pazarchive
	cd build && ld -r -b binary -o assets.o assets.pazarchive
endif

build/assets.pazarchive: $(ASSETS)
	paz-archive -c assets build/assets.pazarchive

build/%_arm64.o: %.cpp
	$(CXX) -arch arm64 $(CXXFLAGS) -c -o $@ $<

build/%_x86_64.o: %.cpp
	$(CXX) -arch x86_64 $(CXXFLAGS) -c -o $@ $<

build/%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

build:
	mkdir -p $@

clean:
	$(RM) -r test build
