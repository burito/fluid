
# Build rules
COMPANY = Daniel Burke
COPYRIGHT = 2013-2020
DESCRIPTION = Fluid Test
BINARY_NAME = fluid
OBJS = main.o version.o log.o global.o 3dmaths.o glerror.o vr.o \
	shader.o fps_movement.o fluid.o mesh.o mesh_gl.o stb_ds.o stb_image.o \
	fluidtest.o octtree.o spacemouse.o vr_helper.o
CFLAGS = -std=c11 -Wall -isystem deps -Ideps/dpb/src -Ideps/dpb/deps/hidapi/hidapi
VPATH = src build deps deps/dpb/src deps/dpb/deps deps/dpb/deps/hidapi/

WIN_LIBS = -lshell32 -luser32 -lgdi32 -lopengl32 -lwinmm -lws2_32 \
	-lxinput9_1_0
LIN_LIBS = -lm -lGL -lX11 -lGLU -lXi -ldl -rpath .
MAC_LIBS = deps/openvr/bin/osx32/libopenvr_api.dylib -framework OpenGL -framework CoreVideo -framework Cocoa -framework IOKit -rpath .

_WIN_OBJS = glew.o win32.o gfx_gl_win.o win32.res windows/hid.o $(OBJS)
_LIN_OBJS = glew.o linux_xlib.o gfx_gl_lin.o linux/hid.o $(OBJS)
_MAC_OBJS = osx.o gfx_gl_osx.o mac/hid.o $(OBJS)


include deps/dpb/src/Makefile

$(BINARY_NAME).exe: $(WIN_OBJS) openvr_api.dll

$(BINARY_NAME): $(LIN_OBJS) libopenvr_api.so

$(BINARY_NAME).bin: $(MAC_OBJS) libopenvr_api.dylib


openvr_api.dll: deps/openvr/bin/win64/openvr_api.dll
	cp $< $@

libopenvr_api.dylib: deps/openvr/bin/osx32/libopenvr_api.dylib
	cp $< $@

libopenvr_api.so: deps/openvr/bin/linux64/libopenvr_api.so
	cp $< $@


# this has to list everything inside the app bundle
$(MAC_CONTENTS)/_CodeSignature/CodeResources : \
	$(MAC_CONTENTS)/MacOS/$(BINARY_NAME) \
	$(MAC_CONTENTS)/Resources/AppIcon.icns \
	$(MAC_CONTENTS)/Frameworks/libopenvr_api.dylib \
	$(MAC_CONTENTS)/Info.plist
	codesign --force --deep --sign - $(BINARY_NAME).app

$(MAC_CONTENTS)/Frameworks/libopenvr_api.dylib: libopenvr_api.dylib
	@mkdir -p $(MAC_CONTENTS)/Frameworks
	cp $< $@

# copies the binary, and tells it where to find libraries
$(MAC_CONTENTS)/MacOS/$(BINARY_NAME): $(BINARY_NAME).bin
	@mkdir -p $(MAC_CONTENTS)/MacOS
	cp $< $@
	install_name_tool -change @loader_path/libopenvr_api.dylib @loader_path/../Frameworks/libopenvr_api.dylib $@
	install_name_tool -add_rpath "@loader_path/../Frameworks" $@
# end build the App Bundle

.PHONY:lint
lint:
	clang-tidy src/* -checks=-clang-analyzer-security.insecureAPI.DeprecatedOrUnsafeBufferHandling -- $(CFLAGS)
