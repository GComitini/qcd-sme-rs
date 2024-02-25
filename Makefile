CARGO = cargo
CARGO_DOC = cargo-doc
CC = gcc
CYTHON = cython
PYTHON_CFLAGS = `pkg-config --cflags python3`
TARGETS = debug release

override BINARIES =
override CARGO_CLIPPY_CMDS =
override CARGO_TEST_CMDS =
override CDYLIBS =
override CYTHON_PXD_HEADERS =
override CYTHON_DYLIB_SRCS =
override CYTHON_DYLIBS =

override DEBUG_DIR = target/debug
override RELEASE_DIR = target/release
override CYTHON_PYX_SRC = python/pyqcd_sme_so.pyx
override RUST_FILES = $(shell find src/ -type f -name '*.rs')

override OS = $(shell uname -s)

ifeq ($(OS),Linux)
	override SHLIB_EXT = so
else
	override SHLIB_EXT = dylib
endif

ifneq (,$(filter $(TARGETS),debug))
	override CARGO_CLIPPY_CMDS += cargo-clippy-debug
	override CARGO_TEST_CMDS += cargo-test-debug
	override CDYLIBS += $(DEBUG_DIR)/libqcd_sme.$(SHLIB_EXT)
	override CYTHON_PXD_HEADERS += $(DEBUG_DIR)/include/qcd_sme.pxd
	override CYTHON_DYLIB_SRCS += $(DEBUG_DIR)/python/pyqcd_sme_so.c
	override CYTHON_DYLIBS += $(DEBUG_DIR)/python/pyqcd_sme_so.$(SHLIB_EXT)
endif
ifneq (,$(filter $(TARGETS),release))
	override CARGO_CLIPPY_CMDS += cargo-clippy-release
	override CARGO_TEST_CMDS += cargo-test-release
	override CDYLIBS += $(RELEASE_DIR)/libqcd_sme.$(SHLIB_EXT)
	override CYTHON_PXD_HEADERS += $(RELEASE_DIR)/include/qcd_sme.pxd
	override CYTHON_DYLIB_SRCS += $(RELEASE_DIR)/python/pyqcd_sme_so.c
	override CYTHON_DYLIBS += $(RELEASE_DIR)/python/pyqcd_sme_so.$(SHLIB_EXT)
endif


.PHONY: all cargo-clippy-debug cargo-clippy-release cargo-doc cargo-test-debug cargo-test-release check clean clippy docs python
all: $(BINARIES) $(CDYLIBS)

cargo-clippy-debug:
	$(CARGO) clippy

cargo-clippy-release:
	$(CARGO) clippy --release

cargo-doc: $(RUST_FILES)
	$(CARGO_DOC)

cargo-test-debug:
	$(CARGO) test

cargo-test-release:
	$(CARGO) test --release

check: $(CARGO_TEST_CMDS)

clean:
	rm -rf target

clippy: $(CARGO_CLIPPY_CMDS)

docs: cargo-doc

python: $(CYTHON_DYLIBS)

$(DEBUG_DIR)/include/qcd_sme.pxd $(DEBUG_DIR)/libqcd_sme.$(SHLIB_EXT) &: $(RUST_FILES)
	$(CARGO) build

$(DEBUG_DIR)/python/pyqcd_sme_so.c: $(CYTHON_PYX_SRC) $(CYTHON_PXD_HEADERS)
	mkdir -p $(DEBUG_DIR)/python
	$(CYTHON) -3 -I$(DEBUG_DIR)/include $(CYTHON_PYX_SRC) -o $(DEBUG_DIR)/python/pyqcd_sme_so.c

$(DEBUG_DIR)/python/pyqcd_sme_so.so: $(CYTHON_PXD_HEADERS) $(CYTHON_DYLIB_SRCS) $(DEBUG_DIR)/libqcd_sme.so
	$(CC) -O3 -Wall -shared -fPIC $(PYTHON_CFLAGS) -I$(DEBUG_DIR)/include -Wl,--no-as-needed -L./$(DEBUG_DIR) -lqcd_sme \
		-Wl,-rpath,$$(pwd)/$(DEBUG_DIR) -o $(DEBUG_DIR)/python/pyqcd_sme_so.so $(DEBUG_DIR)/python/pyqcd_sme_so.c

$(RELEASE_DIR)/include/qcd_sme.pxd $(RELEASE_DIR)/libqcd_sme.$(SHLIB_EXT) &: $(RUST_FILES)
	$(CARGO) build --release

$(RELEASE_DIR)/python/pyqcd_sme_so.c: $(CYTHON_PYX_SRC) $(CYTHON_PXD_HEADERS)
	mkdir -p $(RELEASE_DIR)/python
	$(CYTHON) -3 -I$(RELEASE_DIR)/include $(CYTHON_PYX_SRC) -o $(RELEASE_DIR)/python/pyqcd_sme_so.c

$(RELEASE_DIR)/python/pyqcd_sme_so.so: $(CYTHON_PXD_HEADERS) $(CYTHON_DYLIB_SRCS) $(RELEASE_DIR)/libqcd_sme.so
	$(CC) -O3 -Wall -shared -fPIC $(PYTHON_CFLAGS) -I$(RELEASE_DIR)/include -Wl,--no-as-needed -L./$(RELEASE_DIR) -lqcd_sme \
		-Wl,-rpath,$$(pwd)/$(RELEASE_DIR) -o $(RELEASE_DIR)/python/pyqcd_sme_so.so $(RELEASE_DIR)/python/pyqcd_sme_so.c
