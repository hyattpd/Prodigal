##############################################################################
#   PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
#   Copyright (C) 2007-2016 University of Tennessee / UT-Battelle
#
#   Code Author:  Doug Hyatt
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
##############################################################################

SHELL   = /bin/sh
CC      = gcc

CFLAGS  += -pedantic -Wall -O3 -DSUPPORT_GZIP_COMPRESSED
LFLAGS = -lm $(LDFLAGS) -lz

TARGET  = prodigal
INSTALLDIR  = /usr/local/bin

# Source file groups
CORE_SOURCES = bitmap.c dprog.c gene.c metagenomic.c node.c sequence.c training.c
API_SOURCES  = prodigal_api.c
CLI_SOURCE   = main.c
TEST_SOURCE  = test_api.c

HEADERS = $(shell echo *.h)

# CLI objects (no PRODIGAL_NO_MAIN)
CLI_OBJ = $(CLI_SOURCE:.c=.o)
CLI_CORE_OBJS = $(CORE_SOURCES:.c=.o)
CLI_API_OBJ = $(API_SOURCES:.c=.o)

# Library objects (with PRODIGAL_NO_MAIN to suppress stderr writes)
LIB_CORE_OBJS = $(CORE_SOURCES:.c=.lib.o)
LIB_API_OBJ = $(API_SOURCES:.c=.lib.o)
LIB_OBJS = $(LIB_CORE_OBJS) $(LIB_API_OBJ)

# Default: build CLI binary
all: $(TARGET)

# Objects depend on the compiler, its version and the flags, so changing any
# of them rebuilds everything instead of mixing stale objects.
BUILD_FLAGS := $(CC) $(shell $(CC) --version 2>/dev/null | head -1) $(CFLAGS)
.build-flags: FORCE
	@echo '$(BUILD_FLAGS)' | cmp -s - $@ || echo '$(BUILD_FLAGS)' >$@
FORCE:

# CLI objects: compiled WITHOUT PRODIGAL_NO_MAIN
%.o: %.c $(HEADERS) .build-flags
	$(CC) $(CFLAGS) -c -o $@ $<

# Library objects: compiled WITH PRODIGAL_NO_MAIN (suppresses stderr in core)
%.lib.o: %.c $(HEADERS) .build-flags
	$(CC) $(CFLAGS) -DPRODIGAL_NO_MAIN -c -o $@ $<

# Static library
libprodigal.a: $(LIB_OBJS)
	ar rcs $@ $^

# CLI binary: link CLI main.o + core objects (not using lib to keep stderr)
$(TARGET): $(CLI_OBJ) $(CLI_CORE_OBJS) $(CLI_API_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

# PIC objects for shared library
%.pic.o: %.c $(HEADERS) .build-flags
	$(CC) $(CFLAGS) -fPIC -DPRODIGAL_NO_MAIN -c -o $@ $<

LIB_PIC_OBJS = $(CORE_SOURCES:.c=.pic.o) $(API_SOURCES:.c=.pic.o)

# Shared library
libprodigal.so: $(LIB_PIC_OBJS)
	$(CC) -shared -o $@ $^ -lm

# Test runner (linked against static library)
test_api: $(TEST_SOURCE) libprodigal.a
	$(CC) $(CFLAGS) -DPRODIGAL_NO_MAIN -o $@ $< -L. -lprodigal $(LFLAGS)

test: test_api
	./test_api

# CLI output vs testdata/ground_truth; REGRESS_FLAGS passes --binary or
# --no-binary through (see docs/testing.md)
regress: $(TARGET)
	tests/cli_regress.sh $(REGRESS_FLAGS)

# CLI vs an upstream c1e2d36 build across a flag matrix and corpus; set ECOLI
# to an E. coli FASTA for the genome-scale cases (see docs/testing.md)
oracle: $(TARGET)
	CC="$(CC)" tests/oracle_diff.sh

install: $(TARGET)
	install -d -m 0755 $(INSTALLDIR)
	install -m 0755 $(TARGET) $(INSTALLDIR)

uninstall:
	-rm $(INSTALLDIR)/$(TARGET)

clean:
	-rm -f *.o *.lib.o *.pic.o

distclean: clean
	-rm -f $(TARGET) libprodigal.a libprodigal.so test_api .build-flags

.PHONY: all install uninstall clean distclean test regress oracle FORCE
