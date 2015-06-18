##############################################################################
#   PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
#   Copyright (C) 2007-2014 University of Tennessee / UT-Battelle
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

CFLAGS  += -pedantic -Wall -O3
LFLAGS = -lm $(LDFLAGS)

TARGET  = prodigal
SOURCES = $(shell echo *.c)
HEADERS = $(shell echo *.h)
OBJECTS = $(SOURCES:.c=.o)

EMPTY         =
PANDOC        = pandoc
PREFIX        = /usr/local
INSTALLDIR    = $(PREFIX)/bin
INSTALLMANDIR = $(PREFIX)/share/man/man1

all: $(TARGET) doc

doc: $(TARGET).1

$(TARGET).1:
ifeq ($(EMPTY), $(shell which $(PANDOC) 2>/dev/null))
	@echo "Failed to find $(PANDOC), no documentation"
else
	$(PANDOC) -s -V title=$(TARGET) -t man -f markdown -o $(TARGET).1 README.md
endif

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c -o $@ $<

install: $(TARGET)
	install -d -m 0756 $(INSTALLDIR)
	install -m 0755 $(TARGET) $(INSTALLDIR)
	install -d -m 756 $(INSTALLMANDIR)
	install $(TARGET).1 $(INSTALLMANDIR)
 
uninstall:
	-rm $(INSTALLDIR)/$(TARGET)
	-rm $(INSTALLMANDIR)/$(TARGET).1

clean:
	-rm -f $(OBJECTS)

distclean: clean
	-rm -f $(TARGET) $(TARGET).1

.PHONY: all install uninstall clean distclean
