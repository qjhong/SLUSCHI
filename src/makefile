FC = gfortran
FFLAGS ?=
LDFLAGS ?=
LDLIBS ?=

SRC := $(wildcard *.f90)
TAR := $(SRC:.f90=.x)
SCRIPTS := SLUSCHI StartJob $(wildcard autojob*) getinfo.csh

BINDIR ?= $(HOME)/bin
RCFILE ?= $(HOME)/.sluschi.rc

.PHONY: all install clean

all: $(TAR)

%.x: %.f90
	$(FC) $(FFLAGS) $(LDFLAGS) -o "$@" "$<" $(LDLIBS)

install: all
	mkdir -p "$(BINDIR)"
	cp $(SCRIPTS) "$(BINDIR)/"
	printf 'set sluschipath=%s\n' "$(CURDIR)" > "$(RCFILE)"
	@printf '\nFor Bash or Zsh, add this line once to your shell startup file:\n'
	@printf 'export PATH="$$PATH:%s:%s"\n' "$(CURDIR)" "$(BINDIR)"

clean:
	rm -f $(TAR)
