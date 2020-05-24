
CC = ifort

SRC = $(wildcard *.f90)
TAR = $(SRC:.f90=.x)

.PHONY: all clean

all: $(TAR)
	echo set sluschipath=$(PWD) > $(HOME)/.sluschi.rc
	cp SLUSCHI StartJob autojob* getinfo.csh $(HOME)/bin



%.x: %.f90
	$(CC) -o $@ $<

clean:
	rm -f $(TAR)
	rm $(HOME)/bin/{SLUSCHI,StartJob,autojob_cont,getinfo.csh}

