
CC = ifort

SRC = $(wildcard *.f90)
TAR = $(SRC:.f90=.x)

.PHONY: all clean

all: $(TAR)
	echo set sluschipath=$(PWD) > $(HOME)/.sluschi.rc
	cp SLUSCHI StartJob autojob* getinfo.csh $(HOME)/bin
	@echo "Checking if DIRECTORY is already in PATH..."
	@if grep -q 'export PATH=$$PATH:$(PWD)' ~/.bashrc; then \
		echo "DIRECTORY is already in PATH."; \
	else \
		echo "Adding DIRECTORY to PATH..."; \
		echo "export PATH=$$PATH:$(PWD)" | tee -a $(HOME)/.bashrc; \
		echo "DIRECTORY added to PATH. Please restart your shell."; \
	fi
	@if grep -q 'export PATH=$$PATH:$(PWD)' ~/.zshrc; then \
		echo "DIRECTORY is already in PATH."; \
	else \
		echo "Adding DIRECTORY to PATH..."; \
		echo "export PATH=$$PATH:$(PWD)" | tee -a $(HOME)/.zshrc; \
		echo "DIRECTORY added to PATH. Please restart your shell."; \
	fi

%.x: %.f90
	$(CC) -o $@ $<

clean:
	rm -f $(TAR)
	rm $(HOME)/bin/{SLUSCHI,StartJob,autojob_cont,getinfo.csh}

