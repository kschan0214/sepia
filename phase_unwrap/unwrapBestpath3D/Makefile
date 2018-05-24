NAME = HusseinUnwrapper

CP = /usr/bin/g++
CC = /usr/bin/gcc
CFLAGS = -Wall
SOURCE = $(addsuffix .cpp, $(NAME))
ifeq "$(shell uname -m)" "x86_64"
PROCESSOR = x86_64
else
PROCESSOR = 32
endif	

Compile: $(PROCESSOR)

32: $(SOURCE)
	@echo "Carefull, the 32bit version may not be able to unwrapp any volume"
	$(CP) $(CFLAGS)					 $< -m32 -o $(NAME)

x86_64: $(SOURCE)
	$(CP) $(CFLAGS)					 $< -m64 -o $(NAME)

debug: $(SOURCE)
	$(CP) $(CFLAGS)					 $< -g -m64 -o $(NAME)_debug

.PHONY: clean
clean:
	rm HusseinUnwrapper

help:
	@echo "make         default compilation"
	@echo ".. debug     enable debugging symbols"
	@echo ".. 32        create 32bit version"
	@echo ".. x86_64    create 64bit version"