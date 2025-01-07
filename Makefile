CC         ?= gcc
CFLAGS     ?= -std=c11 -Wall -Wextra -O3 -MMD -MP
SRC         = $(filter-out $(wildcard test_*.c),$(wildcard *.c))
OBJ         = $(SRC:.c=.o)
HDR         = liknorm.h
TEST_SRC    = $(wildcard test_*.c)
TEST_OBJ    = $(TEST_SRC:.c=.o)
TEST_TARGET = $(basename $(TEST_OBJ))


ifeq ($(OS),Windows_NT)
	LIB     ?= liknorm.lib
	PREFIX  ?= C:\Program Files\Common Files
	COPYLIB  = Copy-Item
	COPYHDR  = Copy-Item
else
	LIB     ?= libliknorm.a
	PREFIX  ?= /usr/local
	COPYLIB  = install -m 0755
	COPYHDR  = install -m 0655
endif


all: $(LIB)

$(LIB): $(OBJ)
	ar rcs $@ $^

-include $(SRC:.c=.d)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_TARGET): %: %.o $(LIB)
	$(CC) $(CFLAGS) $< -L. $(LIB) -lm -o $@

check_integration: test_integration
	./test_integration

check_logprod: test_logprod
	./test_logprod

check_specific: test_specific
	./test_specific

check: check_integration check_logprod check_specific

install: $(LIB) $(HDR)
	@mkdir  -p $(PREFIX)\lib $(PREFIX)\include
	$(COPYLIB) $(LIB)        $(PREFIX)\lib
	$(COPYHDR) $(HDR)        $(PREFIX)\include

uninstall:
	rm -f $(PREFIX)/lib/$(LIB) $(HDR:%=$(PREFIX)/include/%)

.PHONY: all clean check
clean:
	rm -f $(OBJ) $(LIB) $(TEST_OBJ) $(TEST_TARGET) *.d
