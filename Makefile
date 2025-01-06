CC   ?= gcc

CFLAGS     ?= -std=c11 -Wall -Wextra -O3 -MMD -MP
PREFIX     ?= /usr/local
SRC         = $(filter-out $(wildcard test_*.c),$(wildcard *.c))
OBJ         = $(SRC:.c=.o)
HDR         = liknorm.h
LIB         = libliknorm.a
TEST_SRC    = $(wildcard test_*.c)
TEST_OBJ    = $(TEST_SRC:.c=.o)
TEST_TARGET = $(basename $(TEST_OBJ))


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
	@mkdir  -p      $(PREFIX)/lib $(PREFIX)/include
	install -m 0755 $(LIB)        $(PREFIX)/lib/
	install -m 0644 $(HDR)        $(PREFIX)/include/

uninstall:
	rm -f $(PREFIX)/lib/$(LIB) $(HDR:%=$(PREFIX)/include/%)

.PHONY: all clean check
clean:
	rm -f $(OBJ) $(LIB) $(TEST_OBJ) $(TEST_TARGET) *.d
