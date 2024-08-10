CXX=g++
IDIR=.
CFLAGS=-I$(IDIR) -Wall -pedantic
CXXFLAGS=-std=c++20# -fsanitize=address
LIBS=-lm

DEPS= $(patsubst %,$(IDIR)/%,$(_DEPS))

ODIR=obj

_OBJ=solver.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cc $(DEPS)
	mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $(CXXFLAGS) $< -o $@ $(LIBS)

# Note: PHONY is important here. Without it, implicit rules will try to build the executable "all", since the prereqs are ".o" files.
.PHONY: all 

all: $(OBJ)

.PHONY: clean

clean:
	rm -r $(ODIR)

# Debugging Makefiles
# Try make -nd all
# make -n all
# make -np all
# man make has a lot more details.
