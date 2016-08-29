#
# TODO:
#

CC := g++ # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity

SRCDIR := src
BUILDDIR := build
TARGET:= bin/BDsim
TARGET-TD := bin/BDsim-td
SRCEXT := cpp

FOLDER := constant-rate
SOURCES := $(shell find $(SRCDIR)/$(FOLDER) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/$(FOLDER)/%,$(BUILDDIR)/$(FOLDER)/%,$(SOURCES:.$(SRCEXT)=.o))

FOLDER-TD := time-dependent
SOURCES-TD := $(shell find $(SRCDIR)/$(FOLDER-TD) -type f -name *$.$(SRCEXT))
OBJECTS-TD := $(patsubst $(SRCDIR)/$(FOLDER-TD)/%,$(BUILDDIR)/$(FOLDER-TD)/%,$(SOURCES-TD:.$(SRCEXT)=.o))

CFLAGS := -g -Wall -O3
LIB := -lgsl
LPATH := /usr/local/lib
INC := -I /usr/local/include

all: $(TARGET) $(TARGET-TD)

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) -L $(LPATH) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(TARGET-TD): $(OBJECTS-TD)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET-TD) -L $(LPATH) $(LIB)"; $(CC) $^ -o $(TARGET-TD) $(LIB)


$(BUILDDIR)/$(FOLDER)/%.o: $(SRCDIR)/$(FOLDER)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)/$(FOLDER)
	@echo " $(CC) $(CFLAGS) -I include/constant-rate $(INC) -L $(LPATH) -c -o $@ $<"; $(CC) $(CFLAGS) -I include/constant-rate $(INC) -c -o $@ $<

$(BUILDDIR)/$(FOLDER-TD)/%.o: $(SRCDIR)/$(FOLDER-TD)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)/$(FOLDER-TD)
	@echo " $(CC) $(CFLAGS) -I include/time-dependent $(INC) -L $(LPATH) -c -o $@ $<"; $(CC) $(CFLAGS) -I include/time-dependent $(INC) -c -o $@ $<

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(TARGET) $(TARGET-TD)"; $(RM) -r $(BUILDDIR) $(TARGET) $(TARGET-TD)

.PHONY: clean
