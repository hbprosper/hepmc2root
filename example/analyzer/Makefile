#----------------------------------------------------------------------------
# Description: Makefile to build analyzers
# Created:     Mon Dec  4 16:20:33 2017 by mkanalyzer.py
# Author:      Shakespeare's ghost
#----------------------------------------------------------------------------
ifndef ROOTSYS
$(error *** Please set up Root)
endif

name    := analyzer

# Sub-directories
srcdir	:= src
tmpdir	:= tmp
libdir	:= lib
incdir	:= include

$(shell mkdir -p tmp)
$(shell mkdir -p lib)

# Set this equal to the @ symbol to suppress display of instructions
# while make executes
ifdef verbose
AT 	:=
else
AT	:= @
endif
#-----------------------------------------------------------------------
# sources and objects
#-----------------------------------------------------------------------
header  := $(incdir)/tnm.h
linkdef := $(incdir)/linkdef.h
cinthdr := $(srcdir)/dictionary.h
cintsrc := $(srcdir)/dictionary.cc

# Construct list of sources to be compiled into applications
appsrcs	:= $(wildcard *.cc)
appobjects	:= $(addprefix $(tmpdir)/,$(appsrcs:.cc=.o))

# Construct list of applications
applications := $(appsrcs:.cc=)

# Construct list of sources to be compiled into shared library
ccsrcs	:= $(filter-out $(cintsrc),$(wildcard $(srcdir)/*.cc))
sources	:= $(ccsrcs) $(cintsrc)
objects	:= $(subst $(srcdir)/,$(tmpdir)/,$(sources:.cc=.o))

# Display list of applications to be built
#say	:= $(shell echo "appsrcs:     $(appsrcs)" >& 2)
#say	:= $(shell echo "appobjects:  $(appobjects)" >& 2)
#say	:= $(shell echo "objects:     $(objects)" >& 2)
#$(error bye!) 

#-----------------------------------------------------------------------
# 	Define which compilers and linkers to use
#-----------------------------------------------------------------------
# If clang++ exists use it, otherwise use g++
COMPILER:= $(shell which clang++)
ifneq ($(COMPILER),)
CXX     := clang++
LINK	:= clang++
else
CXX     := g++
LINK	:= g++
endif 
CINT	:= rootcint

#-----------------------------------------------------------------------
# 	Define paths to be searched for C++ header files (#include ....)
#-----------------------------------------------------------------------
CPPFLAGS:= -I. -I$(incdir) -I$(srcdir) $(shell root-config --cflags) 

# 	Define compiler flags to be used
#	-c		perform compilation step only 
#	-g		include debug information in the executable file
#	-O2		optimize
#	-ansi	require strict adherance to C++ standard
#	-Wall	warn if source uses any non-standard C++
#	-pipe	communicate via different stages of compilation
#			using pipes rather than temporary files

CXXFLAGS:= -c -g -O2 -ansi -Wall -pipe -fPIC

#	C++ Linker
#   set default path to shared library
LD	:= $(LINK) -Wl,-rpath,$(PWD)/$(libdir)

OS	:= $(shell uname -s)
ifeq ($(OS),Darwin)
    LDSHARED	:= $(LD) -dynamiclib
    LDEXT       := .dylib
else
    LDSHARED	:= $(LD) -shared
    LDEXT       := .so
endif

#	Linker flags

LDFLAGS := -g

# 	Libraries

LIBS	:=  $(shell root-config --libs) -L$(libdir) -lMinuit -lMathCore

sharedlib := $(libdir)/libtnm$(LDEXT)

#-----------------------------------------------------------------------
#	Rules
#	The structure of a rule is
#	target : source
#		command
#	The command makes a target from the source. 
#	$@ refers to the target
#	$< refers to the source
#-----------------------------------------------------------------------
all:	$(sharedlib) $(applications) 

bin:	$(applications)

lib:	$(sharedlib)

# Syntax:
# list of targets : target pattern : source pattern

# Make applications depend on shared libraries to force the latter
# to be built first

$(applications)	: %	: $(tmpdir)/%.o  $(sharedlib)
	@echo "---> Linking $@"
	$(AT)$(LD) $(LDFLAGS) $< $(LIBS) -ltnm -o $@

$(appobjects)	: $(tmpdir)/%.o	: %.cc
	@echo "---> Compiling application `basename $<`" 
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ # >& $*.FAILED
	@rm -rf $*.FAILED

$(sharedlib)	: $(objects)
	@echo "---> Linking `basename $@`"
	$(AT)$(LDSHARED) $(LDFLAGS) -fPIC $(objects) $(LIBS) -o $@

$(objects)	: $(tmpdir)/%.o	: $(srcdir)/%.cc
	@echo "---> Compiling `basename $<`" 
	$(AT)$(CXX) $(CXXFLAGS) $(CPPFLAGS)  $< -o $@ # >& $*.FAILED
	$(AT)rm -rf $*.FAILED

$(cintsrc)  : $(header) $(linkdef)
	@echo "---> Generating dictionary `basename $@`"
	$(AT)$(CINT) -f $@ -c -I. -Iinclude -I$(ROOTSYS)/include $+
	$(AT)mv $(srcdir)/*.pcm $(libdir)

# 	Define clean up rules
clean   :
	rm -rf $(tmpdir)/* $(libdir)/* $(srcdir)/dictionary* $(applications)
