#Standalone makefile
PACKAGE = RhoNNO

LIBA	= lib$(PACKAGE).a
LIBSO	= lib$(PACKAGE).so

# Lists of objects to include in library:
NAMELIST = VNeuralNetPlotter VNeuralNet VSupervisedNet VUnsupervisedNet TFD TSGCS TSGNG TPerceptron TMLP TXMLP \
           TNeuralNetCell TGCS TGNG TLVQ TDataServe TNNK

HDRS    = $(addsuffix .h, $(NAMELIST) )
OBJS    = $(addsuffix .o, $(NAMELIST) )
#ROOTCINTTARGETS := $(filter-out $(PACKAGE)_LinkDef.rdl, $(wildcard *.rdl)) $(wildcard $(PACKAGE)_LinkDef.rdl)
#ROOTCINTHDRS :=   $(addsuffix .h, $(basename $(ROOTCINTTARGETS)) )
ROOTCINTHDRS := $(HDRS) $(PACKAGE)_LinkDef.h

# Default action 
all : $(ROOTCINTHDRS) $(LIBA) $(LIBSO)

hdr: $(ROOTCINTHDRS)

shlib: $(LIBSO)

# Static library:
$(LIBA) : $(OBJS)  $(PACKAGE)Cint.o
	rm -f $@
	$(AR) $(ARFLAGS)  $(LIBDIR)/$@ $^

# Shared library:
$(LIBSO) : $(OBJS)  $(PACKAGE)Cint.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(GLIBS) -o $(LIBDIR)/$(LIBSO) $(EXPLLINKLIBS)

# Rules for Dictionary:
$(PACKAGE)Cint.o : $(PACKAGE)Cint.cxx 
	$(CXX) $(CXXFLAGS) -c $<
	@echo "-----$(PACKAGE)Cint.o---------"

$(PACKAGE)Cint.cxx : $(ROOTCINTHDRS)
	rootcint -f $@ -c  -I.. $^
	@echo "-----$(PACKAGE)Cint.cxx---------"

# Static pattern rule for object file dependency on sources:
$(OBJS) : %.o : %.cxx  %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

#$(ROOTCINTHDRS) : %.h : %.rdl
#	cp $< $@

# Some extra dependecies from #includes (needs to be completed,
# if in doubt do a "gmake clean" first):

#
# Complete cleanup:
#
clean :
	rm  -f *.o 
	rm  -f $(PACKAGE)Cint.*
	rm  -f $(LIBA) $(LIBSO)


