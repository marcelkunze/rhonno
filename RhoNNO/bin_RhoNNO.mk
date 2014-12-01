# Please see SoftRelTools/HOWTO-dependency for documentation
# $Id: bin_RhoNNO.mk,v 1.2 2001-05-16 15:45:43 marcel Exp $
#++ NON-STANDARD
        LINKLISTDEPENDS += [LINK_RhoNNO, $(LINK_RhoNNO)]

# Stop gap solution to link RHO (it is not yet in PackageList)
include $(TOPDIR)/RhoMakefiles/link_Rho.mk

override LINK_RHO	+= $(PACKAGE)GNUmakefile
override LINK_ROOT	+= $(PACKAGE)GNUmakefile

