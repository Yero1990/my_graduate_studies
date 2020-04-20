# make file for F2PY and installation into the utilities 
#
# default procedure for compiling a .f file
#
#suffix rule necessary for compiling with absoft
.f.so :
	$(F2PY) $(F2PYFLAGS) $<

HERE          = /data/boeglin.1/user/deut/crosec

CROSEC        = /data/boeglin.1/user/deut/crosec

CERNLIB       =

DEST          = /Users/boeglinw/Documents/Python/Utilities


PRINT	      = 

F2PY          = f2py
F2PYFLAGS     = -c --build-dir ./build/ --f90flags=-m32 --f77flags=-m32 -m

PROGRAM	      = sig_deutpwia

SRCS          = sig_deutpwia.f relkin.f dpolint.f dlocate.f \
		subcc1.f drho_paris.f

all:		$(DEST)/$(PROGRAM).so

$(PROGRAM):     $(SRCS)
		$(F2PY) $(F2PYFLAGS) $(PROGRAM)  $(SRCS)
		touch $(PROGRAM)

clean:	
		rm $(PROGRAM).so
		rm $(PROGRAM)

# copy to the instalation directory 

$(DEST)/$(PROGRAM).so:  $(PROGRAM)
		cp $(PROGRAM).so $(DEST)/$(PROGRAM).so










