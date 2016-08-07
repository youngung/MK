## F2py fortran libraries
f2py -c -m yf_for for.f
## link yld2000 library file
gfortran -shared -fPIC -o yld2000_lib.so yld2000_lib.f
## yld2000 minimalistic version wrapped by f2py
f2py -c -m yld2000 yld2000.f yld2000_lib.so
## yld2000 stand-alone executable
gfortran yld2000.f yld2000_lib.f -o yld2000_sa
