## F2py fortran libraries
f2py -c -m for_lib for.f

## Other items.
# compile func_fld_cy.pyx
python cmp_cython.py build_ext --inplace



## link yld2000
gfortran -shared -fPIC -o yld2000_lib.so yld2000_lib.f
#f2py --opt='-O3' -c -m yld2000 yld2000.f yld2000_lib.so
#f2py  -c -m yld2000 yld2000.f yld2000_lib.so




## plotting vm
python vm_check.py
