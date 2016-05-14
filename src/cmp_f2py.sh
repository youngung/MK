## F2py fortran libraries
f2py -c -m for_lib for.f

## Other items.
# compile func_fld_cy.pyx
python cmp_cython.py func_fld_cy --inplace

## plotting vm
python vm_check.py
