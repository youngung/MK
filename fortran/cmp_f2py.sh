## F2py fortran libraries
f2py -c -m yf_for for.f

## yld2000 stand-alone executable
gfortran yld2000.f yld2000_lib.f -o yld2000_sa



# ## link yld2000 library file
gfortran -shared -fPIC -o yld2000_lib.so yld2000_lib.f
# ## yld2000 minimalistic version wrapped by f2py
f2py -c -m yf_yld2000 yld2000.f yld2000_lib.so


## One might want to copy
# yf_for.so; yld2000_sa; yf_yld2000.so; yld2000_lib.so

## blow is manual-copy used in Pal
cp yf_for.so /home/younguj/.local/lib/python2.7/site-packages/
cp yld2000_sa /home/younguj/.local/lib/python2.7/site-packages/
cp yld2000_lib.so /home/younguj/.local/lib/python2.7/site-packages/
cp yf_yld2000.so /home/younguj/.local/lib/python2.7/site-packages/



## one might want to change LD_LIBRARY_PATH to use the 'shared' library object yld2000_lib.so when yf_yld2000 is loaded.
## I added below line to my .bashrc
## export LD_LIBRARY_PATH="/home/younguj/.local/lib/python2.7/site-packages/:$LD_LIBRARY_PATH"
