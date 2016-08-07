#!/bin/bash

## PBS Arugments
#PBS -m abe
#PSB -M youngungj@gmail.com
#PBS -N MK
#PBS -l select=1:ncpus=24:mem=24gb,walltime=48:00:00
#PBS -l place=free
#PBS -j eo

cd ~/repo/mk/mk/

## VPSC-Hardening Function

## Palmetto
fnpickle=/home/younguj/repo/mk/matDatabase/IFsteel/yf_collection.dll
fnpickle_vpsc_hard=/home/younguj/repo/mk/matDatabase/IFsteel/vpsc_FC_rhop.dll

## Mac
fnpickle=~/repo/mk/matDatabase/IFsteel/yf_collection.dll
fnpickle_vpsc_hard=~/repo/mk/matDatabase/IFsteel/vpsc_FC_rhop.dll


## running conditions
r0=-0.5
r1=2.5
nr=15

#  f0s=("0.980" "0.985" "0.990" "0.992" "0.994"  "0.996")
f0s=("0.980")
# f0s=("0.985")
# f0s=("0.990")
# f0s=("0.992")
# f0s=("0.994")
# f0s=("0.996")
for f0 in "${f0s[@]}"
do
    python mk_run_pickles.py --f0 $f0 --r0 $r0 --r1 $r1 --nr $nr --fnpickle $fnpickle  --fnpickle_vpsc_hard $fnpickle_vpsc_hard
done


# python mk_run_pickles.py --f0 0.990 --r0 0 --r1 2.5 --nr 31 --fnpickle /home/younguj/repo/mk/matDatabase/IFsteel/yf_collection.dll
# python mk_run.py --f0 0.995 --r0 -0.5 --r1 2.5 --nr 31 --mat 6  ## H48R
# cp minFLD.txt minFLD_H48R.txt
# cp mk_fld_pp.pdf mk_fld_pp_HH48R.pdf

# python mk_run.py --f0 0.995 --r0 -0.5 --r1 2.5 --nr 31 --mat 7  ## H48YR0
# cp minFLD.txt minFLD_H48YR0.txt
# cp mk_fld_pp.pdf mk_fld_pp_HH48YR0.pdf

# python mk_run.py --f0 0.995 --r0 -0.5 --r1 2.5 --nr 31 --mat 8  ## H48YR90
# cp minFLD.txt minFLD_H48YR90.txt
# cp mk_fld_pp.pdf mk_fld_pp_HH48YR90.pdf

# python mk_run.py --f0 0.995 --r0 -0.5 --r1 2.5 --nr 31 --mat 9  ## HH48R_vpsc
# cp minFLD.txt minFLD_H48R_vpsc.txt
# cp mk_fld_pp.pdf mk_fld_pp_HH48R_vpsc.pdf

# python mk_run.py --f0 0.995 --r0 -0.5 --r1 2.5 --nr 31 --mat 10 ## HH48YR0_vpsc
# cp minFLD.txt minFLD_H48YR0_vpsc.txt
# cp mk_fld_pp.pdf mk_fld_pp_HH48YR0_vpsc.pdf

# python mk_run.py --f0 0.995 --r0 -0.5 --r1 2.5 --nr 31 --mat 11 ## HH48YR90_vpsc
# cp minFLD.txt minFLD_H48YR90_vpsc.txt
# cp mk_fld_pp.pdf mk_fld_pp_HH48YR90_vpsc.pdf

# python mk_run.py --f0 0.995 --r0 -0.5 --r1 2.5 --nr 31 --mat 12 ## yld1
# cp minFLD.txt minFLD_yld1.txt
# cp mk_fld_pp.pdf mk_fld_pp_yld1.pdf

# python mk_run.py --f0 0.995 --r0 -0.5 --r1 2.5 --nr 31 --mat 13 ## yld2
# cp minFLD.txt minFLD_yld2.txt
# cp mk_fld_pp.pdf mk_fld_pp_yld2.pdf
