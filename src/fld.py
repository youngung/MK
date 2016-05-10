"""
*** NEWTON-RAPHSON AND OTHER METHOD - USE RECOMMENDED VALUES
1D-10			EPSILON (CONVERGENCE CRITERION, 1D-10)
100			ITERATIONS MAX (100)
*** INTEGRATION - USE RECOMMENDED VALUES
200000			MAXIMUM NUMBER OF STEPS (200000)
100			FREQUENCY OF OUTPUT (100)
0.001			STEP SIZE (0.001) 
1.5	 		EFFECTIVE STRAIN UPPER LIMIT (2.0)
*** MATERIAL - TITLE NEXT LINE
3    "Isotropic a=6" (1: VM, 2: YLD2000, 3: Hill48
"yld-alfas.txt"	YLD INPUT FILE NAME
1			STRAIN HARDENING LAW (1: Swift - 2: Voce - 3: Wh2000)
"Swift ***"     sigma = K*(eps+eps_0)**2
500.0			K
0.5			n
0.00001			eps_0
"Voce ***"  sigma = A - B0*EXP(-C*eps) + B1*eps
403.200000			A 		
273.400000			B0	
11.950000			C
191.300000			B1 
"KMEKB dislocation ***" See paper
0.107343E+03 	Hardening parameter 1 (sigma_0)
0.109382E+03	Hardening parameter 2 (sigma_y)
0.666667E+08	Hardening parameter 3 (F_y)
0.183326E+11	Hardening parameter 4 (F_s)
0.153798E+01	Hardening parameter 5 (k_2)
0.135387E+02	Hardening parameter 6 (k_1)
*** STRAIN RATE SENSITIVITY m=m0=constant when tau0 and f0 are both equal to 0.
0.00			tau0 Relaxation time associated with diffusion
0.00			f0 Maximum stress increase produced by DSA 
0.05			m0 initial strain rate sensitivity
1000D0			Strain rate ratio (E_a.dot / E_0.dot)
*** IMPERFECTION
0.999			IMPERFECTION PARAMETER F0 (0.996) 
*** FLD OPTION
5				OPTION (see below)
*** FOR OPTIONS 2 AND 4
0				FIRST ANGLE
90				SECOND ANGLE
5				ANGLE INCREMENT (in degrees)
### OPTIONS
1 ONE LINEAR STRAIN PATH (prompt), ONE ANGLE (prompt)
2 ONE LINEAR STRAIN PATH (prompt), SEVERAL ANGLES (as defined above) 
3 SEVERAL LINEAR STRAIN PATH (usually biaxial stretching), ONE ANGLE (0)
4 SEVERAL LINEAR STRAIN PATH (usually biaxial stretching), SEVERAL ANGLES
5 FULL FLD, ALL LINEAR STRAIN PATHS, ALL REASONABLE ANGLES
"""
##
EPS   = 1e-10
ITMAX = 100
NBPAS = 2e6
#FREQ = frequence of writing
DELTAT = 0.001 ## stepsize
TMAX   = 1.5   ## effective strain upper limit

ANG1   = 0.
ANG2   = 90.
ANGINC = 5.

## strain rate sensitivity parameter
P7  = 0.    # tau0 relaxatio time associated with diffusion
P8  = 0.    # f0 Maximum stress increase produced by DSA
P9  = 0.05  # m0 initial strain rate sensitivity
P10 = 1e3.  # strin rate ratio (e_a.dot / e_0.dot)
f0  = 0.999   ## impoerfaction parameter



