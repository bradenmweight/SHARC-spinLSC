***,Fulvene freqency calculation
memory,1000,m

file,2,Fulvene.wf,new;

basis=6-31++G**
gthresh,energy=1.d-8

geomtyp=xyz
symmetry,nosym
geometry={
12
	Fulvene CARTESIAN COORDINATES
C	0.13437	1.19472	-0.00000
C	1.41722	0.74879	0.00000
C	1.41713	-0.74891	0.00000
C	0.13423	-1.19470	-0.00000
H	-0.21822	2.23122	-0.00000
H	2.32763	1.36037	0.00000
H	2.32748	-1.36057	0.00000
H	-0.21853	-2.23114	-0.00001
C	-0.77571	0.00007	0.00000
C	-2.12821	0.00006	0.00000
H	-2.70629	0.93407	0.00000
H	-2.70628	-0.93415	0.00000
}

hf;accu,14
optg;coord,3n;

{frequencies,analytic
thermo,sym=auto
print,thermo}

mp2 ! {dfunc,pbe}
optg;coord,3n
{frequencies
thermo,
print,thermo}
put,molden,Fulvene.molden;