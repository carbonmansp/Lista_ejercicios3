clc
clear all
A=[-0.21,0.20;
 0.20,-0.21]; 

B=[0.01,0;
 0,0.01]; 

C=[1,0;
 0,1]; 

D=[0,0;
 0,0]; 

//U=[xi;eta;zeta]; 
s=poly(0,"s");
//S=s*eye(6))
S=[s,0; 0,s]; 

// The transfer-function of the system is: 
// H(s)= C(SI-A)^-1*B+D

// Calculate (SI-A)^*1
SIAinv= inv(S-A)
// Calculate the transfer-function
H=C*SIAinv*B +D
// Display the transfer-function between the lateral

//Calculate and Plot of Singular Values
G = syslin('c', A, B, C, D);
GN=clean(ss2tf(G));
w = logspace(-2,2);
scf(2)
sv = svplot(G,w);
plot2d("ln", w, [20*log(sv')/log(10)],leg="Max Singular Values@Min Singular Values")
xtitle('Valores Singulares','w(rads/s)','dB')
//En la primera fila se toma en cuenta el valor maximo con max(sv(1)), mientras que de la segunda fila se considera el valor minimo con min(sv(2)).
