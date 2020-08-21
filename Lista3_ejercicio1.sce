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
G = syslin('c', A, B, C, D);
GN=clean(ss2tf(G));
// Para poder observar la funcion de transferencia, se debe colocar en console "GN"!!!!!!!!!!!!!!!!


