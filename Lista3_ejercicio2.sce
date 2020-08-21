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


//Plot Tansmission Zeros
scf(1)
plzr(A,B,C,D)
