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
// acceleration, ay, and the rudder deflection,?.
//H(5,3)

//Plot Tansmission Zeros
//scf(1)
//plzr(A,B,C,D)
//Calculate and Plot of Singular Values
G = syslin('c', A, B, C, D);
GN=clean(ss2tf(G));
w = logspace(-2,2);
scf(2)
sv = svplot(G,w);
plot2d("ln", w, [20*log(sv')/log(10)],leg="Max Singular Values@Min Singular Values")
//LQR
Q=diag(1,1); //Weights on states
R=1; //Weight on input
//[K,X]=lqr(G,Q,R);

ap=A;
bp=B;
cp=C;
dp=D;
//PLANTA CON INTEGRADOR//

[ns,nc]=size(bp); //ns= numero de entradas; nc=numero de controles
Ai=[ap             bp;
    0*ones(nc,ns) 0*ones(nc,nc)];

Bi=[0*ones(ns,nc); eye(nc,nc)];
    
Ci=[cp 0*ones(nc,nc)];

Di=0*ones(nc,nc);

sysi=syslin('c',Ai,Bi,Ci,Di);

I=eye(nc);
q = Ci'*Ci;          //Matriz de ponderación del estado
rho = 1e-9;        //Parámetro de recuperación de control barato

r = rho*eye(nc,nc)                 //Matriz de ponderación de control
         //how we calculate B
B=Bi*inv(r)*Bi';
A=Ai;
        //Solv the ricatti equation
X=riccati(A,B,q,'c','eigen');
        //matriz de ganacia
G_1=inv(r)*Bi'*X;
