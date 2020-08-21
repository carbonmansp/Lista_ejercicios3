// CONTROL INTERCAMBIADOR DE CALOR
// Alumno: Anthony Alvarez Oviedo
// Autor del codigo: juan C. Cutipa-Luque


//Codigo de un ejemplo Mostrado en Clase
clf();         // close current figure
clear          // clear all pasta variables
xdel(winsid()) // close all windows

A=    [-0.21,0.20;
 0.20,-0.21];
 
 
B=[0.01,0;
 0,0.01];

C=[1,0;
 0,1]; 

D=[0,0;
 0,0];
 
ap_s=A;  
bp_s=B;
cp_s=C;
dp_s=D;

// Controllability and Observability
// Cc=[B, AB, A^2 B,..., A^(n-1) B]
Cc = cont_mat(ap_s,bp_s)
rankCc=rank(Cc)
//
// O=[C; CA; CA^2;...; CA^(n-1) ]
O = obsv_mat(ap_s, cp_s)
rankO=rank(O)
// verify if the rank of Cc is n, dimension of a
// verify if the rank of O is n, dimension of a

//      Singular values of LTI the model          //
G = syslin('c', ap_s, bp_s, cp_s, dp_s);

w = logspace(-3,3);
sv = svplot(G,w);

scf(0);
plot2d("ln", w, 20*log(sv')/log(10))
xgrid(12)
xtitle("Grafica de Valores Singulares","Frecuencia (rad/s)", "Amplitud (dB)");


ms=1.6;// 0.3;%1.5;    % guarantee overshot Mp < 6dB = 20*log10(2) 
wbs=0.23;//0.05;%0.23;
ee=1e-3;//1e-4
ki=1; // used to give more accurate adjustment to the cut-off frequency wbs
      // by default set it to 1
//           --------     WT Data    ------------
mt=1.3;//1.00;    % guarantee overshot Mp < 2dB = 20*log10(1.26)
wbt=4.1;//9.1;%4.1;
ee=1e-3;//1e-4

//           --------     WS     ------------

s=poly(0,'s');
ws1=(s/ms+wbs)/(s+wbs*ee),
ws2=ws1;
ws=[ws1,0;0,ws2]
//Ws=syslin('c',ws)
Ws=blockdiag(ws1,ws2)

//           --------     WT     ------------

s=poly(0,'s');
wt1=(s+wbt/mt)/(ee*s+wbt),
wt2=wt1;
wt=[wt1,0;0,wt2]
Wt=syslin('c',wt)
//Wt=blockdiag(wt1,wt2)


//           --------     WR     ------------
s=poly(0,'s');
wr1=s/s,
wr2=wr1;
wr=[wr1,0;0,wr2]

// ------------------ Plot weighting functions

svs = svplot(Ws,w);
svt = svplot(Wt,w);
scf(2);
plot2d("ln", w, [-20*log(svs')/log(10) -20*log(svt')/log(10)])
xgrid(12)
xtitle("Singular values plot inv(Ws) and inv(Wt)","Frequency (rad/s)", "Amplitude (dB)");
