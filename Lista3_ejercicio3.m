clc;
clear;

w=logspace(-1,3,400);
x1=10*ones(1,60); x2=zeros(1,140); x3=-10*ones(1,200);
xtotal=[x1 x2 x3];
x4=[5*ones(1,100) 0*zeros(1,300)];
% Para lograr hallar la barrera de estabilidad se necesita:

semilogx(w, xtotal,'r')
ylim([-60,60])
grid;
hold on;
semilogx(w,x4,'b')
title('Barreras de Estabilidad')
xlabel('w(rad/s)')
ylabel('dB')
hold off;
