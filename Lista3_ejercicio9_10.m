
clear
clc
ap= [   -0.21   0.20    
         0.20  -0.21 ]

bp = [0.01   0    
      0        0.01]

cp = [    1   0
         0   1 ]

dp = 0*ones(2,2)


w = logspace(-2,3,100);
sv = sigma(ss(ap, bp, cp, dp),w)
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Valores singulares de la planta')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


%
% Scaling Matrices
%
%  unew = su uold
%  xnew = sx xold
%  ynew = sy yold
%
%su = diag( [1/110, 1/22] )
%sx = diag( [1/250,  1/28] ) %1/350, 
%sy = diag( [1/250,  1/28] )

%
% Scaled F404 Dynamics
%
%
% u = [ W_f/110    A_8/22 ]
% x = [ N_2/250  N_25/350  T_45/28 ]
% y = [ N_2/250            T_45/28 ]
%
%ap = sx*ap*inv(sx)
%bp = sx*bp*inv(su)
%cp = sy*cp*inv(sx)
%dp = sy*dp*inv(su)
    
%


%
% Augment Plant with Integrators at Plant Input and Plot Singular Values
%
[ns nc] = size(bp);                      % ns = number of inputs;  nc = number of controls;   
a = [ ap             bp
      0*ones(nc,ns)    0*ones(nc,nc) ]

b = [ 0*ones(ns,nc)
      eye(nc)      ]

c = [ cp  0*ones(nc,nc) ]

d = 0*ones(nc,nc)
sv = sigma(ss(a, b, c, d),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Design Plant Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%return


%
%*****
%


%
% Design of Target Loop Singular Values Using Kalman Filter
%
ll =  inv(cp*inv(-ap)*bp + dp);     % Choose ll and lh to match singular values at all frequencies
lh = -inv(ap)*bp*ll;
l = [lh 
     ll];                           % ll, lh - for low and high frequency loop shaping

sv = sigma(ss(a, l, c, d),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Filter Open Loop (G_{FOL}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause



pnint = eye(nc)                                    % Process Noise Intensity Matrix
mu = 0.01;                                         % Measurement Noise Intesity; Used to adjust Kalman Filter Bandwidth
                                                   % Small mu - expensive sensor   - large bandwidth
                                                   % Large mu - inexpensive sensor - small bandwidth
mnint = mu*eye(nc)                                 % Measurement Noise Intensity Matrix 
% sysKal=ss(a, [b l], c, [d 0*ones(nc,nc)]);
% [kest, h, sig]= kalman(sysKal,pnint, mnint);  % Compute Filter Gain Matrix h
%[sig, poles, g1, rr] = care(a',c',l*l', mnint);                          
[sig, poles, g1] = care(a',c',l*l', mnint);                          
% Alternate Method for Computing h
h = g1';
sv = sigma(ss(a, h, c, d),w);
tsv = 20*log10(sv);
semilogx(w, tsv)
%clear sv
title('Target Loop (G_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


sv = sigma(ss(a-h*c, h, -c, eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Target Sensitivity (S_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


sv = sigma(ss(a-h*c, h, c, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv, w, 20*log10(10./w))
%clear sv
title('Target Complementary (T_{KF}) Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


%
% Recover Target Loop By Solving Cheap LQR Problem
%
q = c'*c;                                            % State Weighting Matrix
rho = 1e-3;                                          % Cheap control recovery parameter;
                                                     % The smaller the parameter, the better the recovery.
r = rho*eye(nc)                                      % Control Weigthing Matrix
%[k, poles, g, rr] = care(a,b,q,r);                   % Compute Control Gain Matrix G
[k, poles, g] = care(a,b,q,r);                   % Compute Control Gain Matrix G


%
% Compensator Analysis
%
ak = [ a-b*g-h*c  0*ones(ns+nc,nc)
       g          0*ones(nc,nc) ]

bk = [ h
       0*ones(nc,nc) ]

ck = [0*ones(nc, ns+nc) eye(nc,nc) ]


%cpoles = eig(ak)                               % Compensator Poles
%czeros = tzero(a, h, g, 0*ones(nc,nc))         % Compensator Zeros
[cpoles, czeros] = pzmap(ss(a, h, g, 0*ones(nc,nc)))         % Compensator Zeros
%zerocheck = tzero(ak, bk, ck, 0*ones(nc,nc))   % Check Compensator Zeros
[polecheck, zerocheck] = pzmap(ss(ak, bk, ck, 0*ones(nc,nc)))   % Check Compensator Zeros

sv = sigma(ss(ak, bk, ck, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Compensator Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause


%
% Open Loop Analysis
%
al = [ ap                     bp*ck
       0*ones(ns+nc+nc,ns)    ak    ]

bl = [ 0*ones(ns,nc)
       bk ]
    
cl = [ cp  0*ones(nc,ns+nc+nc) ]
    
%olpoles = eig(al)                          % Open Loop Poles
%olzeros = tzero(al,bl,cl,0*ones(nc,nc))    % Open Loop Zeros
[olpoles, olzeros] = pzmap(ss(al,bl,cl,0*ones(nc,nc)))    % Open Loop Zeros    
sv = sigma(ss(al, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv, w, tsv)
%clear sv
title('Open Loop Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause   
%%%%%%%%%%%%%% Para hallar las curvas de sensibilidad S y T
sv = sigma(ss(al-bl*cl, bl, -cl, eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause 

sv = sigma(ss(al-bl*cl, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);
semilogx(w, sv)
%clear sv
title('Complementary Sensitivity Singular Values')
grid
xlabel('Frequency (rad/sec)')
ylabel('Singular Values (dB)')
pause

%
% Closed Loop Analysis
%
clpoles = eig(al-bl*cl)           % Closed Loop Poles
clpkf = eig(a - h*c)              % Closed Loop Poles Due to Kalman Filter
clpreg = eig(a - b*g)             % Closed Loop Poles Due to Regulator

%
% Step Response in Closed Loop 
%

[y,t] = step(ss(al-bl*cl, bl, cl, 0*eye(nc)));
plot(t,y(:,1,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')
pause
%return

%
plot(t,y(:,1,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 2')
pause
%return

plot(t,y(:,2,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 response caused by input 1')
pause
%return

%
plot(t,y(:,2,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 response caused by input 2')
pause
%return
disp('You can see a good tracking and a good disturbance rejection')
