clear
clc
ap= [   -0.21   0.20    
         0.20  -0.21 ]

bp = [0.01   0    
      0        0.01]

cp = [  1   0
         0   1 ]

dp = 0*ones(2,2)  


w = logspace(-2,3,100);
sv = sigma(ss(ap, bp, cp, dp),w);
sv = 20*log10(sv);


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
%return

%
% Design of Target Loop Singular Values Using Kalman Filter
%
ll =  inv(cp*inv(-ap)*bp + dp);     % Choose ll and lh to match singular values at all frequencies
lh = -inv(ap)*bp*ll;
l = [lh 
     ll];                           % ll, lh - for low and high frequency loop shaping

pnint = eye(nc)                                    % Process Noise Intensity Matrix
mu = 0.01;                                         
mnint = mu*eye(nc)                                              
[sig, poles, g1] = care(a',c',l*l', mnint);                          
% Alternate Method for Computing h
h = g1';
sv = sigma(ss(a, h, c, d),w);
tsv = 20*log10(sv);


sv = sigma(ss(a-h*c, h, -c, eye(nc)),w);
sv = 20*log10(sv);


sv = sigma(ss(a-h*c, h, c, 0*eye(nc)),w);
sv = 20*log10(sv);


%
% Recover Target Loop By Solving Cheap LQR Problem
%
q = c'*c;                                            % State Weighting Matrix
rho = 1e-3;                                         
r = rho*eye(nc)                                      % Control Weigthing Matrix
[k, poles, g] = care(a,b,q,r);                % Compute Control Gain Matrix G
%
% Compensator Analysis
%
ak = [ a-b*g-h*c  0*ones(ns+nc,nc)
       g          0*ones(nc,nc) ]

bk = [ h
       0*ones(nc,nc) ]

ck = [0*ones(nc, ns+nc) eye(nc,nc) ]

sv = sigma(ss(ak, bk, ck, 0*eye(nc)),w);
sv = 20*log10(sv);


%
% Open Loop Analysis
%
al = [ ap                     bp*ck
       0*ones(ns+nc+nc,ns)    ak    ]

bl = [ 0*ones(ns,nc)
       bk ]
    
cl = [ cp  0*ones(nc,ns+nc+nc) ]
sv = sigma(ss(al-bl*cl, bl, -cl, eye(nc)),w);
sv = 20*log10(sv);


sv = sigma(ss(al-bl*cl, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);

   
%olpoles = eig(al)                          % Open Loop Poles
%olzeros = tzero(al,bl,cl,0*ones(nc,nc))    % Open Loop Zeros
[olpoles, olzeros] = pzmap(ss(al,bl,cl,0*ones(nc,nc)))    % Open Loop Zeros    
sv = sigma(ss(al, bl, cl, 0*eye(nc)),w);
sv = 20*log10(sv);

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
subplot(2,2,1)
plot(t,y(:,1,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 1')

subplot(2,2,2)
plot(t,y(:,1,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 1 response caused by input 2')

subplot(2,2,3)
plot(t,y(:,2,1))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 response caused by input 1')

subplot(2,2,4)
plot(t,y(:,2,2))
xlabel('Time (s)')
ylabel('Amplitude')
title('output 2 response caused by input 2')
