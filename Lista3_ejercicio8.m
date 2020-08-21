ap= [   -0.21   0.20    
         0.20  -0.21 ]

bp = [0.01   0    
      0        0.01]

cp = [    1   0
         0   1 ]

dp = 0*ones(2,2)


%*************************
%
% Natural Modes: Poles (Eigenvalues), Eigenvectors
%
[evec,eval] = eig(ap)   % evec contains eigenvectors
                        % eval contains poles or eigenvalues
                 
%*************************
%
% Unscaled Modal Analysis: We want to examine the natural modes (tendencies) of the F404 engine.
%
%
t = [0:0.1:12];
u = [0*t' 0*t'];         % Set input u to zero for all time in order to generate zero input response;
                         % i.e. response to an initial condition x_o.
                         % Transmission Zeros
%
%z = tzero(ss(ap,bp,cp,dp))                % transmission zeros
[p,z] = pzmap(ss(ap,bp,cp,dp))                % transmission zeros
% F404 Engine TRANSFER FUNCTIONS: From u_i to x_j
%
sys = zpk(ss(ap,bp,eye(2,2),0*ones(2,2))) % Zeros, Poles, and Gains fron u_i to x_j



%*************************
%
% Controllability 
%
cm = [bp ap*bp (ap^2)*bp];  % Controllability Matrix
rcm= rank(cm)           ;   % Rank of Controllability Matrix


%*************************
%
% Observability
%
om = [cp; 
      cp*ap;
      cp*(ap^2) ]     ;     % Observability Matrix
rom = rank(om)       ;      % Rank of Observability Matrix

%*************************
%
%
% F404 FREQUENCY RESPONSE: Unscaled Singular Values
%
% u = [ W_f (pph)  A_8 (sq in) ]
% x = [ N_2 (rpm)  N_25 (rpm)  T_45 (deg Fah) ]
% y = [ N_2 (rpm)              T_45 (deg Fah) ]
%
w = logspace(-2,3,100);
sv = sigma(ss(ap, bp, cp, dp),w);
sv = 20*log10(sv);

pause
% F404 SVD Analysis at w = 0.1 rad/sec
%
%s = j*0.1
%g = cp*inv(s*eye(3)-ap)*bp + dp
%[u, s, v ] = svd(g)

%*************************



[ns nc] = size(bp);                      % ns = number of inputs;  nc = number of controls;   
a = [ ap             bp
      0*ones(nc,ns)    0*ones(nc,nc) ]

b = [ 0*ones(ns,nc)
      eye(nc)      ]

c = [ cp  0*ones(nc,nc) ]

d = 0*ones(nc,nc)
sv = sigma(ss(a, b, c, d),w);
sv = 20*log10(sv);






% Design of Target Loop Singular Values Using Kalman Filter
%
ll =  inv(cp*inv(-ap)*bp + dp);     % Choose ll and lh to match singular values at all frequencies
lh = -inv(ap)*bp*ll;
l = [lh 
     ll];                           % ll, lh - for low and high frequency loop shaping

sv = sigma(ss(a, l, c, d),w);
sv = 20*log10(sv);

%%%%%%%%%%%%%%%%%%%%%%%%





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

%%%%%%%%%%%%%


% Augment Plant with Integrators at Plant Input and Plot Singular Values
%
[ns nc] = size(bp);                      % ns = number of inputs;  nc = number of controls;   
a = [ ap             bp
      0*ones(nc,ns)    0*ones(nc,nc) ];

b = [ 0*ones(ns,nc)
      eye(nc)      ];

c = [ cp  0*ones(nc,nc) ];

d = 0*ones(nc,nc);
sv = sigma(ss(a, b, c, d),w);
sv = 20*log10(sv);


% Recover Target Loop By Solving Cheap LQR Problem
%
q = c'*c;                                            % State Weighting Matrix
rho = 1e-3;                                          % Cheap control recovery parameter;
                                                     % The smaller the parameter, the better the recovery.
r = rho*eye(nc)                                      % Control Weigthing Matrix
%[k, poles, g, rr] = care(a,b,q,r);                   % Compute Control Gain Matrix G
[k, poles, g] = care(a,b,q,r);                   % Compute Control Gain Matrix G





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




