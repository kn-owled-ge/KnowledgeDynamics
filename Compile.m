%%%%% Run.m %%%%%%

% Initial cleanup
clear all; format short; format compact; beep on; close all;
warning('on','all');

%      [1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16     17    ]
%      [theta  omega  rhoZ   muZ    sigZ   delK   lambK  delN   lambN  alpha  beta   gama   tauC   tauN   Pext   MQent  SQent ] 
over = [0.728  0.500  0.933  0.000  0.384  0.110  0.377  0.300  0.400  0.001  0.892  0.070  0.350  0.000  0.095  2.000  1.000 ];
fid  = 1; % {1 = screen, NaN = file}

%%%%% ModelF.m %%%%%%

% Calibrate model
[P, Qz, Qx, Qw, k, kp, n, np, z, x, w, fail] = calibF(over,fid);

%%%%% MainF.m %%%%%%

% initializing the value function
if (exist(['Saves\Model ' P.sMod '\vfi3Mat.mat'],'file') && P.Vload)
   load (['Saves\Model ' P.sMod '\vfi3Mat.mat']);
else
   V = -10.*ones(P.nk*P.nn,P.nz*P.nx*P.nw);
end

% initialize panels
kpd  = ones(P.Ns, 1).*k(ceil(P.nk/2));
npd  = ones(P.Ns, 1).*n(ceil(P.nn/2));
zpd  = P.muZ + (P.sigZ/sqrt(1-P.rhoZ^2)).*randn(P.Ns,1);
%xpd = P.muX + (P.sigX/sqrt(1-P.rhoX^2)).*randn(P.Ns,1);
%wpd = P.muW + (P.sigW/sqrt(1-P.rhoW^2)).*randn(P.Ns,1);
xpd  = zeros(P.Ns,1);
wpd  = zeros(P.Ns,1);
