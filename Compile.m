%%%%% Run.m %%%%%%

% Initial cleanup
clear all; format short; format compact; beep on; close all;
warning('on','all');

%      [1      2      3      4      5      6      7      8      9      10     11     12 
%      [theta  rhoZ   muZ    sigZ   lamb1k lamb2k lamb3k beta   tauC   pi     MQent  SQent ] 
over = [0.607  0.956  0.704  0.375  0.700  0.807  0.446  0.861  0.355  0.095  2.000  1.000 ];
fid  = nan; % {1 = screen, NaN = file}

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
