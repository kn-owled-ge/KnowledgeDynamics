% ----------------------------------------------------
% MAINF:
%   Construct the industry equilibrium by applying Krusell and Smith [1998] method
% 
%   h' = alp.b0 + alp.b1*(h-alp.b0)
%
% By Robert Parham, 2015
% Based on code by Lu Zhang, 2001
%
% ----------------------------------------------------
function [V, optK, optN, kpd, npd, zpd, xpd, wpd, fail] = mainF(P, Qz, Qx, Qw, k, kp, n, np, z, x, w)

% initializing the value function
if (exist('vfi3Mat.mat','file') && P.Vload)
   load ('vfi3Mat.mat');
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

% for given parameters solve for the value function and optimal decision rules
[V, optK, optN, vfail] = vfiCC_mex(P, Qx, Qz, Qw, k, kp, n, np, z, x, w, V);
if (vfail)
    fprintf(P.fid,'\nVFI failure! Aborting.\n\n');
    fail=1;
    return;
end

% simulate the economy over a long time period
[Kft, Nft, Zft, Xft, Wft, Vft, Gft, kpd, npd, zpd, xpd, wpd, sfail] = simCC_mex(P, k, n, z, x, w, V, optK, optN, kpd, npd, zpd, xpd, wpd, P.Ts1);
if (sfail)
    fprintf(P.fid,'\nSimulation failure! Aborting.\n\n');
    fail=1;
    return;
end

if (P.SIMfin)
   % visualize
   eps     = 0.05;
   tot     = double(P.Ns*(P.Ts1-P.Tb));
   exiters = sum(sum(Gft(:,P.Tb+1:end)~=0))/tot;
   minVft  = min(min(Vft(:,P.Tb+1:end)));
   maxVft  = max(max(Vft(:,P.Tb+1:end)));
   atMinK  = sum(sum(Kft(:,P.Tb+1:end)<(k(1  )*(1+eps))))/tot;
   atMaxK  = sum(sum(Kft(:,P.Tb+1:end)>(k(end)*(1-eps))))/tot;
   atMinN  = sum(sum(Nft(:,P.Tb+1:end)<(n(1  )*(1+eps))))/tot;
   atMaxN  = sum(sum(Nft(:,P.Tb+1:end)>(n(end)*(1-eps))))/tot;
   atMinZ  = sum(sum(Zft(:,P.Tb+1:end)<(z(1  )*(1+eps))))/tot;
   atMaxZ  = sum(sum(Zft(:,P.Tb+1:end)>(z(end)*(1-eps))))/tot;
   atMinX  = sum(sum(Xft(:,P.Tb+1:end)<(x(1  )*(1+eps))))/tot;
   atMaxX  = sum(sum(Xft(:,P.Tb+1:end)>(x(end)*(1-eps))))/tot;
   atMinW  = sum(sum(Wft(:,P.Tb+1:end)<(w(1  )*(1+eps))))/tot;
   atMaxW  = sum(sum(Wft(:,P.Tb+1:end)>(w(end)*(1-eps))))/tot;
   atMinV  = sum(sum(Vft(:,P.Tb+1:end)<(minVft*(1+eps))))/tot;
   atMaxV  = sum(sum(Vft(:,P.Tb+1:end)>(maxVft*(1-eps))))/tot;

   fprintf(P.fid,'       min | at min |   max  | at max |   mean \n');
   fprintf(P.fid,'K: %8.3f %8.3f %8.3f %8.3f %8.3f\n', min(min(Kft(:,P.Tb+1:end))), atMinK, max(max(Kft(:,P.Tb+1:end))), atMaxK, mean(mean(Kft(:,P.Tb+1:end))));
   fprintf(P.fid,'N: %8.3f %8.3f %8.3f %8.3f %8.3f\n', min(min(Nft(:,P.Tb+1:end))), atMinN, max(max(Nft(:,P.Tb+1:end))), atMaxN, mean(mean(Nft(:,P.Tb+1:end))));
   fprintf(P.fid,'Z: %8.3f %8.3f %8.3f %8.3f %8.3f\n', min(min(Zft(:,P.Tb+1:end))), atMinZ, max(max(Zft(:,P.Tb+1:end))), atMaxZ, mean(mean(Zft(:,P.Tb+1:end))));
   fprintf(P.fid,'X: %8.3f %8.3f %8.3f %8.3f %8.3f\n', min(min(Xft(:,P.Tb+1:end))), atMinX, max(max(Xft(:,P.Tb+1:end))), atMaxX, mean(mean(Xft(:,P.Tb+1:end))));
   fprintf(P.fid,'W: %8.3f %8.3f %8.3f %8.3f %8.3f\n', min(min(Wft(:,P.Tb+1:end))), atMinW, max(max(Wft(:,P.Tb+1:end))), atMaxW, mean(mean(Wft(:,P.Tb+1:end))));
   fprintf(P.fid,'V: %8.3f %8.3f %8.3f %8.3f %8.3f\n', minVft                     , atMinV, maxVft                     , atMaxV, mean(mean(Vft(:,P.Tb+1:end))));
   fprintf(P.fid,'Exit: %6.4f\n', exiters);
end

% save matrices for next use
if (P.Vsave)
   save ('vfi3Mat', 'V', 'optK', 'optN');
end

fail = 0;

end
