% ----------------------------------------------------
% ANALYZEF:
%   Simulate Panel Data from the Cross-Sectional Distribution of Firms
%   return the simulated moments
%
% By Robert Parham, 2015
% Based on code by Lu Zhang, 2001
%
% ----------------------------------------------------
function [mom, Data] = analyzeF(P, k, n, z, x, w, V, optK, optN, kpd, npd, zpd, xpd, wpd)

% reinitialize time periods
Ts = P.Ts2 + 1;

% prepare zeroed out moms and averaging constant
mom   = zeros(P.Nmoms,1);
reps1 = 1/double(P.Reps);

for rep=1:P.Reps
   % generate panel
   [Kft, Nft, Zft, Xft, Wft, Vft, Gft, ~, ~, ~, ~, ~, sfail] = simCC_mex(P, k, n, z, x, w, V, optK, optN, kpd, npd, zpd, xpd, wpd, Ts);
   if (sfail)
      error('Simulation failed where it should not!\n');
   end

   % create derived panels
   % DV will have NaN in the cells where there's exit
   [GP, IN, OI, IK, FCF, DV] = panCC1(P, Kft, Nft, Zft, Xft, Wft, Vft, Gft);

   % Use derived panels to calculate moments
   [momt, Data] = momCC(P, Kft, Nft, Zft, Xft, Wft, Vft, Gft, GP, IN, OI, IK, FCF, DV);
   mom  = mom + reps1*momt;
end
end