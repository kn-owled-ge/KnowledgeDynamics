% ----------------------------------------------------
% MODELF:
%   Run a model and return all moments
%
% over  - a row vector of overriding parameters to calibF. Only numerical
% params that may be adjusted during SMM are allowed in over
%
% By Robert Parham, 2015
%
% ----------------------------------------------------
function [Smom, Sdata, Dmom] = modelF(over, fid, Ddata)

% deal with missing fid
if (nargin<2)
   fid = NaN;
end

% set the random number generator for the estimation
rng(1337);

% Calibrate model and verify parametrization makes sense
[P, Qz, Qx, Qw, k, kp, n, np, z, x, w, fail] = calibF(over,fid);

% Solve model
if(~fail)
   [V, optK, optN, kpd, npd, zpd, xpd, wpd, fail] = mainF(P, Qz, Qx, Qw, k, kp, n, np, z, x, w);
end

% Test for grid failure while solving model
if (fail)
   if (isnan(fid))
      fclose(P.fid);
   end
   Smom  = nan(P.Nmoms,1);
   Dmom  = nan(P.Nmoms,1);
   Sdata = NaN;
   return;
end

% Create model moments
[Smom, Sdata] = analyzeF(P, k, n, z, x, w, V, optK, optN, kpd, npd, zpd, xpd, wpd);

% Create data moments
if nargin>2
   Dmom = distmoms(P, Ddata);
else
   Dmom = nan(P.Nmoms,1);
end

% Print data moments
if (P.DMOMprn)
   fprintf(P.fid,'Data moments:\n');
   for i=1:numel(Dmom)
      fprintf(P.fid,'%8.4f\n', Dmom(i));
   end
end

% Print simulated moments
if(P.SMOMprn)
   fprintf(P.fid,'Simulated moments:\n');
   for i=1:numel(Smom)
      fprintf(P.fid,'%8.4f\n', Smom(i));
   end
end

% Print relative error
if(P.EMOMprn)
   fprintf(P.fid,'Relative error in moments:\n');
   for i=1:numel(Smom)
      fprintf(P.fid,'%8.4f\n', (Smom(i)-Dmom(i))./abs(Dmom(i)));
   end
end

% Close write file if one was opened
if (isnan(fid))
   fclose(P.fid);
end

% make random numbers random again
rng('shuffle');

end
