% ----------------------------------------------------
% PANCC1:
%   Create derived panels from simulation panels.
%
% By Robert Parham, 2015
% Based on code by Lu Zhang, 2001
%
% ----------------------------------------------------
function [GP, IN, OI, IK, FCF, DV] = panCC1(P, Kft, Nft, Zft, Xft, Wft, Vft, Gft)

% move from scale to level
eKft = exp(Kft);
eZft = exp(Zft);

% define matrices
kmtrx  = eKft(:,1:end-1);
kpmtrx = eKft(:,2:end);
zmtrx  = eZft(:,1:end-1);
mInf   = 1E30.*ones(size(kmtrx));

IN = zeros(size(kmtrx));
GP = zeros(size(kmtrx));

% physical investment
IK = ((max(kpmtrx./kmtrx - (1-P.delK/P.lambK),1e-3).*P.lambK./(P.delK^(1-P.lambK))).^(1/P.lambK)).*kmtrx;


% NaN IK where there was exit
IK(Gft(:,1:end-1)~=0)=NaN;

% cash flow
OI   = zmtrx.*kmtrx.^P.theta;
FCF  = (1-P.tauC).*OI + P.tauC.*P.delK.*kmtrx - IK;

if (P.doFin)
   FUB     = -(1/(2*P.gama)).*kmtrx;
   FCFpos  = FCF>=0;
   FCFneg  = ~FCFpos;
   FCFlow  = FCF<FUB;
   FCF1    = -(1/P.gama).*kmtrx + sqrt(max((kmtrx./P.gama).^2 + (2*kmtrx.*FCF./P.gama),0));
   if(P.FinEx)
      DV = FCFpos.*FCF + FCFneg.*FCF1 - FCFlow.*mInf;
   else
      DV = FCFpos.*FCF + FCFneg.*FCF1;
   end
else
   DV = FCF;
end

end
