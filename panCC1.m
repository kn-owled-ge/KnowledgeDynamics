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

IN = zeros(size(kmtrx));
GP = zeros(size(kmtrx));

% physical investment
IK = ((max(kpmtrx./kmtrx - P.lamb1k,1e-3)./P.lamb2k).^(1/P.lamb3k)).*kmtrx;


% NaN IK where there was exit
IK(Gft(:,1:end-1)~=0)=NaN;

% cash flow
OI   = zmtrx.*kmtrx.^P.theta;
FCF  = (1-P.tauC).*(OI-IK);
DV   = FCF;

end
