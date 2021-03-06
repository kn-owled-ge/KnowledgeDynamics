function mom = distmoms_all(P, Data)
% distmoms - recieves Data = [gvkey fyear ind age Z L_Z X L_X W L_W F_N N L_N F_K K L_K F_V V L_V F_GP GP L_GP F_IN IN L_IN F_OI OI L_OI F_IK IK L_IK F_CF CF L_CF F_D D L_D]
% returns moments vector

warning('off','stats:robustfit:RankDeficient');

% prepare vars
mom  = zeros(P.Nmoms,1);
pos  = 0;

fyear = Data(:,2);
ind   = Data(:,3);
age   = Data(:,4);
Z     = Data(:,5);
L_Z   = Data(:,6);
X     = Data(:,7);
L_X   = Data(:,8);
W     = Data(:,9);
L_W   = Data(:,10);
F_N   = Data(:,11);
N     = Data(:,12);
L_N   = Data(:,13);
F_K   = Data(:,14);
K     = Data(:,15);
L_K   = Data(:,16);
F_V   = Data(:,17);
V     = Data(:,18);
L_V   = Data(:,19);
F_GP  = Data(:,20);
GP    = Data(:,21);
L_GP  = Data(:,22);
F_IN  = Data(:,23);
IN    = Data(:,24);
L_IN  = Data(:,25);
F_OI  = Data(:,26);
OI    = Data(:,27);
L_OI  = Data(:,28);
F_IK  = Data(:,29);
IK    = Data(:,30);
L_IK  = Data(:,31);
F_CF  = Data(:,32);
CF    = Data(:,33);
L_CF  = Data(:,34);
F_DV  = Data(:,35);
DV    = Data(:,36);
L_DV  = Data(:,37);

%{
% Verify firm size is within reasonable bounds
trlK = 3;
truK = 8;
lowK  = sum(log(K)<trlK)/Ndt;
highK = sum(log(K)>truK)/Ndt;
if (lowK > 0.6 || highK > 0.6)
   mom  = nan(P.Nmoms,1);
   return;
end
%}

try
   % prepare variables
   lL_K  = real(log(L_K));
   lK    = real(log(K));
   lF_K  = real(log(F_K));
   lL_V  = real(log(L_V));
   lV    = real(log(V));
   lF_V  = real(log(F_V));
   lL_GP = real(log(L_GP));
   lGP   = real(log(GP));
   lF_GP = real(log(F_GP));
   lL_IN = real(log(L_IN));
   lIN   = real(log(IN));
   lF_IN = real(log(F_IN));
   lL_IK = real(log(L_IK));
   lIK   = real(log(IK));
   lF_IK = real(log(F_IK));

   F_Q = lF_V  - lF_K;
   Q   = lV    - lK;
   L_Q = lL_V  - lL_K;
   F_R = lF_GP - lF_K;
   R   = lGP   - lK;
   L_R = lL_GP - lL_K;
   F_S = lF_IN - lF_K;
   S   = lIN   - lK;
   L_S = lL_IN - lL_K;
   F_T = lF_IK - lF_K;
   T   = lIK   - lK;
   L_T = lL_IK - lL_K;
   F_U = lF_IN - lF_V;
   U   = lIN   - lV;
   L_U = lL_IN - lL_V;

   F_OIK = F_OI./F_K;
   OIK   = OI./K;
   L_OIK = L_OI./L_K;
   F_OIV = F_OI./F_V;
   OIV   = OI./V;
   L_OIV = L_OI./L_V;
   F_DVK = F_DV./F_K;
   DVK   = DV./K;
   L_DVK = L_DV./L_K;
   F_DVV = F_DV./F_V;
   DVV   = DV./V;
   L_DVV = L_DV./L_V;

   rF_K  = lF_K  - lK;
   rK    = lK    - lL_K;
   rF_V  = lF_V  - lV;
   rV    = lV    - lL_V;
   rF_GP = lF_GP - lGP;
   rGP   = lGP   - lL_GP;
   rF_IN = lF_IN - lIN;
   rIN   = lIN   - lL_IN;
   rF_IK = lF_IK - lIK;
   rIK   = lIK   - lL_IK;
   
   rF_Q = F_Q - Q;
   rQ   = Q   - L_Q;
   rF_R = F_R - R;
   rR   = R   - L_R;
   rF_S = F_S - S;
   rS   = S   - L_S;
   rF_T = F_T - T;
   rT   = T   - L_T;
   rF_U = F_U - U;
   rU   = U   - L_U;
   
   rF_OIK = F_OIK - OIK;
   rOIK   = OIK   - L_OIK;
   rF_OIV = F_OIV - OIV;
   rOIV   = OIV   - L_OIV;
   rF_DVK = F_DVK - DVK;
   rDVK   = DVK   - L_DVK;
   rF_DVV = F_DVV - DVV;
   rDVV   = DVV   - L_DVV;
   
   % pre-processing
   lb  = 1e-2;
   goodVar = lF_K>lb & lK>lb & lL_K>lb & lF_V>lb & lV>lb & lL_V>lb;
   %doGraph = [1 1 1 0];
   doGraph = zeros(1,4);
   doAll0  = [1 1 0 1 1 1 1 0 0];
   doAll1  = [1 1 0 1 1 1 1 1 1];
   doAll2  = [1 1 1 1 1 1 1 1 1];
   
   % features: lK, lV, lGP, lIN, lIK
   [mom1, lK1] = features(lF_K,lK,goodVar,nan,ind,'normal','logK',doGraph,doAll0);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   
   [mom1,lV1] = features(lF_V,lV,goodVar,lK1,ind,'normal','logV',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,lGP1] = features(lF_GP,lGP,goodVar,lK1,ind,'normal','logGP',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,lIN1] = features(lF_IN,lIN,goodVar,lK1,ind,'normal','logIN',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,lIK1] = features(lF_IK,lIK,goodVar,lK1,ind,'normal','logIK',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;


   % features: Q, R, S, T, U

   [mom1,Q1] = features(F_Q,Q,goodVar,lK1,ind,'normal','Q',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,R1] = features(F_R,R,goodVar,lK1,ind,'normal','R',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,S1] = features(F_S,S,goodVar,lK1,ind,'normal','S',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,T1] = features(F_T,T,goodVar,lK1,ind,'normal','T',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,U1] = features(F_U,U,goodVar,lK1,ind,'normal','U',doGraph,doAll1);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   
   % features: OIK, OIV, DVK, DVV

   [mom1,OIK1] = features(F_OIK,OIK,goodVar,lK1,ind,'tdist','OIK',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,OIV1] = features(F_OIV,OIV,goodVar,lK1,ind,'tdist','OIV',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,DVK1] = features(F_DVK,DVK,goodVar,lK1,ind,'tdist','DVK',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,DVV1] = features(F_DVV,DVV,goodVar,lK1,ind,'tdist','DVV',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   
   % features: rK, rV, rGP, rIN, rIK
   
   [mom1,rK1] = features(rF_K,rK,goodVar,lK1,ind,'tdist','rK',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,rV1] = features(rF_V,rV,goodVar,lK1,ind,'tdist','rV',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,rGP1] = features(rF_GP,rGP,goodVar,lK1,ind,'tdist','rGP',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,rIN1] = features(rF_IN,rIN,goodVar,lK1,ind,'tdist','rIN',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,rIK1] = features(rF_IK,rIK,goodVar,lK1,ind,'tdist','rIK',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   
   
   % features: rQ, rR, rS, rT, rU
   
   [mom1,rQ1] = features(rF_Q,rQ,goodVar,lK1,ind,'tdist','rQ',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   
   [mom1,rR1] = features(rF_R,rR,goodVar,lK1,ind,'tdist','rR',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   
   [mom1,rS1] = features(rF_S,rS,goodVar,lK1,ind,'tdist','rS',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   
   [mom1,rT1] = features(rF_T,rT,goodVar,lK1,ind,'tdist','rT',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   
   [mom1,rU1] = features(rF_U,rU,goodVar,lK1,ind,'tdist','rU',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   

   % features: rOIK, rOIV, rDVK, rDVV

   [mom1,rOIK1] = features(rF_OIK,rOIK,goodVar,lK1,ind,'tdist','rOIK',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   
   [mom1,rOIV1] = features(rF_OIV,rOIV,goodVar,lK1,ind,'tdist','rOIV',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   [mom1,rDVK1] = features(rF_DVK,rDVK,goodVar,lK1,ind,'tdist','rDVK',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;
   
   [mom1,rDVV1] = features(rF_DVV,rDVV,goodVar,lK1,ind,'tdist','rDVV',doGraph,doAll2);
   pos1 = size(mom1,1);
   mom(pos+(1:pos1)) = mom1;
   pos  = pos + pos1;

   
   % regression coeffs
   
   [b, st] = robustfit([lK1,lGP1,lIN1,lIK1],lV1);
   mom(pos+(1:4)) = b(2:5);
   mom(pos+5) = st.robust_s;
   pos = pos + 5;
   
   [b, st] = robustfit([Q1,R1],S1);
   mom(pos+(1:2)) = b(2:3);
   mom(pos+3) = st.robust_s;
   pos = pos + 3;

   [b, st] = robustfit([Q1,R1],T1);
   mom(pos+(1:2)) = b(2:3);
   mom(pos+3) = st.robust_s;
   pos = pos + 3;
   
   
   % correlations
   Dat = [lK1, lV1, lGP1, lIN1, lIK1, Q1, R1, S1, T1, U1, OIK1, OIV1, DVK1, DVV1, rK1, rV1, rGP1, rIN1, rIK1, rQ1, rR1, rS1, rT1, rU1, rOIK1, rOIV1, rDVK1, rDVV1];
   cr  = corr(Dat);
   crv = cr(tril(true(size(cr)),-1));
   mom(pos+(1:size(crv,1))) = crv;
   
catch ex
   disp(ex);
   disp(mom);
   mom  = nan(P.Nmoms,1);
end   

end


