function mom = distmoms(P, Data)
% distmoms - recieves Data = [gvkey fyear ind age Z L_Z X L_X W L_W F_N N L_N F_K K L_K F_V V L_V F_GP GP L_GP F_IN IN L_IN F_OI OI L_OI F_IK IK L_IK F_CF CF L_CF F_D D L_D]
% returns moments vector

warning('off','stats:robustfit:RankDeficient');

% prepare vars
mom  = zeros(P.Nmoms,1);
pos  = 0;

%fyear = Data(:,2);
ind   = Data(:,3);
%age   = Data(:,4);
%Z     = Data(:,5);
%L_Z   = Data(:,6);
%X     = Data(:,7);
%L_X   = Data(:,8);
%W     = Data(:,9);
%L_W   = Data(:,10);
%F_N   = Data(:,11);
%N     = Data(:,12);
%L_N   = Data(:,13);
F_K   = Data(:,14);
K     = Data(:,15);
L_K   = Data(:,16);
F_V   = Data(:,17);
V     = Data(:,18);
L_V   = Data(:,19);
F_GP  = Data(:,20);
GP    = Data(:,21);
%L_GP  = Data(:,22);
F_IN  = Data(:,23);
IN    = Data(:,24);
%L_IN  = Data(:,25);
%F_OI  = Data(:,26);
OI    = Data(:,27);
%L_OI  = Data(:,28);
F_IK  = Data(:,29);
IK    = Data(:,30);
%L_IK  = Data(:,31);
%F_CF  = Data(:,32);
%CF    = Data(:,33);
%L_CF  = Data(:,34);
%F_DV  = Data(:,35);
DV    = Data(:,36);
%L_DV  = Data(:,37);

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
   lGP   = real(log(GP));
   lF_GP = real(log(F_GP));
   lIN   = real(log(IN));
   lF_IN = real(log(F_IN));
   lIK   = real(log(IK));
   lF_IK = real(log(F_IK));

   F_Q = lF_V  - lF_K;
   Q   = lV    - lK;
   F_R = lF_GP - lF_K;
   R   = lGP   - lK;
   F_S = lF_IN - lF_K;
   S   = lIN   - lK;
   F_T = lF_IK - lF_K;
   T   = lIK   - lK;
   F_U = lF_IN - lF_V;
   U   = lIN   - lV;

   OIK   = OI./K;
   OIV   = OI./V;
   DVK   = DV./K;
   DVV   = DV./V;

   % pre-processing
   lb  = 1e-2;
   GV  = lF_K>lb & lK>lb & lL_K>lb & lF_V>lb & lV>lb & lL_V>lb;
   FE  = ind(GV); 
   uFE = unique(FE)';
   rFE = 1*(FE*ones(size(uFE)) == ones(size(FE))*uFE);
   
   
   %%% logK %%%
   VR  = lK;
   FVR = lF_K;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   FVR1  = FVR0(good1);
   prc2  = prctile(FVR1,[0.1,99.9]);
   good2 = FVR1>prc2(1) & FVR1<prc2(2);

   [b1,~] = robustfit(VR1(good2),FVR1(good2));

   mom(pos+1) = median(VR1);
   mom(pos+2) = mad(VR1);
   mom(pos+3) = b1(2);
   pos = pos+3;
   
   
   %%% logV %%%
   VR  = lV;
   FVR = lF_V;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   FVR1  = FVR0(good1);
   prc2  = prctile(FVR1,[0.1,99.9]);
   good2 = FVR1>prc2(1) & FVR1<prc2(2);

   [b1,~] = robustfit(VR1(good2),FVR1(good2));

   mom(pos+1) = median(VR1);
   mom(pos+2) = mad(VR1);
   mom(pos+3) = b1(2);
   pos = pos+3;


   %%% Q %%%
   VR  = Q;
   FVR = F_Q;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   FVR1  = FVR0(good1);
   prc2  = prctile(FVR1,[0.1,99.9]);
   good2 = FVR1>prc2(1) & FVR1<prc2(2);
   
   [b1,~] = robustfit(VR1(good2),FVR1(good2));

   mom(pos+1) = median(VR1);
   mom(pos+2) = mad(VR1);
   mom(pos+3) = b1(2);
   pos = pos+3;
   
   % save Q results for latet
   uQ = VR0;
   gQ = good1;

   
   %%% R %%%
   VR  = R;
   FVR = F_R;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   FVR1  = FVR0(good1);
   prc2  = prctile(FVR1,[0.1,99.9]);
   good2 = FVR1>prc2(1) & FVR1<prc2(2);

   [b1,~] = robustfit(VR1(good2),FVR1(good2));

   mom(pos+1) = median(VR1);
   mom(pos+2) = mad(VR1);
   mom(pos+3) = b1(2);
   pos = pos+3;

   % save R results for latet
   uR = VR0;
   gR = good1;

   
   %%% S %%%
   VR  = S;
   %FVR = F_S;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      %FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      %FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   %FVR1  = FVR0(good1);
   %prc2  = prctile(FVR1,[0.1,99.9]);
   %good2 = FVR1>prc2(1) & FVR1<prc2(2);

   mom(pos+1) = median(VR1);
   mom(pos+2) = mad(VR1);
   pos = pos+2;

   % save S results for latet
   uS = VR0;
   gS = good1;


   %%% T %%%
   VR  = T;
   %FVR = F_T;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      %FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      %FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   %FVR1  = FVR0(good1);
   %prc2  = prctile(FVR1,[0.1,99.9]);
   %good2 = FVR1>prc2(1) & FVR1<prc2(2);

   mom(pos+1) = median(VR1);
   mom(pos+2) = mad(VR1);
   pos = pos+2;
   
   % save T results for later
   uT = VR0;
   gT = good1;


   %%% OIK %%%
   VR  = OIK;
   %FVR = F_OIK;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      %FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      %FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   %FVR1  = FVR0(good1);

   mom(pos+1) = median(VR1);
   pos = pos+1;
   
   
   %%% OIV %%%
   VR  = OIV;
   %FVR = F_OIV;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      %FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      %FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   %FVR1  = FVR0(good1);

   mom(pos+1) = median(VR1);
   pos = pos+1;
   
   
   %%% DVK %%%
   VR  = DVK;
   %FVR = F_DVK;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      %FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      %FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   %FVR1  = FVR0(good1);

   mom(pos+1) = median(VR1);
   pos = pos+1;
   
   
   %%% DVV %%%
   VR  = DVV;
   %FVR = F_DVV;

   % take ind fixed effects
   if(size(uFE,2)>1)
      bVR  = robustfit(rFE,VR(GV),'bisquare',4.685,'off');
      mVR  = median(VR(GV));
      VR0  = VR(GV) - rFE*bVR + mVR;
      %FVR0 = FVR(GV) - rFE*bVR + mVR;
   else
      VR0  = VR(GV);
      %FVR0 = FVR(GV);
   end

   % cut off outliers
   prc1  = prctile(VR0,[0.5,99.5]);
   good1 = VR0>prc1(1) & VR0<prc1(2);
   VR1   = VR0(good1);
   %FVR1  = FVR0(good1);

   mom(pos+1) = median(VR1);
   pos = pos+1;
   
%{   
   %%% reg [Q,R],S %%%
   [b,~] = robustfit([uQ(gQ&gR) uR(gQ&gR)], uS(gQ&gR));
   mom(pos+1:pos+2) = b(2:3);
   pos = pos+2;

   
   %%% reg [Q,R],T %%%
   [b,~] = robustfit([uQ(gQ&gR) uR(gQ&gR)], uT(gQ&gR));
   mom(pos+1:pos+2) = b(2:3);
   pos = pos+2;
%}

catch ex
   disp(ex);
   disp(mom);
   mom  = nan(P.Nmoms,1);
end   

end


