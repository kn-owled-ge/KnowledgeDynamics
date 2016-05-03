% ----------------------------------------------------
% SIMCC:
%   Simulate the economy
%
% By Robert Parham, 2015
% Based on code by Lu Zhang, 2001
%
% ----------------------------------------------------
function [K, N, Z, X, W, V, G, kd, nd, zd, xd, wd, sfail] = simCC(P, k, n, z, x, w, V0t, optKt, optNt, kd, nd, zd, xd, wd, Ts)

% initialize
sfail = false;
nk    = numel(k);
nn    = numel(n);
nz    = numel(z);
nx    = numel(x);
nw    = numel(w);
kpd   = kd;
npd   = nd;
zpd   = zd;
xpd   = xd;
wpd   = wd;
K     = zeros(P.Ns,Ts);
N     = zeros(P.Ns,Ts);
Z     = zeros(P.Ns,Ts);
X     = zeros(P.Ns,Ts);
W     = zeros(P.Ns,Ts);
V     = zeros(P.Ns,Ts);
G     = zeros(P.Ns,Ts);
O     = ones(P.Ns,1);

% reshape to facilitate interpolation later
optK = reshape(optKt, [nk nn nw nx nz]);
optN = reshape(optNt, [nk nn nw nx nz]);
V0   = reshape(V0t,   [nk nn nw nx nz]);

% start iterations on time
for t = 1:Ts
   % For each firm, interrogate policy and value guides (5D interpolation)
   for j = 1:P.Ns
      % interpolate data (5D interpolation)
      [kpd(j),npd(j),V(j,t)] = m_interp5(k,kd(j),n,nd(j),w,wd(j),x,xd(j),z,zd(j),optK,optN,V0);

      % impose upper bound
      
   end
   
   % do exit and entry
   if (P.ExEnt)
      % define exits
      exent1   = V(:,t)<=P.Vexit+1E-3;
      exent2   = kd    <=P.Kexit+1E-3;
      pexit    = rand(P.Ns,1);
      exent3   = (pexit <= P.Pext);
      %lVK      = V(:,t)-K(:,t); %note both already in log terms
      %Pexta = 4;     % parameter of sigmoid exit function (from data)
      %Pextb = 1.44;  % parameter of sigmoid exit function (from data)
      %Pextc = 0.25;  % probability of exit per period
      %exent3   = (pexit <=(Pextc.*(1./(1+exp(Pexta.*lVK+Pextb))) + P.Pext));
      exent    = exent1 | exent2 | exent3;
      bankrupt = exent & kd>P.Bexit;
      G(:,t)   = exent1*1 + exent2*10 + exent3*100 + bankrupt*1000;
      
      % update panel, with new stochastic states from current distribution
      % make sure distribution of Q upon entry fits observed one
      tocpy      = randi(P.Ns,P.Ns,1);
      toQ        = P.MQent + P.SQent.*randn(P.Ns,1);
      toQ        = toQ(exent);
      zdrnd      = zd(tocpy);
      %kdrnd      = kd(tocpy);
      %cktmp      = kdrnd(exent);
      ztmp       = zdrnd(exent);
      ktmp       = ones(size(ztmp));
      ntmp       = ones(size(ztmp));
      
      for i1 = 1:numel(ztmp)
         currdist = 100;
         for i2 = 1:nk
            uk = k(i2);
            %if (uk>cktmp(i1))
            %   continue;
            %end
            un = n;
            [~,~,uv] = m_interp5(k,uk,n,un,w,wd(1),x,xd(1),z,ztmp(i1),optK,optN,V0);
            newdist = abs(uv-uk - toQ(i1));
            if (newdist<currdist)
               ktmp(i1) = uk;
               ntmp(i1) = un;
               currdist = newdist;
            end
         end
      end
      
      zpd(exent) = ztmp;
      kpd(exent) = ktmp;
      npd(exent) = ntmp;
   else
      zpd = zd;
      xpd = xd;
      wpd = wd;
   end
   
   % save current period states
   K(:, t) = kd;
   N(:, t) = nd;
   Z(:, t) = zd;
   X(:, t) = xd;
   W(:, t) = wd;
   
   % update next period capital stocks, and make sure on grid
   kd = min(max(kpd,k(1)),k(nk));
   nd = min(max(npd,n(1)),n(nn));
   
   % update idiosyncratic shocks in continuous state space
   shckz = randn(P.Ns, 1);
   %shckx = randn(P.Ns, 1);
   %shckw = randn(P.Ns, 1);
   zd    = P.rhoZ.*zpd + (1-P.rhoZ).*P.muZ + P.sigZ.*shckz;
   %xd    = P.rhoX.*xpd + (1-P.rhoX).*P.muX + P.sigX.*shckx;
   %wd    = P.rhoW.*wpd + (1-P.rhoW).*P.muW + P.sigW.*shckw;
   xd = xpd;
   wd = wpd;
   
   % visuals
   %if (mod(t, P.SIMprn) == 0)
   %   fprintf(P.fid,'Periods: %8d\n', t);
   %end
end

end


%-------------------------------------------------------------------------------------------------------
% Linear interpolation routine dealing with the staggared optK,optN,V0 structures
% 
% OUTPUT:
%   kp - next period k
%   np - next period n
%   v  - current period v
% 
% INPUT: 
%   k,kd,n,nd,x,xd,w,wd,z,zd  - grids and locations within them
%   optK,optN,V0              - matrices to interpolate
% 
%------------------------------------------------------------------------------------------------------
function [kp,np,v] = m_interp5(k,kd,n,nd,w,wd,x,xd,z,zd,optK,optN,V0)

   kpos = sum(k <= kd);
   npos = sum(n <= nd);
   wpos = sum(w <= wd);
   xpos = sum(x <= xd);
   zpos = sum(z <= zd);

   nk   = numel(k);
   nn   = numel(n);
   nz   = numel(z);
   nx   = numel(x);
   nw   = numel(w);

   pos  = [kpos npos wpos xpos zpos];
   down = max([pos ; ones(size(pos))]);
   up   = min([pos+1 ; nk nn nw nx nz]);
   both = [down ; up];

   % Get outer edges of hypercube we're in
   Kpt = zeros(2,2,2,2,2);
   Npt = zeros(2,2,2,2,2);
   Vpt = zeros(2,2,2,2,2);
   for i1 = 1:2
      for i2=1:2
         for i3=1:2
            for i4=1:2
               for i5=1:2
                  Kpt(i1,i2,i3,i4,i5) = optK(both(i1,1),both(i2,2),both(i3,3),both(i4,4),both(i5,5));
                  Npt(i1,i2,i3,i4,i5) = optN(both(i1,1),both(i2,2),both(i3,3),both(i4,4),both(i5,5));
                  Vpt(i1,i2,i3,i4,i5) = V0  (both(i1,1),both(i2,2),both(i3,3),both(i4,4),both(i5,5));
               end
            end
         end
      end
   end

   % get interpolation weights per axis
   kwgt = min(max((kd - k(down(1)))/(k(up(1))-k(down(1))),0),1);
   nwgt = min(max((nd - n(down(2)))/(n(up(2))-n(down(2))),0),1);
   wwgt = min(max((wd - w(down(3)))/(w(up(3))-w(down(3))),0),1);
   xwgt = min(max((xd - x(down(4)))/(x(up(4))-x(down(4))),0),1);
   zwgt = min(max((zd - z(down(5)))/(z(up(5))-z(down(5))),0),1);

   % interpolate successively
   Kpt = squeeze(kwgt*Kpt(2, :, :, :, :) + (1 - kwgt)*Kpt(1, :, :, :, :));
   Npt = squeeze(kwgt*Npt(2, :, :, :, :) + (1 - kwgt)*Npt(1, :, :, :, :));
   Vpt = squeeze(kwgt*Vpt(2, :, :, :, :) + (1 - kwgt)*Vpt(1, :, :, :, :));
   Kpt = squeeze(nwgt*Kpt(2, :, :, :)    + (1 - nwgt)*Kpt(1, :, :, :));
   Npt = squeeze(nwgt*Npt(2, :, :, :)    + (1 - nwgt)*Npt(1, :, :, :));
   Vpt = squeeze(nwgt*Vpt(2, :, :, :)    + (1 - nwgt)*Vpt(1, :, :, :));
   Kpt = squeeze(wwgt*Kpt(2, :, :)       + (1 - wwgt)*Kpt(1, :, :));
   Npt = squeeze(wwgt*Npt(2, :, :)       + (1 - wwgt)*Npt(1, :, :));
   Vpt = squeeze(wwgt*Vpt(2, :, :)       + (1 - wwgt)*Vpt(1, :, :));
   Kpt = squeeze(xwgt*Kpt(2, :)          + (1 - xwgt)*Kpt(1, :));
   Npt = squeeze(xwgt*Npt(2, :)          + (1 - xwgt)*Npt(1, :));
   Vpt = squeeze(xwgt*Vpt(2, :)          + (1 - xwgt)*Vpt(1, :));

   % set next state and value funtion
   kp = zwgt*Kpt(2) + (1 - zwgt)*Kpt(1);
   np = zwgt*Npt(2) + (1 - zwgt)*Npt(1);
   v  = zwgt*Vpt(2) + (1 - zwgt)*Vpt(1);
   
end


