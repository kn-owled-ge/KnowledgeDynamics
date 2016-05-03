% ----------------------------------------------------
% VFICC:
%   For given parameters solve for the value function and optimal decision rules
%
% By Robert Parham, 2015
% Based on code by Lu Zhang, 2001
%
% ----------------------------------------------------
function [V, optK, optN, vfail] = vfiCC(P, Qx, Qz, Qw, k, kp, n, np, z, x, w, V0)

% single-ize vars, for less memory use
Qx  = single(Qx);
Qz  = single(Qz);
Qw  = single(Qw);
lk  = single(k);
k   = exp(lk);
lkp = single(kp);
kp  = exp(lkp);
ln  = single(n);
n   = exp(ln);
lnp = single(np);
np  = exp(lnp);
z   = single(z);
x   = single(x);
w   = single(w);
V0  = single(exp(V0));

% initialize size variables
nk       = numel(k);
nkp      = numel(kp);
nn       = numel(n);
nnp      = numel(np);
nz       = numel(z);
nx       = numel(x);
nw       = numel(w);

% calculate compound sizes
nkn      = nk*nn;
nknp     = nkp*nnp;
nzxw     = nz*nx*nw;

% compound the (z x w) space: w changes the fastes, followed by x, and z changes the slowest
Qxw  = single(kron(Qx, Qw ));
Qzxw = single(kron(Qz, Qxw));

% State-space matrices
lkmtrx = repmat(lk, P.nn, P.nz*P.nx*P.nw);                               % k changes the fastest
lnmtrx = kron(repmat(ln, 1, P.nz*P.nx*P.nw), ones(P.nk,1));              % n changes the slowest
%wmtrx  = repmat(repmat(w, P.nz*P.nx, 1)', P.nk*P.nn, 1);                 % w changes the fastest
%xmtrx  = repmat(repmat(kron(x, ones(P.nw, 1)), P.nz, 1)', P.nk*P.nn, 1); % x changes the second fastest
zmtrx  = repmat(kron(z, ones(P.nx*P.nw, 1))', P.nk*P.nn, 1);             % z changes the slowest

% gross profits (independent of firm choices)
GP = (1-P.tauC).*exp(zmtrx + P.theta.*lkmtrx);

% Initializing VFI
errK     = single(Inf);
errN     = single(Inf);
errV     = single(Inf);
iter     = 0;
vfail    = 0;
V        = V0;
Vold     = V0;
optK     = single(zeros(nkn, nzxw));
optN     = single(zeros(nkn, nzxw));
optKold  = optK + 1;
optNold  = optN + 1;
Obj      = single(zeros(nkn, nzxw));

% Start VFI
while (1)
   if ((errV < P.Verr && errK < P.Kerr && errN < P.Nerr) || iter>=P.VIter)
      if (iter>=P.VIter)
         vfail=1;
      end
      V = log(V);
      break
   else
      iter = iter + 1;
      for jz = 1:nz
         for jx = 1:nx
            for jw = 1:nw
               % Expectation vector position given current z,x,w
               pos = (jz-1)*nx*nw + (jx-1)*nw + jw;
               
               % next period value function by linear interpolation: nknp by 1 matrix
               % note interpolation is linear in logs
               EMV = exp(m_interp2(lk, ln, log(P.beta.*Vold*Qzxw(:,pos)), lkp, lnp));
               
               % firms exit with some probability (possibly dependent on V/K)
               if (P.ExEnt)
                  %Pexta = 4;     % parameter of sigmoid exit function (from data)
                  %Pextb = 1.44;  % parameter of sigmoid exit function (from data)
                  %Pextc = 0.25;  % probability of exit per period
                  %lVK = log(EMV./kmtrx);
                  %pVK = Pextc.*(1./(1+exp(Pexta.*lVK+Pextb))) + P.Pext;
                  % firms that exit at Pext shock are only to make room for
                  % entrants - they get full continuation value (e.g. going
                  % privte) so no need to account for them in the value 
                  % function 
                  % pVK = P.Pext;
                  % EMV = EMV.*(1-pVK);
               end
               
               for jn = 1:nn
                  IN    = 0;
                  tIN   = kron(IN, ones(nkp,1));
                  
                  for jk = 1:nk
                     IK  = ((max(kp./k(jk) - P.lamb1k,0)./P.lamb2k).^(1/P.lamb3k)).*k(jk);
                     tIK = kron(ones(nnp,1), IK);
                     
                     tObj = -(1-P.tauC).*tIK + EMV;
                     
                     % direct discrete maximization
                     [Obj((jn-1)*nk+jk, pos), optJ]  = max(tObj);
                     
                     %optJ is 1 by 1, but its value need parsing
                     optJN = floor((optJ-1)/nkp)+1;
                     optJK = optJ-(optJN-1)*nkp;
                     optN((jn-1)*nk+jk, pos) = log(np(optJN));
                     optK((jn-1)*nk+jk, pos) = log(kp(optJK));
                  end
               end
            end
         end
      end
      
      % update value function
      V     = max(GP + Obj,exp(-10));

      % convergence criteria
      errK  = max(max(abs(optK   - optKold)));
      errN  = max(max(abs(optN   - optNold)));
      errV  = max(max(abs(log(V) - log(Vold))));

      % revise Initial Guess
      Vold    = V;
      optKold = optK;
      optNold = optN;

      % visuals
      if (mod(iter, P.VFIprn) == 0)
         %fprintf(P.fid,'errV: %10.7f, errK: %10.7f, errN: %10.7f\n', errV, errK, errN);
      end
   end
end

% visuals
%fprintf(P.fid,'errV: %10.7f, errK: %10.7f, errN: %10.7f\n', errV, errK, errN);
   
optK = double(optK);
optN = double(optN);
V    = double(V);
end



%-------------------------------------------------------------------------------------------------------
% Linear interpolation routine dealing with the staggared k,n structure
% 
% OUTPUT:
%   v - function values on non-grid points (n1*n2 by col matrix)
% 
% INPUT: 
%   x1  - grid1 (m1 by one vector) 
%   x2  - grid2 (m2 by one vector)
%   y   - vector function defined on the grid x1 X x2 (m1*m2 by col matrix)
%   u1  - non-grid points on which y(x) is to be interpolated (n1 by one vector)
%   u2  - non-grid points on which y(x) is to be interpolated (n2 by one vector)
% 
%------------------------------------------------------------------------------------------------------
function v = m_interp2(x1, x2, y, u1, u2)

m1    = numel(x1);
m2    = numel(x2);
n1    = numel(u1);
n2    = numel(u2);
col   = size(y,2);

tv      = single(zeros(n1*m2,col));
v       = single(zeros(n1*n2, col));

% First interp1 all the x1 to u1 for each x2, yielding a n1*m2 matrix
for i = 1:m2
   tv(((i-1)*n1+1):i*n1,:) = interp1(x1, y(((i-1)*m1+1):i*m1,:), u1, 'linear', 'extrap');
end

% Next interp1 all the x2 to u2, yielding a n1*n2 matrix to be returned
if (m2>1)
   for i = 1:n1
      subind1 = i:n1:n1*m2;
      subind2 = i:n1:n1*n2;
      v(subind2,:) = interp1(x2, tv(subind1,:), u2, 'linear', 'extrap');
   end
else
   assert(m2==1 & n2==1);
   v = tv;
end

end
