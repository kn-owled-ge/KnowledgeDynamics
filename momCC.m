function [mom, Data] = momCC(P, Kft, Nft, Zft, Xft, Wft, Vft, Gft, GP, IN, OI, IK, FCF, DV)

% Tidy incoming vars
GP  = GP';
IN  = IN';
OI  = OI';
IK  = IK';
FCF = FCF';
DV  = DV';
K  = Kft(:,1:P.Ts2)';
N  = Nft(:,1:P.Ts2)';
Z  = Zft(:,1:P.Ts2)';
X  = Xft(:,1:P.Ts2)';
W  = Wft(:,1:P.Ts2)';
V  = Vft(:,1:P.Ts2)';
G  = Gft(:,1:P.Ts2)';

% scale firm state
K  = exp(K);
N  = exp(N);
V  = exp(V);

% Generate firm age panels
Tvar = ((1:size(V,1))')*ones(1,size(V,2));
Jvar = ones(size(Tvar));
Avar = ones(size(G));
Ivar = ones(size(G));
Ivar(1,:) = (1:P.Ns).*P.Ts2;
O    = ones(1,size(G,2));
for i=2:size(G,1)
   Avar(i,:) = Avar(i-1,:)+O;
   Avar(i,G(i-1,:)~=0) = O(1,G(i-1,:)~=0);
   Ivar(i,:) = Ivar(i-1,:);
   Ivar(i,G(i-1,:)~=0) = Ivar(i,G(i-1,:)~=0) + O(1,G(i-1,:)~=0);
end

% generate data matrix
Data = [vec1(Ivar) vec1(Tvar) vec1(Jvar) vec1(Avar) vec2(Z) vec2(X) vec2(W) vec3(N) vec3(K) vec3(V) vec3(GP) vec3(IN) vec3(OI) vec3(IK) vec3(FCF) vec3(DV)];

% remove all exits (marked by NaN in DV)
mask = isnan(sum(Data(:,end-2:end),2));
Data = Data(~mask,:);

% get moms
mom = distmoms_all(P, Data);
end

function ret = vec1(val)
   r   = reshape(val(2:end-1,:),[],1);
   ret = [r];
end

function ret = vec2(val)
   r   = reshape(val(2:end-1,:),[],1);
   L_r = reshape(val(1:end-2,:),[],1);
   ret = [r L_r];
end

function ret = vec3(val)
   F_r = reshape(val(3:end  ,:),[],1);
   r   = reshape(val(2:end-1,:),[],1);
   L_r = reshape(val(1:end-2,:),[],1);
   ret = [F_r r L_r];
end