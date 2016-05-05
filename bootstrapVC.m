function VC = bootstrapVC(P, Data, VCname)

nobs   = size(Data,1);
reps   = P.VCreps;
inf    = zeros(reps,P.Nmoms);

rng(1337);
warning('off','all');
parfor i=1:reps
   Dati     = Data(randi(nobs,nobs,1),:);
   mom      = distmoms_all(P,Dati);
   inf(i,:) = mom';
   if (mod(i,10)==1)
      disp(i);
   end
end
warning('on','all');

meaninf = mean(inf,1);
inf1    = inf - ones(reps,1)*meaninf;
VC      = inf1'*inf;

save (VCname, 'VC');

end