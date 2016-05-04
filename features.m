%
% features of a distribution
%
% for normal   : mu, sig
% for tdist    : mu, sig, nu
% for alaplace : m, lambda, kappa
% for all      : rho, sigeps, rho-, rho+, slope on lK, iqr on lK
function [mom,V0] = features(FV,V,GV,lK,ind,type,name,doGraph,doMoms)

   doSave = false;
   
   % take ind fixed effects
   FE  = ind(GV); 
   uFE = unique(FE)';
   if(size(uFE,2)>1)
      rFE = 1*(FE*ones(size(uFE)) == ones(size(FE))*uFE);
      bV  = robustfit(rFE,V(GV),'bisquare',4.685,'off');
      mV  = median(V(GV));
      V0  = V(GV) - rFE*bV + mV;
      FV0 = FV(GV) - rFE*bV + mV;
   else
      V0  = V(GV);
      FV0 = FV(GV);
   end

   % cut off outliers
   prc1  = prctile(V0,[0.5,99.5]);
   good1 = V0>prc1(1) & V0<prc1(2);
   V1    = V0(good1);
   FV1   = FV0(good1);
   prc2  = prctile(FV1,[0.1,99.9]);
   good2 = FV1>prc2(1) & FV1<prc2(2);
   mom   = zeros(1,9);
   
   % declare which graphs to print
   % 1:stat dist, 2:rho, 3:k slope, 4:iqr slope
   if(isnan(doGraph))
      doGraph=[0 0 0 0];
   end

   % declare which moms to return
   % 1:dist1, 2:dist2, 3:dist3, 4:rho, 5:sigeps, 6:rho-, 7:rho+, 8:scale, 9:iqr
   if(isnan(doMoms))
      doMoms=ones(1,9);
      if (strcmp(type,'normal'))
         doMoms(3) = 0;
      end
      if (sum(isnan(lK))>0)
         doMoms(8) = 0;
         doMoms(9) = 0;
      end
   end

   if (doMoms(1) || doMoms(2) || doMoms(3) || doMoms(6) || doMoms(7))
      % statics of distribution
      switch type
         case 'normal'
            pdf_tn = @(x,mu,sigma) pdf('normal',x,mu,sigma) ./ (cdf('normal',prc1(2),mu,sigma)-cdf('normal',prc1(1),mu,sigma));
            pr     = mle(V1, 'pdf',pdf_tn, 'start',[mean(V1) std(V1)], 'lower',[prc1(1) 0], 'upper',[prc1(2) Inf]);
            mom(1) = pr(1); %mu
            mom(2) = pr(2); %sigma
            if (doGraph(1))
               figure('name',[name ' - static distribution']);
               hold on
               [~,cen]  = hist(V1,100);
               hist(V1,100);
               t = pdf_tn(cen,pr(1),pr(2));
               n = size(V1,1);
               plot(cen,t./sum(t).*n,'-g','LineWidth',2);
               hold off
               loglik = sum(log(pdf_tn(V1,pr(1),pr(2))));
               disp(loglik);
               if (doSave)
                  saveas(gcf,[name '1.jpg']);
               end
            end
         case 'tdist'
            pdf_tt = @(x,mu,sigma,nu) pdf('tLocationScale',x,mu,sigma,nu) ./ (cdf('tLocationScale',prc1(2),mu,sigma,nu)-cdf('tLocationScale',prc1(1),mu,sigma,nu));
            pr     = mle(V1, 'pdf',pdf_tt, 'start',[mean(V1) std(V1) 4], 'lower',[prc1(1) 0 1E-1], 'upper',[prc1(2) Inf Inf]);
            mom(1) = pr(1); %mu
            mom(2) = pr(2); %sigma
            mom(3) = pr(3); %nu
            if (doGraph(1))
               figure('name',[name ' - static distribution']);
               hold on
               [~,cen]  = hist(V1,100);
               hist(V1,100);
               t = pdf_tt(cen,pr(1),pr(2),pr(3));
               n = size(V1,1);
               plot(cen,t./sum(t).*n,'-g','LineWidth',2);
               hold off
               loglik = sum(log(pdf_tt(V1,pr(1),pr(2),pr(3))));
               disp(loglik);
               if (doSave)
                  saveas(gcf,[name '1.jpg']);
               end
            end
         case 'alaplace'
            pdf_al = @(x,m,lam,kap) ((lam/(kap+(1/kap))).*exp(-(x-m).*lam.*sign(x-m).*kap.^sign(x-m)));
            cdf_al = @(x,m,lam,kap) 0.5*(1+sign(x-m)) - sign(x-m).*(kap.^-sign(x-m)/lam).*pdf_al(x,m,lam,kap);
            pdf_tl = @(x,m,lam,kap) pdf_al(x,m,lam,kap) ./ (cdf_al(prc1(2),m,lam,kap)-cdf_al(prc1(1),m,lam,kap));
            med1   = median(V1);
            mad1   = mean(abs(V1-med1));
            pr     = mle(V1, 'pdf',pdf_tl, 'start',[med1 mad1 1], 'lower',[-Inf 1E-3 1E-3], 'upper',[Inf 50 Inf]);
            mom(1) = pr(1); %m
            mom(2) = pr(2); %lambda
            mom(3) = pr(3); %kappa
            if (doGraph(1))
               figure('name',[name ' - static distribution']);
               hold on
               [~,cen]  = hist(V1,100);
               hist(V1,100);
               t = pdf_tl(cen,pr(1),pr(2),pr(3));
               n = size(V1,1);
               plot(cen,t./sum(t).*n,'-g','LineWidth',2);
               hold off
               loglik = sum(log(pdf_tl(V1,pr(1),pr(2),pr(3))));
               disp(loglik);
               if (doSave)
                  saveas(gcf,[name '1.jpg']);
               end
            end
         otherwise
            error('wrong type in features');
      end
   end
   
   % dynamics of distribution
   if (doMoms(4) || doMoms(5))
      [b1,st1] = robustfit(V1(good2),FV1(good2));
      mom(4)   = b1(2);
      mom(5)   = st1.robust_s;
      if (doGraph(2))
         figure('name',[name ' - rho']);
         scatter(V1(good2),FV1(good2),'.');
         if (doSave)
            saveas(gcf,[name '2.jpg']);
         end
      end
   end
   if (doMoms(6))
      if(sum(V1<mom(1))>100)
         [b2,~  ] = robustfit(V1(good2 & V1<mom(1)),FV1(good2 & V1<mom(1)));
         mom(6)   = b2(2);
      else
         warning(['problem with rho- at ' name]);
      end
   end
   if (doMoms(7))
      if(sum(V1>mom(1))>100)
         [b3,~  ] = robustfit(V1(good2 & V1>mom(1)),FV1(good2 & V1>mom(1)));
         mom(7)   = b3(2);
      else
         warning(['problem with rho+ at ' name]);
      end
   end
   
   % distribution slope on lK and iqr slope on lK
   if (doMoms(8))
      lK1    = lK(good1);
      b4     = robustfit(lK1,V1);
      mom(8) = b4(2);
      if (doGraph(3))
         figure('name',[name ' - k slope']);
         scatter(lK1,V1,'.');
         if (doSave)
            saveas(gcf,[name '3.jpg']);
         end
      end
   end
   if (doMoms(9))
      lK1    = lK(good1);
      b5     = scaler(lK1,V1,5,20,30,doGraph(4)); % see scaler.m for options description
      mom(9) = b5(2);
   end

   % return correct column vector
   mom = mom(doMoms==1)';
end
