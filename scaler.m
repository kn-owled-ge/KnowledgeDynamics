% ----------------------------------------------------
% scaler:
%   calculate the way Y's scale changes with respect
%   to X.
%   method: which method to use for scale
%   1 = iqr, 2 = mean absolute deviations
%   3 = median absolute deviations, 4 = standard deviation
%   5-8 = repeat 1-4, but in logs
%
%
% By Robert Parham, 2015
%
% ----------------------------------------------------

function [b, st] = scaler(X, Y, method, step, minbin, doGraph)

% set default parameter values
switch (nargin)
   case 2
      method  = 1;  % 1=iqr, 2=mean ad, 3=median ad, 4=std, 5-8=with log
      step    = 20;
      minbin  = 30;
      doGraph = 0;
   case 3
      step    = 20;
      minbin  = 30;
      doGraph = 0;
   case 4
      minbin  = 30;
      doGraph = 0;
   case 5
      doGraph = 0;
end

% verify data conforms
assert (sum(isnan(X) | isinf(X))==0)
assert (sum(isnan(Y) | isinf(Y))==0)
assert (minbin>0 & step>0 & method>=1 & method<=8)

% prepare vars
bins  = round(X.*step)./step;
Ubins = unique(bins)';
use   = (bins*ones(size(Ubins)) == ones(size(bins))*Ubins)*1;
use(~use) = nan;
Y1 = Y*ones(size(Ubins));

switch method
   case 1
      lmet = iqr(Y1.*use);
   case 2
      lmet = mad(Y1.*use);
   case 3
      lmet = mad(Y1.*use,1);
   case 4
      lmet = std(Y1.*use);
   case 5
      lmet = log(iqr(Y1.*use));
   case 6
      lmet = log(mad(Y1.*use));
   case 7
      lmet = log(mad(Y1.*use,1));
   case 8
      lmet = log(std(Y1.*use));
end
count = sum(~isnan(use));

% only use bins with sufficient obs in regression
Cbins = Ubins(count>=minbin);
Clmet = lmet (count>=minbin);

% run regression
[b,st] = robustfit(Cbins,Clmet);
if (doGraph)
   figure('name','scaler');
   scatter(Cbins,Clmet,'.');
end

end

