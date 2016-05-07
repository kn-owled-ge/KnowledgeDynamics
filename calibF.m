% ----------------------------------------------------
% CALIBF:
%   Calibrate the economy in the yearly frequency
% 
% By Robert Parham, 2015
% Based on code by Lu Zhang, 2001
%
% ----------------------------------------------------
function [P, Qz, Qx, Qw, k, kp, n, np, z, x, w, fail] = calibF(over,fid)

% Open log file (for appending) or write to screen/file provided externally
P.sMod    = '1';
if (isnan(fid))
   P.fid = fopen(['log' P.sMod '.txt'], 'a');  
else
   P.fid = fid;
end
clear fid;

% Print params being used
%{
fprintf(P.fid,'Params: theta  omega  rhoZ   sigZ   tau    lamb1k lamb2k lamb3k lamb1n lamb2n lamb3n alpha  beta   pi     MQent  SQent\n');
for i=1:numel(over)
   fprintf(P.fid,'%7.3f\n', over(i));
end
%}

% Technical estimation values
P.VIter   = 300;                                      % stop VFI after this many iterations
P.VFIprn  = 20;                                       % print VFI progress every X iterations
P.SIMprn  = 1000;                                     % print simulation progress every X iterations
P.SIMfin  = false;                                    % print end of simulation statistics
P.SMOMprn = false;                                    % print simulated moments
P.DMOMprn = false;                                    % print data moments
P.EMOMprn = false;                                    % print relative moment errors
P.Vload   = true;                                     % load previous value function?
P.Vsave   = false;                                     % save current value function?

% technology parameter values
P.theta  = over(1);                                  % capital's scale in income
P.rhoZ   = over(2);                                  % persistence of income proc
P.muZ    = over(3);                                  % mean of income proc is hard coded zero
P.sigZ   = over(4);                                  % SD of income proc shock
P.delK   = over(5);                                  % depreciation of physical capital
P.lambK  = over(6);                                  % adjustment elasticity
P.beta   = over(7);                                  % Discount factor=1/(1+r*) (to match data going from DV to V)
P.tauC   = over(8);                                  % wedge from OI to FCF (tax rate*)
P.gama   = over(9);                                  % financing friction coefficient
P.Pext   = over(10);                                 % Probability of random exit per year = pi
P.MQent  = over(11);                                 % Mean of Q at entry
P.SQent  = over(12);                                 % SD of Q at entry

% Financial friction
P.doFin  = true;                                     % Do financial friction?
P.FinEx  = false;                                     % Do firms violating the constraint exit?

% exit values
P.ExEnt = false;                                      % is there exit (and entry) in the model?
P.Vexit = 0;                                          % value scale too low to remain
P.Kexit = 0;                                          % capital scale too low to remain
P.Bexit = 5;                                          % minimal value scale to register as major bankruptcy

% tolerance values
P.Verr   = 1e-2;                                      % Maximum change in value function to stop VFI
P.Nerr   = 1e-1;                                      % Maximum change in N policy to stop VFI
P.Kerr   = 1e-1;                                      % Maximum change in K policy to stop VFI

% simulation parameters
P.Ns     = int32(6000);                               % number of firms 
P.Ts1    = int32(400);                                % number of periods in initial simulation
P.Tb     = int32(300);                                % number of burn-in periods
P.Ts2    = int32(30);                                 % number of periods in moment simulations
P.Reps   = int32(1);                                  % number of repeats for moment simulations
%P.Nmoms  = int32(311);                                % number of moments returned by moment func
P.Nmoms  = int32(15);                                 % number of moments returned by moment func
P.VCreps = int32(100);                                % number of repeats for VC bootstrap

% Construct grid for physical stock
kmin  = 0;
kmax  = 15;
P.nk  = 31;
k     = linspace(kmin,kmax,P.nk)';
P.nkp = P.nk*5;
kp    = linspace(kmin,kmax,P.nkp)';
clear kmin kmax;

% Construct grid for knowledge stock
%nmin   = 0;
%nmax   = 15;
%P.nn   = 31;
%n      = linspace(nmin,nmax,P.nn)';
%P.nnp  = P.nn*5;
%np     = linspace(nmin,nmax,P.nnp)';
%clear nmin nmax;
P.nn  = 1;
P.nnp = 1;
n     = 0;
np    = 0;

% Construct grid for Z
P.nz     = 15;
zdev     = 2*P.sigZ/sqrt((1 - P.rhoZ^2)*(P.nz - 1)); % do ''help rouwTrans'' to understand the formula
[Qz, z]  = rouwTrans(P.rhoZ, P.muZ, zdev, P.nz);
clear zdev;

% Construct grid for X
%P.nx     = 19;
%xdev     = 2*P.sigX/sqrt((1 - P.rhoX^2)*(P.nx - 1));
%[Qx, x]  = rouwTrans(P.rhoX, P.muX, xdev, P.nx);
%clear xdev;
P.nx = 1;
Qx   = 1;
x    = 0;

% Construct grid for W
%P.nw     = 17;
%wdev     = 2*P.sigW/sqrt((1 - P.rhoW^2)*(P.nw - 1));
%[Qw, w]  = rouwTrans(P.rhoW, P.muW, wdev, P.nw);
%clear wdev;
P.nw = 1;
Qw   = 1;
w    = 0;

% Verify model params make sense
fail = false;
r    = (1-P.beta)/P.beta;

% FCF at max(z) for stagnant scale 10 firm is positive, and value is
% within reasonable scale
lK = 10;
OI = exp(z(end) + P.theta*lK);
dK = P.delK*exp(lK);
FC = (1-P.tauC)*OI + P.tauC*P.delK*exp(lK) - dK;
TV = FC/r;
fail = fail | (TV<exp(lK-1)) | (TV>exp(lK+8));

%if(fail)
%   fprintf(P.fid,'Failure when verifying params!\n');
%end

end