
%set pwd
disp('starting');
cd('/gpfs/fs1/home/rparham2/Documents/MATLAB');

% setup cluster
pc=parcluster('local');
JOB_ID=getenv('SLURM_JOBID');
CPUS=str2num(getenv('SLURM_CPUS_PER_TASK'));
pc.JobStorageLocation=strcat('/local_scratch/',JOB_ID);
parpool(pc,CPUS);

%%

% Prepare params
D.method   = 'pattern';                          % 'gradient'/'anneal'/'pattern'/'multis'
D.useVC    = false;                              % use bootstrapped VC matrix?

%          [1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16     17    ]
%          [theta  omega  rhoZ   muZ    sigZ   delK   lambK  delN   lambN  alpha  beta   gama   tauC   tauN   Pext   MQent  SQent ] 
D.initX1 = [0.900  0.305  0.908  0.000  0.291  0.112  0.331  0.276  0.171  0.001  0.876  0.070  0.350  0.000  0.095  2.000  1.000 ];

%         [theta  omega  rhoZ   sigZ   lambK  delN   lambN  beta ]
D.initX = [0.900  0.305  0.908  0.291  0.331  0.276  0.171  0.876];
D.lbX   = [0.700, 0.100, 0.800, 0.050, 0.050, 0.100, 0.050, 0.800];
D.ubX   = [0.990, 0.900, 0.990, 0.900, 0.950, 0.700, 0.950, 0.970];

%%
%{
% Prepare starting points
N    = 1000;
L    = size(D.initX,2);
b    = zeros(N,L);
good = true(N,1);
over = D.initX1;
fid  = 1;

for i=1:N
   b(i,:) = D.lbX + rand(1,L).*(D.ubX-D.lbX);
   
   over(1)  = b(i,1);
   over(2)  = b(i,2);
   over(3)  = b(i,3);
   over(5)  = b(i,4);
   over(7)  = b(i,5);
   over(8)  = b(i,6);
   over(9)  = b(i,7);
   over(11) = b(i,8);

   [~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, fail] = calibF(over,fid);
   good(i) = ~fail;
end
sum(good)
ptmatrix = b(good,:);
ptmatrix = [ptmatrix ; D.initX];
save ('ptmatrix', 'ptmatrix');
clear good over N L b 
%}

%%

% Get a P structure for all the defs (and open log)
fid    = 1;
P      = calibF(D.initX1,fid);

% get data
Dfile = 'data_160508BOTH.csv';          % which data file to use to create data moments
D.data = csvread(Dfile);
D.Dmom = distmoms(P, D.data);

% Set D members
D.fid  = P.fid;
D.moms = 1:P.Nmoms;

% get weight matrix
if (D.useVC)
   VCname = ['VC.mat'];
   load(VCname)
   D.W = inv(VC(D.moms,D.moms));
else
   %wgts = ones(size(D.moms));
   %D.W  = diag((abs(D.Dmom(D.moms)).^-2)'.*wgts)./sum(wgts);
   D.W  = diag([1.3 2.1 7.9   1.2 2.1 7.9   2.7 4.0 7.3   2.5 5.3 7.9   2.1 4.7   5.0 4.3   7.9 7.9 7.9 7.9]);
   %clear wgts;
end

% prepare function pointer for fmincon
Qfunc = @(b) Qobj3(D, b);
%%
warning('off','all');
warning('off','stats:robustfit:RankDeficient');

% run optimization procedure
if (strcmp(D.method,'gradient'))
   options = optimoptions('fmincon');
   options = optimoptions(options, 'Display'       , 'iter'           );  % 'iter'/'final'
   options = optimoptions(options, 'MaxFunEvals'   , 5000             );
   options = optimoptions(options, 'FinDiffRelStep', 0.01             );
   options = optimoptions(options, 'Algorithm'     , 'interior-point' );  % 'interior-point'/'active-set'/'sqp'
   options = optimoptions(options, 'TolX'          , 1E-3             );
   options = optimoptions(options, 'TypicalX'      , D.initX          );
   options = optimoptions(options, 'UseParallel'   , true);
   [X,Q] = fmincon(Qfunc,D.initX,[],[],[],[],D.lbX,D.ubX,[],options);

elseif (strcmp(D.method,'anneal'))
   options = saoptimset('simulannealbnd');
   options = saoptimset(options, 'Display'           , 'iter'           );
   options = saoptimset(options, 'MaxFunEvals'       , 10000            );
   options = saoptimset(options, 'InitialTemperature', 10               );
   options = saoptimset(options, 'TimeLimit'         , Inf              );
   options = saoptimset(options, 'ReannealInterval'  , 125              );
   options = saoptimset(options, 'TemperatureFcn'    , @temperatureexp  );
   options = saoptimset(options, 'StallIterLim'      , 1000             );
   options = saoptimset(options, 'TolFun'            , 1E-3             );
   [X,Q] = simulannealbnd(Qfunc,D.initX,D.lbX,D.ubX,options);

elseif (strcmp(D.method,'pattern'))
   options = psoptimset(@patternsearch);
   options = psoptimset(options, 'Display'            , 'iter');
   options = psoptimset(options, 'MaxFunEvals'        , 10000 );
   options = psoptimset(options, 'TolX'               , 1E-3);
   options = psoptimset(options, 'TolFun'             , 1E-3);
   options = psoptimset(options, 'TolMesh'            , 1E-4);
   options = psoptimset(options, 'MeshAccelerator'    , 'on');
   options = psoptimset(options, 'SearchMethod'       , @GPSPositiveBasisNp1);
   options = psoptimset(options, 'PollMethod'         , 'GPSPositiveBasis2N');
   options = psoptimset(options, 'TimeLimit'          , 3*24*60*60); %3 days
   options = psoptimset(options, 'CompletePoll'       , 'on');
   options = psoptimset(options, 'CompleteSearch'     , 'on');
   options = psoptimset(options, 'InitialMeshSize'    , 0.125);
   options = psoptimset(options, 'UseParallel'        , true);
   [X,Q] = patternsearch(Qfunc,D.initX,[],[],[],[],D.lbX,D.ubX,[],options);

elseif (strcmp(D.method,'multis'))
   options = optimoptions('fmincon');
   options = optimoptions(options, 'Display'       , 'off'            );  % 'off'/'iter'/'final'
   options = optimoptions(options, 'MaxFunEvals'   , 50               );
   options = optimoptions(options, 'FinDiffRelStep', 0.01             );
   options = optimoptions(options, 'Algorithm'     , 'interior-point' );  % 'interior-point'/'active-set'/'sqp'
   options = optimoptions(options, 'TolX'          , 1E-3             );
   options = optimoptions(options, 'TypicalX'      , D.initX          );
   prob    = createOptimProblem('fmincon', 'x0', D.initX, 'objective', Qfunc, 'lb', D.lbX, 'ub', D.ubX, 'options', options);
   ms      = MultiStart('Display'            , 'iter',      ...
                        'FunctionTolerance'  , 1E-2,        ...
                        'XTolerance'         , 1E-2,        ...
                        'MaxTime'            , 3*24*60*60,  ...
                        'StartPointsToRun'   , 'bounds',    ...
                        'UseParallel'        , true            );
   load('ptmatrix.mat');
   tpoints = CustomStartPointSet(ptmatrix);
   [X,Q,flagm,outptm,allmins] = run(ms,prob,tpoints);
   disp(flagm);
   disp(outptm);
   disp(allmins);
else
   error('other methods unimplemented');
end

warning('on','all');

disp(X);
disp(Q);

if (D.fid ~= 1)
   fclose(D.fid);
end


delete(gcp('nocreate'))
exit
