%%
% Initial cleanup
clear all; format short; format compact; beep on; close all;
warning('on','all');

%%

%      [1      2      3      4      5      6      7      8      9      10     11     12     13     14     15     16     17    ]
%      [theta  omega  rhoZ   muZ    sigZ   delK   lambK  delN   lambN  alpha  beta   gama   tauC   tauN   Pext   MQent  SQent ] 
over = [0.728  0.500  0.933  0.000  0.384  0.110  0.377  0.300  0.400  0.001  0.892  0.070  0.350  0.000  0.095  2.000  1.000 ];
fid  = 1; % {1 = screen, NaN = file}

% get data
Dfile = 'data_160508BOTH.csv';          % which data file to use to create data moments
Ddata = csvread(Dfile);
P     = calibF(over,fid);
Data  = Ddata;


%%
% Run model
[Smom, Sdata, Dmom] = modelF(over,fid,Ddata);

beep;

%%
Data = Sdata;

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


%%
%{
% get VC, s.e. and QR decomposition
VC      = bootstrapVC(P, Ddata, 'VC.mat');
se      = sqrt(diag(VC));
CR      = VC./((se*ones(1,P.Nmoms)).*(ones(P.Nmoms,1)*se'));
[~,~,E] = qr(CR,0);
[~,IX]  = sort(E);
disp(IX');

% choose a subset
Nuse = 1:size(VC,1);
CR1 = CR(Nuse,Nuse);
[~,R1,E1] = qr(CR1,0);
[~,IX1]  = sort(E1);
disp(' ');
disp(IX1');

% show weights
wgt = abs(diag(R1));
wgt = wgt./sum(wgt);
disp(' ');
disp(wgt);
%}
