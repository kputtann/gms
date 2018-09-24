% MATLAB code for mu-synthesis using DK-iteration technique 
% Uses MATLAB's dksyn command in Robust Control Toolbox

clear; close all; 

% Transfer function variable
s = tf('s');

%% Plant model
P = 1/s * [10 9; 9 8]; 
[n_e,n_u] = size(P);

%% Define uncertainty blocks
InputUnc = ultidyn('InputUnc',[2 2],'Bound',1);
OutputUnc = ultidyn('OutputUnc',[2 2],'Bound',1);

%% Weighting Functions
% Weight on Divisive Uncertainty at Input
Wi = 1*eye(n_u);
% Weight on Divisive Uncertainty at Output
Wo = 1*eye(n_e);
% Weight on Se
rho = 1;
Eps=0.01;
Ms=2; wb=0.1;
W1 = tf([1/Ms wb], [1 wb*Eps])*eye(n_e);
% Weight on KSe
W2 = ss(0.1)*eye(n_u); 

%% Generalized plant
systemnames = 'P W1 W2 Wi Wo';
inputvar = '[di(2); do(2); w(2); u(2)]';
outputvar = '[Wi; Wo; W1; W2; w-P+do]';  
input_to_P = '[u-di]';
input_to_W1 = '[w-P+do]';
input_to_W2 = '[u-di]';
input_to_Wi = '[u-di]';
input_to_Wo = '[P-do]';
sysoutname = 'GenP';
cleanupsysic = 'yes';
sysic; 
% GenP = minreal(ss(GenP));
GenP_unc = lft([InputUnc zeros(2); zeros(2) OutputUnc],GenP);

%% dksyn command
[K_dksyn,CL_dksyn,Bnd_dksyn,Info_dksyn] = dksyn(GenP_unc,2,2);

%% mu-analysis
% del_I, del_O for robust stability (RS)
BlockStructureRS = [2 2; 2 2];   
% % del_I, del_O for NP
% BlockStructureNP = [2 2];
% % del_I, del_O, del_P for RP
% BlockStructureRP = [BlockStructureRS; BlockStructureNP];  

wvec1=logspace(-3,3,1000);  % freq vec for Mu-analysis
N = lft(GenP,K_dksyn); 
Nf = frd(N,wvec1);  
% mu for RS: 
MuData=mussv(Nf(1:sum(BlockStructureRS(:,1)),1:sum(BlockStructureRS(:,1))),BlockStructureRS);  % Pick the channels corresponding to del_I and del_O. Reject the channel corresponding to del_P
muRS(1:length(wvec1))=MuData(1,1).ResponseData(1,1,:);
% % mu for RP: 
% MuData=mussv(Nf,BlockStructureRP);
% muRP(1:length(wvec1))=MuData(1,1).ResponseData(1,1,:);
% % mu for NP: 
% MuData=mussv(Nf(end-sum(BlockStructureNP(:,1))+1:end,end-sum(BlockStructureNP(:,1))+1:end),BlockStructureNP);  % Pick the channels corresponding to del_I and del_O. Reject the channel corresponding to del_P
% muNP(1:length(wvec1))=MuData(1,1).ResponseData(1,1,:);

figure; 
semilogx(wvec1,mag2db(muRS)); grid on;
title('\mu for Robust Stability');
ylabel('\mu');
xlabel('Frequency (rad/s)');
plot_axis
% plot_legend(hL,hObj)
ylim([-40 20]);

%% CL maps
[Lo,Li,So,Si,To,Ti,KS,PS] = f_CLTFM(P,K_dksyn);

wvec=logspace(-3,3,10000);  % freq vec for plotting

figure; sigma(So,wvec); grid on; hold on; sigma(inv(W1),wvec);
title('S_e'); 
[hL,hObj]=legend('S_e','W_1^{-1}');
plot_axis
plot_legend(hL,hObj)
ylim([-100 20]);

figure; sigma(KS,wvec); grid on; hold on; sigma(inv(W2),wvec);
title('KS_e'); 
[hL,hObj]=legend('KS_e','W_2^{-1}');
plot_axis
plot_legend(hL,hObj)
ylim([-100 20]);
