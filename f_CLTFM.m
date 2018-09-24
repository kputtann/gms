function [Lo,Li,So,Si,To,Ti,KS,SP] = f_CLTFM(P,K)

% OL and CL frequency responses (ss based)
% Works for MIMO P and K
% Inputs: 
% P: Plant in state space form
% K: Control in state space form
% Outputs:
% Lo, Li: Open loop tfs in ss
% So,Si,To,Ti,KS,SP: Closed loop tfs in ss

[Ap, Bp, Cp, Dp] = ssdata(P);
n_e = size(P,1);
n_u = size(P,2);
n_p = size(P,'order');
[Ak, Bk, Ck, Dk] = ssdata(K);
n_k = size(K,'order');

%% Lo = PK
A_Lo = [Ap Bp*Ck; zeros(n_k,n_p) Ak]; 
B_Lo = [Bp*Dk; Bk];
C_Lo = [Cp Dp*Ck]; 
D_Lo = Dp*Dk;
Lo = ss(A_Lo,B_Lo,C_Lo,D_Lo);

%% Li = KP
A_Li = [Ak Bk*Cp; zeros(n_p,n_k) Ap]; 
B_Li = [Bk*Dp; Bp];
C_Li = [Ck Dk*Cp]; 
D_Li = Dk*Dp;
Li = ss(A_Li,B_Li,C_Li,D_Li);

%% Mo
Mo = inv(eye(n_e)+Dp*Dk);
%% Mi
Mi = inv(eye(n_u)+Dk*Dp);

%% So = inv(I+PK)
A_So = [Ap-Bp*Dk*Mo*Cp Bp*Ck-Bp*Dk*Mo*Dp*Ck; -Bk*Mo*Cp Ak-Bk*Mo*Dp*Ck]; 
B_So = [Bp*Dk*Mo; Bk*Mo];
C_So = [-Mo*Cp -Mo*Dp*Ck]; 
D_So = Mo;
So = ss(A_So,B_So,C_So,D_So);

%% Si = inv(I+KP)
A_Si = [Ak-Bk*Dp*Mi*Ck Bk*Dp*Mi*Dk*Cp-Bk*Cp; Bp*Mi*Ck Ap-Bp*Mi*Dk*Cp]; 
B_Si = [-Bk*Dp*Mi; Bp*Mi];
C_Si = [Mi*Ck -Mi*Dk*Cp]; 
D_Si = Mi;
Si = ss(A_Si,B_Si,C_Si,D_Si);

%% To = PKinv(I+PK)
A_To = [Ap-Bp*Dk*Mo*Cp Bp*Ck-Bp*Dk*Mo*Dp*Ck; -Bk*Mo*Cp Ak-Bk*Mo*Dp*Ck]; 
B_To = [Bp*Dk*Mo; Bk*Mo];
C_To = [Mo*Cp Mo*Dp*Ck]; 
D_To = Mo*Dp*Dk;
To = ss(A_To,B_To,C_To,D_To);

%% Ti = inv(I+KP)KP
A_Ti = [Ak-Bk*Dp*Mi*Ck Bk*Dp*Mi*Dk*Cp-Bk*Cp; Bp*Mi*Ck Ap-Bp*Mi*Dk*Cp]; 
B_Ti = [-Bk*Dp*Mi; Bp*Mi];
C_Ti = [Mi*Ck -Mi*Dk*Cp]; 
D_Ti = -Dk*Dp*Mi;
Ti = ss(A_Ti,B_Ti,C_Ti,D_Ti);

%% KS 
A_ks = [Ap-Bp*Dk*Mo*Cp Bp*Ck-Bp*Dk*Mo*Dp*Ck; -Bk*Mo*Cp Ak-Bk*Mo*Dp*Ck]; 
B_ks = [Bp*Dk*Mo; Bk*Mo];
C_ks = [-Dk*Mo*Cp Ck-Dk*Mo*Dp*Ck]; 
D_ks = Dk*Mo;
KS = ss(A_ks,B_ks,C_ks,D_ks);

%% SP 
A_sp = [Ak-Bk*Dp*Mi*Ck Bk*Dp*Mi*Dk*Cp-Bk*Cp; Bp*Mi*Ck Ap-Bp*Mi*Dk*Cp]; 
B_sp = [-Bk*Dp*Mi; Bp*Mi];
C_sp = [Mo*Dp*Ck Mo*Cp];  
D_sp = Mo*Dp;
SP = ss(A_sp,B_sp,C_sp,D_sp);