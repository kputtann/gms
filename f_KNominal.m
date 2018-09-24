function [Ko,F,L]=f_KNominal(P_ss)
% Nominal Controller 
[Ap, Bp, Cp, Dp] = ssdata(P_ss);
n_x=size(Ap,1); n_e=size(Cp,1); n_u=size(Bp,2);
F = lqr(Ap, Bp, 1e0*eye(n_x), 1.5e1*eye(n_u));
L = lqr(Ap',Cp',1e0*eye(n_x), 1.5e1*eye(n_e));
L=L';
Ko = ss(Ap-Bp*F-L*Cp+L*Dp*F, -L, -F, 0);