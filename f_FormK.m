function K=f_FormK(P_ss,Q,F,L)
% Form the controller from Q
% Uses Youla parameterization

Ap=P_ss.a; Bp=P_ss.b; Cp=P_ss.c; Dp=P_ss.d;
Aq = Q.a; Bq = Q.b; Cq = Q.c; Dq = Q.d; % Q - Parameter
n_x=size(Ap,1); n_e=size(Cp,1); n_u=size(Bp,2);

Delta = eye(n_u) - Dq*Dp;
invDelta = inv(Delta);
Ak11 = (Ap-L*Cp)-(Bp-L*Dp)*invDelta*(-Dq*Cp+F);
Ak12 = -(Bp-L*Dp)*invDelta*Cq;
Ak21 = -Bq*Cp+Bq*Dp*invDelta*(-Dq*Cp+F);
Ak22 = Aq+Bq*Dp*invDelta*Cq;
Ak = [Ak11 Ak12; Ak21 Ak22];
Bk = [L-(Bp+L*Dp)*invDelta*Dq;
    Bq+Bq*Dp*invDelta*Dq];
Ck = [invDelta*(-Dq*Cp+F) invDelta*Cq];
Dk = invDelta*Dq;
K = ss(Ak, Bk, Ck, Dk);