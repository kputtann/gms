function [Lo,Li,So,Si,To,Ti,Tru,PS,Tniy,Tniu]=f_CLMapInnerOuter_generic...
    (P,Ki,Ko,Mi)

% Computes the open and closed loop maps for inner-outer loop 
% configuration
% Mi is a matrix with following properties
% num of col = num of states
% num of row = num of states being fed back in inner loop
% Eg: feeding back states 3 and 4:
% Mi = [0 0 1 0 0 0;
%       0 0 0 1 0 0];

[Ap, Bp, Cp, Dp] = ssdata(P);
[Ai, Bi, Ci, Di] = ssdata(Ki);
[Ao, Bo, Co, Do] = ssdata(Ko);

% Open loop
%% Lo
A_Lo=[Ap+Bp*Di*Mi                  Bp*Co                        -Bp*Ci; 
      zeros(size(Ao,1),size(Ap,2)) Ao                           ...
      zeros(size(Ao,1),size(Ai,2)); 
      Bi*Mi                        zeros(size(Ai,1),size(Ao,2)) Ai];
  
B_Lo=[Bp*Do;
      Bo;
      zeros(size(Bi,1),size(Bo,2))];
C_Lo = [Cp-Dp*Di*Mi  Dp*Co -Dp*Ci];
D_Lo = Dp*Do;
        
Lo=ss(A_Lo,B_Lo,C_Lo,D_Lo);

%% Li
A_Li=[ Ap       zeros(size(Ap,1),size(Ao,2)) zeros(size(Ap,1),size(Ai,2)); 
      -Bo*Cp    Ao                           zeros(size(Ao,1),size(Ai,2)); 
      Bi*Mi     zeros(size(Ai,1),size(Ao,2)) Ai                         ];
  
B_Li=[Bp;
      -Bo*Dp;
      zeros(size(Ai,1),size(Bp,2))];

C_Li=[-Do*Cp-Di*Mi Co -Ci];
D_Li= -Do*Dp;
Li=ss(A_Li,B_Li,C_Li,D_Li);

%% Closed loop
Q =inv(eye(size(Do,1)) + Do*Dp);

Amat=[Ap+Bp*Q*(Do*Dp*Di*Mi-Do*Cp)-Bp*Di*Mi          Bp*Q*Co          ...
    -Bp*Ci+Bp*Q*Do*Dp*Ci     ; 
      Bo*Dp*Di*Mi-Bo*Cp-Bo*Dp*Q*(Do*Dp*Di*Mi-Do*Cp) Ao-Bo*Dp*Q*Co    ...
      Bo*Dp*Ci-Bo*Dp*Q*Do*Dp*Ci;
      Bi*Mi                                         zeros(size(Ai,1),...
      size(Ao,2)) Ai                       ];

%% To
B_To = [Bp*Q*Do                      ;
        Bo-Bo*Dp*Q*Do                ;
        zeros(size(Ai,1),size(Bo,2))];
C_To=[Cp+Dp*Q*(Do*Dp*Di*Mi-Do*Cp)-Dp*Di*Mi       Dp*Q*Co         ...
    Dp*Q*Do*Dp*Ci];
D_To = Dp*Q*Do;
To=ss(Amat,B_To,C_To,D_To);

%% So
D_So=eye(size(Dp,1),size(Bo,2))-Dp*Q*Do;
So=ss(Amat,B_To,-C_To,D_So);

%% KS (Tru, strictly speaking)
C_KS = [Do*Dp*Di*Mi-Do*Cp-Di*Mi Co Do*Dp*Ci-Ci];
D_KS = Q*Do;
Tru=ss(Amat,B_To,C_KS,D_KS);

%% Si
B_Si=[Bp - Bp*Q*Do*Dp              ; 
      -Bo*Dp+Bo*Dp*Q*Do*Dp         ;
      zeros(size(Bi,1),size(Dp,2))];
D_Si=eye(size(Q,1))- Q*Do*Dp;
Si=ss(Amat,B_Si,C_KS,D_Si);

%% Ti
D_Ti = -Q*Do*Dp ;
Ti=ss(Amat,B_Si,C_KS,D_Ti);

%% PS
D_PS = Dp-Dp*Q*Do*Dp ;
PS=ss(Amat,B_Si,C_To,D_PS);

%% Tniy
B_Tniy=[Bp*Q*Do*Dp*Di-Bp*Di;
    Bo*Dp*Di-Bo*Dp*Q*Do*Dp*Di;
    Bi];
D_Tniy=[Dp*Q*Do*Dp*Di-Dp*Di];
Tniy=ss(Amat,B_Tniy,C_To,D_Tniy);

%% Tniu
D_Tniu=Q*Do*Dp*Di-Di;
Tniu=ss(Amat,B_Tniy,C_KS,D_Tniu);

