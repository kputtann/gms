function [Lo,Li,So,Si,To,Ti,Tru,PS,Tniy,Tniu]=f_CLMapInnerOuter_BigK(P,...
    K,Mi)

% Computes the open and closed loop maps for inner-outer loop configuration
% This code can be used when both Ko and Ki have same A and C matrices.
% In other words, when big K is formed which has both Ko and Ki in it.
% This code must be used rather than f_CLMapInnerOuter_generic.m in order
% to avoid nonminimality in the system
% Mi is a matrix with following properties
% num of col = num of states
% num of row = num of states being fed back in inner loop
% Eg: feeding back states 3 and 4:
% Mi = [0 0 1 0 0 0;
%       0 0 0 1 0 0];

[Ap, Bp, Cp, Dp] = ssdata(P);
[Ak, Bk, Ck, Dk] = ssdata(K);

Bo=K.b(:,1:size(Cp,1)); % Number of inputs to Ko is same as num of plant...
% output
Do=K.d(:,1:size(Cp,1));
Bi=K.b(:,size(Cp,1)+1:end); % Inputs to Ki are the remaining inputs to K...
% Also equal to num of rows in Mi
Di=K.d(:,size(Cp,1)+1:end);
% 
% Bo=K.b(:,1:size(Mi,1)); % Number of inputs same as num of rows in Mi
% Do=K.d(:,1:size(Mi,1));
% Bi=K.b(:,size(Mi,1)+1:end); % Inputs are the remaining inputs to K...
% Also equal to num of plant output
% Di=K.d(:,size(Mi,1)+1:end);

% Open loop
%% Lo
A_Lo=[Ap-Bp*Di*Mi                  Bp*Ck; 
      -Bi*Mi                       Ak  ];
B_Lo=[Bp*Do;
      Bo  ];
C_Lo = [Cp-Dp*Di*Mi Dp*Ck];
D_Lo = Dp*Do;
        
Lo=ss(A_Lo,B_Lo,C_Lo,D_Lo);

%% Li
A_Li=[ Ap           zeros(size(Ap,1),size(Ak,2));
      -Bo*Cp-Bi*Mi  Ak                         ];
B_Li=[Bp;
      -Bo*Dp];
C_Li=[-Do*Cp-Di*Mi Ck];
D_Li= -Do*Dp;
Li=ss(A_Li,B_Li,C_Li,D_Li);

%% Closed loop
Q =inv(eye(size(Dp,1)) + Dp*Do); % 

Amat=[Ap-Bp*Do*Q*Cp+Bp*Do*Q*Dp*Di*Mi-Bp*Di*Mi  Bp*Ck-Bp*Do*Q*Dp*Ck; 
      -Bo*Q*Cp+Bo*Q*Dp*Di*Mi-Bi*Mi              Ak-Bo*Q*Dp*Ck     ];
      
%% To
B_To = [Bp*Do*Q;
        Bo*Q  ];
C_To=[Cp-Dp*Do*Q*Cp+Dp*Do*Q*Dp*Di*Mi-Dp*Di*Mi  Dp*Ck-Dp*Do*Q*Dp*Ck];
D_To = Dp*Do*Q;
To=ss(Amat,B_To,C_To,D_To);

%% So
D_So=eye(size(Dp,1))-Dp*Do*Q;
So=ss(Amat,B_To,-C_To,D_So);

%% KS (Tru, strictly speaking)
C_KS = [-Do*Q*Cp+Do*Q*Dp*Di*Mi-Di*Mi Ck-Do*Q*Dp*Ck];
D_KS = Do*Q;
Tru=ss(Amat,B_To,C_KS,D_KS);

%% Si
B_Si=[Bp-Bp*Do*Q*Dp;
      -Bo*Q*Dp    ];
D_Si=eye(size(Do,1))-Do*Q*Dp;
Si=ss(Amat,B_Si,C_KS,D_Si);

%% Ti
D_Ti =-Do*Q*Dp ;
Ti=ss(Amat,B_Si,C_KS,D_Ti);

%% PS
D_PS = Dp-Dp*Do*Q*Dp ;
PS=ss(Amat,B_Si,C_To,D_PS);

%% Tniy
B_Tniy=[Bp*Do*Q*Dp*Di-Bp*Di;
    Bo*Q*Dp*Di-Bi         ];
D_Tniy=-Dp*Di+Dp*Do*Q*Dp*Di;
Tniy=ss(Amat,B_Tniy,C_To,D_Tniy);

%% Tniu
D_Tniu=Do*Q*Dp*Di-Di;
Tniu=ss(Amat,B_Tniy,C_KS,D_Tniu);
