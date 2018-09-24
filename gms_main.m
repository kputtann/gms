% Main File for Control Design using Generalized Mixed Sensitivity(GMS)
%
% Computes a H-infinity based Feedback Controller based on
% multiobjective constrained convex optimization.
%
% Outline of steps for GMS problem setup:
%   - Form the design plant:
%       - Define the original plant
%       - Integrator augmentation if needed
%       - Bilinear transformation values if needed
%
%   - Select weighting functions:
%       - Tradeoff param rho
%       - W for obj
%       - W for constraint
%
%   - Select optimization params:
%       - LB and UB
%       - Init point
%       - Maximum number of iterations
%
%   - Select Youla/Zames parametrization:
%       - Select Youla or Zames
%       - Initial controller
%
%   - Finite Dimensionality
%       - Basis params
%
%   - Objective function:
%           - sum/max/stacking
%
%   - Find initial controller (Ko, F, L)
%
%   - Youla parameterization
%
%   - Find Initial Q parameter using initial controller (Ko, F, L)
%
%   - Extract required data from problem setup
%
%   - Vectorize the optimization problem
%
%   - Optimization process
%       - define how subgradient is picked based on sum/max/stacking
%
%   - form Q using the optimized variables and bases
%
%   - form Controller K using the obtained Q
%
%   - Inverse bilinear transformation if needed
%
%   - Inverse of integrator augmentation if needed
%
%   - Compute OL and CL maps

%% Initial Code Setup
clear;
close all;
% clc;
% warning off;

% Transfer function variable
s = tf('s');

%% Design Plant

% % ------- Select one of the following PlntLabel ------- %
% PlntLabel='SISO_Stable'; Bilinear=0; AugInteg=0;
PlntLabel='acad_2by2'; Bilinear=0; AugInteg=0;
% PlntLabel='hsv_io'; Bilinear=1; AugInteg=0; AugTwoChannel=1;
% PlntLabel='1bys'; Bilinear=0; AugInteg=0;

% %-------------- Plant --------------%
switch PlntLabel
    case 'SISO_Stable'
        P_tf = tf([1],[1 1]);
        P_ss = ss(P_tf);
        [Ap, Bp, Cp, Dp] = ssdata(P_ss);
        
    case 'acad_2by2'
        s=tf('s');
        P_tf = 1/s * [10 9; 9 8]; % Doyle Example
        P_orig = P_tf;
        %         rolloff = 100/(s+100);
        %         P_tf = series(rolloff,P_tf);
        P_design = P_tf;
        P_ss=ss(P_tf);
        [Ap, Bp, Cp, Dp] = ssdata(P_ss);
        
    case 'hsv_io'
        %         % %------------ Flexible model ------------%
        %         A = [
        %             -0.0008 -0.0006  0.0000  0.0001 -0.0005 -0.0000 -0.0012...
        %             0  -0.0009  0
        %             0.0398  -0.1068 -0.0000  0.1068 -0.0565  0.0002 -0.0502...
        %             0  -0.1122  0
        %             -8.1371 -6.4875 -0.0018  6.4875 -2.6934 -0.0500 -13.9072...
        %             0  -3.5889  0
        %             0       0       1       0       0       0       0       ...
        %             0   0       0
        %             0       0       0       0       0       1       0       ...
        %             0   0       0
        %             97.7307 -175.4193  0  175.4193 -486.4654 -0.7903 -62.4896...
        %             0 -194.3275  0
        %             0       0       0       0       0       0       0       ...
        %             1   0       0
        %             -29.74  -8.29    0       8.29    6.187  0       1.578   ...
        %             0  -8995 -3.796
        %             ];
        %         B = [
        %             0.07128  -0.0006543
        %             0.242     0.01766
        %             -34.65    -9.551
        %             0         0
        %             0         0
        %             -20.3    39.96
        %             0         0
        %             176.1   -25.53
        %             0         0
        %             -94.96   -4.327
        %             ];
        %         C = [
        %             1    0    0    0    0    0    0    0    0    0
        %             0    1    0    0    0    0    0    0    0    0
        %             ];
        %         D = [
        %             0    0
        %             0    0
        %             ];
        %
        %         P_ss=ss(Ap,Bp,Cp,Dp);
        %         P_ss_TwoChannel=P_ss; % backup the original plant (2-output)
        %
        %         % Parameters for bilinear transformation
        %         p2 = -1e20; p1 = -0.001;
        %
        %         % Augment integrator at output, in first two channels
        %         AugTwoChannel=1;
        %         if AugTwoChannel==1
        %             s=tf('s');
        %             [Ap,Bp,Cp,Dp]=ssdata(P_ss);
        %
        %             % Make it 3-channel.
        %             Cps=[P_ss_TwoChannel.c; 0 0 0 1 zeros(1,6)];
        %             % First two states are of integrator
        %             P_ss_ThreeChannel=ss(Ap,Bp,Cps,[]); % backup
        %             M=[0 0 0 1 zeros(1,6)];
        %
        %             % Augment integrator in first two output channels of plant
        %             Ap=blkdiag(Ap,zeros(size(Cp,1)));
        %             Ap(end-size(Cp,1)+1:end,1:size(Cp,2))=Cp;
        %             Bp=[Bp; zeros(size(Cp,1),size(Bp,2))];
        %             Cp=[zeros(2,10), eye(2)];
        %             Cp=[Cp; zeros(1,3) 1 zeros(1,8)];
        %             P_ss=ss(Ap,Bp,Cp,[]);
        %             P_ss_ThreeChannel_AugInteg=P_ss; % Backup
        %         else
        %             % Make it 3-channel
        %             Cps=[P_ss_TwoChannel.c; 0 0 0 1 zeros(1,6)];
        %             P_ss=ss(Ap,Bp,Cps,[]); [Ap,Bp,Cp,Dp]=ssdata(P_ss);
        %             P_ss_ThreeChannel=P_ss; % Backup
        %         end
        
        % %------------ Rigid model ------------%
        Ap=[-0.0008659  -0.0004395   9.981e-09  -0.0001174
            0.02865    -0.08627  -2.467e-06     0.08627
            -8.706      -5.512   -0.001827       5.512
            0           0           1           0];
        Bp =[  0.07122  -0.0006823
            0.242     0.01353
            -35.54      -9.621
            0           0];
        Cp =[1   0   0   0
            0   1   0   0
            0   0   0   1];
        Dp =[0   0
            0   0
            0   0];
        P_ss=ss(Ap,Bp,Cp,Dp);
        P_ss_TwoChannel=P_ss(1:2,:); % backup original plant (2-output)
        
        % Parameters for bilinear transformation
        p2 = -1e20; p1 = -0.001;
        
        % Augment integrator at output, in first two channels
        if AugTwoChannel==1
            s=tf('s');
            [Ap,Bp,Cp,Dp]=ssdata(P_ss_TwoChannel);
            
            % Make it 3-channel.
            P_ss_ThreeChannel=P_ss; % backup
            M=[0 0 0 1];
            
            % Augment integrator in first two output channels of plant
            Ap=blkdiag(Ap,zeros(size(Cp,1)));
            Ap(end-size(Cp,1)+1:end,1:size(Cp,2))=Cp;
            Bp=[Bp; zeros(size(Cp,1),size(Bp,2))];
            Cp=[zeros(2,4), eye(2)];
            Cp=[Cp; zeros(1,3) 1 zeros(1,2)];
            P_ss=ss(Ap,Bp,Cp,[]);
            P_ss_ThreeChannel_AugInteg=P_ss; % Backup
            
        else
            % Make it 3-channel
            Cps=[P_ss_TwoChannel.c; 0 0 0 1];
            P_ss=ss(Ap,Bp,Cps,[]); [Ap,Bp,Cp,Dp]=ssdata(P_ss);
            P_ss_ThreeChannel=P_ss; % Backup 3-output
        end
        
        P_ss=ss(Ap,Bp,Cp,[]);
        P_tf=tf(P_ss);
        
    case '1bys'
        P_tf = tf([1],[1 0]);%*[1 0.1; 0.1 1];
        P_ss = ss(P_tf);
        % s=tf('s'); P_tf=P_tf/s; P_ss=series(ss(0,1,1,0),P_ss);
        [Ap, Bp, Cp, Dp] = ssdata(P_ss);
        
end

% %-------- Integrator augmentation at output if needed -------%
if AugInteg==1
    Integ=ss(0,1,1,0);
    P0=P_ss;
    P_ss=series(Integ,P_ss);
    [Ap, Bp, Cp, Dp] = ssdata(P_ss);
end

% %------------ Size of design plant ------------%
[n_e, n_u] = size(P_ss);

%% Bilinear Transformation
if  Bilinear==1
    P_ss_BeforeBilin=P_ss; % Backup plant before bilin transform
    [Ap,Bp,Cp,Dp]=bilin(P_ss.a,P_ss.b,P_ss.c,P_ss.d,1,'Sft_jw',[p2 p1]);
    P_ss=ss(Ap,Bp,Cp,Dp);
    P_ss_BilinPlnt=P_ss; % backup bilin transformed plant
end

%% Objective Weighting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design parameters: Weights for multiobjective function
% mu corresponds to weight on properties at plant output,
% rho corresponds to weight on properties at plant input
% If eta is defined, it corresponds to weight on properties at sensor noise
% for the inner-outer loop case

mu1=1; mu2=1; mu3=1; rho1=1; rho2=1; rho3=1;
% eta1=1e-1; eta3=0;

switch PlntLabel
    case {'SISO_Stable'}
        Eps=0.01;
        Ms=1.5; wb=10;
        W1 = tf([1/Ms wb], [1 wb*Eps]);
        % % W1 = (1/Ms)*tf([1 nthroot(Ms,k1)*wb],[1 nthroot(Eps,k1)*wb])^k1;
        Mu=1/30; wbu=750;
        W2 = [tf([1 wbu*Mu],[Eps wbu])];
        % % W2 = (1/Eps)*tf([1 wbu*nthroot(Mu,k2)],[1 wbu/nthroot(Eps,k2)])^k2;
        My=30; wbc=1000;
        % W3 = tf([1 wbc/My], [Eps wbc]);
        % % W3 = (1/Eps)*tf([1 wbc/nthroot(My,k3)],[1 wbc/nthroot(Eps,k3)])^k3;
        W3 = ss(1);
        Wd1=W1(1,1);
        wd21=0.1; wd22=1; wd23=10; s=tf('s');
        Wd2=((wd21/(s+wd21))*((s+wd22)/wd22)^2*(wd23/(s+wd23)));
        Wd3=W3(1,1);
        W1=mu1*W1; W2=mu2*W2; W3=mu3*W3;
        W1=ss(W1); W2=ss(W2); W3=ss(W3);
        Wd1=rho1*Wd1; Wd2=rho2*Wd2; Wd3=rho3*Wd3;
        Wd1=ss(Wd1); Wd2=ss(Wd2); Wd3=ss(Wd3);
        
    case {'acad_2by2'}
        Eps=0.01; Ms=2; wb1=1; wb2=1;
        W1 = [tf([1/Ms wb1], [1 wb1*Eps]) 0; 0 tf([1/Ms wb2], [1 wb2*Eps])];
        Eps=0.01; Ms=2; wb1=2; wb2=2;
        Wd1 = [tf([1/Ms wb1], [1 wb1*Eps]) 0; 0 tf([1/Ms wb2], [1 wb2*Eps])];
        W2 = ss(eye(2));
        W3 = ss(eye(2));
        Wd2 = ss(eye(2));
        Wd3 = ss(eye(2));
        W1=mu1*W1; W2=mu2*W2; W3=mu3*W3;
        W1=ss(W1); W2=ss(W2); W3=ss(W3);
        Wd1=rho1*Wd1; Wd2=rho2*Wd2; Wd3=rho3*Wd3;
        Wd1=ss(Wd1); Wd2=ss(Wd2); Wd3=ss(Wd3);
        
    case 'hsv_io'
        % Weights for case: Two loop-breaking points
        Eps1=0.01; Ms1=1.08; wb1=0.01; Ms2=1.08; wb2=0.01;
        W1 = [tf([1/Ms1 wb1], [1 wb1*Eps1]) 0; 0 tf([1/Ms2 wb2],...
            [1 wb2*Eps1])];
        Eps2=1; Mu1=0.1; wbu1=1000; Mu2=0.1; wbu2=1000;
        W2 = [tf([1 wbu1*Mu1],[Eps2 wbu1]) 0; 0 tf([1 wbu2*Mu2],...
            [Eps2 wbu2])];
        Eps3=0.01; My=1.3; wbc=100;
        W3 = tf([1 wbc/My], [Eps3 wbc])*eye(2);
        Epsd1=0.01; Msd1=1.05; wbd1=0.18;
        Wd1=tf([1/Msd1 wbd1], [1 wbd1*Epsd1])*eye(n_u);
        wd21=1; wd22=1; wd23=1;
        Wd2=((wd21/(s+wd21))*((s+wd22)/wd22)^2*(wd23/(s+wd23)))*eye(2);
        Epsd3=0.01; Myd=1.3; wbcd=1000;
        Wd3=tf([1 wbcd/Myd], [Epsd3 wbcd])*eye(n_u);
        W1=mu1*W1; W2=mu2*W2; W3=mu3*W3;
        W1=ss(W1); W2=ss(W2); W3=ss(W3);
        Wd1=rho1*Wd1; Wd2=rho2*Wd2; Wd3=rho3*Wd3;
        Wd1=ss(Wd1); Wd2=ss(Wd2); Wd3=ss(Wd3);
        
        %         % Weights for case: Three loop-breaking points
        %         Eps1=0.01; Ms1=1.08; wb1=0.01; Ms2=1.08; wb2=0.01;
        %         W1 = [tf([1/Ms1 wb1], [1 wb1*Eps1]) 0; 0 tf([1/Ms2 wb2],...
        %             [1 wb2*Eps1])];
        %         Eps2=1; Mu1=0.1; wbu1=2000; Mu2=0.1; wbu2=2000;
        %         W2 = [tf([1 wbu1*Mu1],[Eps2 wbu1]) 0; 0 tf([1 wbu2*Mu2],...
        %             [Eps2 wbu2])];
        %         Eps3=0.01; My=1.3; wbc=100;
        %         W3 = tf([1 wbc/My], [Eps3 wbc])*eye(2);
        %         Epsd1=0.01; Msd1=1.0; wbd1=0.1;
        %         Wd1=tf([1/Msd1 wbd1], [1 wbd1*Epsd1])*eye(n_u);
        %         wd21=1; wd22=1; wd23=1;
        %         Wd2=((wd21/(s+wd21))*((s+wd22)/wd22)^2*(wd23/(s+wd23)))*eye(2);
        %         Epsd3=0.01; Myd=1.2; wbcd=1000;
        %         Wd3=tf([1 wbcd/Myd], [Epsd3 wbcd])*eye(n_u);
        %         Epsni1=0.1;
        %         Mu1=0.001; wbu1=520; Mu2=0.001; wbu2=520;
        %         Wni1 = [(1/sqrt(Epsni1)*tf([1 wbu1*sqrt(Mu1)],[1 wbu1...
        %             /sqrt(Epsni1)]))^2 0; 0 (1/sqrt(Epsni1)*tf([1 ...
        %             wbu2*sqrt(Mu2)],[1 wbu2/sqrt(Epsni1)]))^2];
        %         Wni1 = Wni1(1,1)*eye(n_u);
        %         Wni3 = Wd1(1,1);
        %         W1=mu1*W1; W2=mu2*W2; W3=mu3*W3;
        %         W1=ss(W1); W2=ss(W2); W3=ss(W3);
        %         Wd1=rho1*Wd1; Wd2=rho2*Wd2; Wd3=rho3*Wd3;
        %         Wd1=ss(Wd1); Wd2=ss(Wd2); Wd3=ss(Wd3);
        %         Wni1=eta1*Wni1; Wni1=ss(Wni1);
        %         Wni3=eta3*Wni3; Wni3=ss(Wni3);
        
    case '1bys'
        Eps = 0.00001;
        Ms = 10; wb = 1; k1 = 1;
        W1 = (1/Ms)*(s+nthroot(Ms,k1)*wb)^k1/((s+wb*Eps*0.1/0.1^k1)...
            *(s+0.1*wb)^(k1-1));
        W2 = ss(1);
        My = 10; wbc = 1; k3 = 1;
        W3 = (1/Eps)*(s+wbc/nthroot(My,k3))^k3/((s+10*wbc)^(k3-1)*...
            (s+wbc*10/(Eps*10^k3)));
        Wd1 = W1;
        Wd2 = W2;
        Wd3 = W3;
        
        W1=mu1*W1; W2=mu2*W2; W3=mu3*W3;
        W1=ss(W1); W2=ss(W2); W3=ss(W3);
        Wd1=rho1*Wd1; Wd2=rho2*Wd2; Wd3=rho3*Wd3;
        Wd1=ss(Wd1); Wd2=ss(Wd2); Wd3=ss(Wd3);
end


%% Constraint Weighting Functions

constr_flag = 1;  % 0 = unconstrained, 1 = constrained.

if constr_flag == 0
    W1c=[];
    W2c=[];
    W3c=[];
    Wd1c=[];
    Wd2c=[];
    Wd3c=[];
    
elseif constr_flag == 1
    W1c=[];
    W2c=[];
    W3c=[];
    Wd1c=[];
    Wd2c=[];
    Wd3c=[];
    
    W2c{1}.tfm = ss(1)*eye(n_u);            % Constraint Weigthing
    W2c{1}.Fun = 'f_Hinf';                  % Constraint Type
    W2c{1}.Val = 20;                        % Constraint Value
    
    %     W2c{1}.tfm = ss(1)*eye(n_u);      % Constraint Weigthing
    %     W2c{1}.Fun = 'f_Linf';            % Constraint Type
    %     W2c{1}.Val = 6;                   % Constraint Value
    %
    %     W2c{1}.tfm = ss(1)*eye(n_u);      % Constraint Weigthing
    %     W2c{1}.Fun = 'f_Linf';            % Constraint Type
    %     W2c{1}.Val = [Inf Inf; Inf 24];   % Constraint Value
    
else
    disp('Set constr_flag to indicate unconstrained or constrained')
end

%% Weighting functions data structure
weights.W1 = W1;
weights.W2 = W2;
weights.W3 = W3;
weights.Wd1 = Wd1;
weights.Wd2 = Wd2;
weights.Wd3 = Wd3;

weights.W1c = W1c;
weights.W2c = W2c;
weights.W3c = W3c;
weights.Wd1c = Wd1c;
weights.Wd2c = Wd2c;
weights.Wd3c = Wd3c;

if (exist('Wni1','var'))
    weights.Wni1 = Wni1;
end
if (exist('Wni3','var'))
    weights.Wni3 = Wni3;
end

%% Finite Dimensionality: Basis parameters

Basis.n    = 5;
Basis.type = 2;
Basis.p    = 10;
Basis.z    = 10;

N = Basis.n;

%% Youla/Zames Parameterization
% Youla=1 or Zames=0; Type of parameterization.
% Zames only for stable plant, zero (0) initial controller
YoulaOrZames=1;

%% Optimization Parametrs: Bounds and initial point
xmax1 = 100; xmin1 = -100;
x01 = 1;  % Initial point for optimization xk'
MaxIter = 100;

% Form the vector of LB, UB and initial point
x0=x01*ones(N*n_u*n_e,1);
xmax=xmax1*ones(N*n_u*n_e,1); xmin=xmin1*ones(N*n_u*n_e,1);

%% Optimization to find Controller
SumOrMax =1;  % 1=Max, 2=Sum.

%% Nominal Controller
[Ko,F,L]=f_KNominal(P_ss);

%% Coprime factorization
if YoulaOrZames==1
    % --------- Classic P*K structure (no inner-outer) ---------%
    [T11rz, T12rz, T21rz,T11dz, T12dz, T21dz]=f_CoprFac(P_ss,...
        F,L, weights);
    
    %     %--------- Added for Hypersonic inner-outer ---------%
    %     [T11rz, T12rz, T21rz,T11dz, T12dz, T21dz]=f_CoprFac_hsvio(P_ss,...
    %         F,L, weights);
    
    %     %--------- Added for Hypersonic inner-outer WITH Tniu ---------%
    %     [T11rz, T12rz, T21rz, T11dz, T12dz, T21dz,T11niz,T12niz,...
    %         T21niz]=f_CoprFac_hsvio_Tniu(P_ss,F,L, weights);
    
    %
else
    [T11rz, T12rz, T21rz,T11dz, T12dz, T21dz]=f_CoprFac_ZamesParam...
        (P_ss,F,L, weights);
end

%% Initial Q-parameter
n_x=size(Ap,1); n_e=size(Cp,1); n_u=size(Bp,2);
N = Basis.n;
q = f_Basis(N, Basis.p, Basis.z, Basis.type);
% x0 = x01*ones(N*n_u*n_e,1);
Q = f_FormQN(x0, q, n_u, n_e, N);
xk=x0;

%% Problem Data

% %--------- Classic P*K structure (no inner-outer) ---------%
[n_e, n_u, ProblemDatarz, ProblemDatadz] = f_GenData(P_ss, weights);

% %--------- Added for Hypersonic inner-outer ---------%
% % P_ss1=P_ss(1:2,1:2); n_u=2; n_e=2;% Added for HSV inner-outer
% [n_e, n_u, ProblemDatarz, ProblemDatadz] = f_GenData_hsvio(P_ss, weights);

% --------- Added for Hypersonic inner-outer WITH Tniu ---------%
% [n_e, n_u, ProblemDatarz, ProblemDatadz, ProblemDataniz] = f_GenData_hsvio_Tniu(P_ss, weights);


Datarz=ProblemDatarz; Datadz=ProblemDatadz;
% --------- Added for Hypersonic inner-outer WITH Tniu ---------%
% Datarz=ProblemDatarz; Datadz=ProblemDatadz;  Dataniz=ProblemDataniz;

%% Vectorization
[Mrz, Mobjrz, Mconrz]=f_Vectorize(T11rz,T12rz,T21rz,q,N,n_u,n_e,ProblemDatarz);
[Mdz, Mobjdz, Mcondz]=f_Vectorize(T11dz,T12dz,T21dz,q,N,n_u,n_e,ProblemDatadz);

% [Mrzc, Mobjrzc, Mconrzc] = f_Vectorize(T11rz1c, T11rz1c, T11rz1c, q, N, n_u, n_e, ProblemDatarzc);

% --------- Added for Hypersonic inner-outer WITH Tniu ---------%
% [Mniz, Mobjniz, Mconniz]=f_Vectorize(T11niz,T12niz,T21niz,q,N,n_u,n_e,ProblemDataniz);


%% Optimization process
NQ=N;
N = length(xk);         % Dimension of problem

algo = 2; % 1=ACCPM, 2=Kelley's CPM, 3=SolvOpt.

if algo == 1
    % -------- ACCPM -------- %
    switch SumOrMax
        case 1
            % % Weighted Minmax
            
            [xk,fx,iter_cnt,perf_meas]=...
                f_ACCPM_GenMixSens_Optimizer(N,NQ,xk,Mobjrz,Mobjdz,...
                Mconrz,Mcondz,T11rz, T12rz, T21rz,T11dz, T12dz, ...
                T21dz,Datarz,Datadz,Q,q,n_u,n_e,xmax,xmin,MaxIter);
            
%             %--------- Added for HSV IO WITH Tniu ---------%
%             [xk,fx,iter_cnt,perf_meas]=...
%                 f_ACCPM_GenMixSens_Optimizer_With_Tniu...
%                 (N,NQ,xk,Mobjrz,Mobjdz,Mobjniz,Mconrz,Mcondz,...
%                 Mconniz,T11rz, T12rz, T21rz,T11dz, T12dz, T21dz,...
%                 T11niz,T12niz,T21niz, Datarz, Datadz,Dataniz, Q,q,...
%                 n_u,n_e,xmax,xmin,MaxIter);
%             
        case 2
            % % Weighted Sum
            [xk,fx,iter_cnt,perf_meas]=...
                f_ACCPM_GenMixSens_Optimizer_Sum...
                (N,NQ,xk,Mobjrz,Mobjdz,Mconrz,Mcondz,...
                T11rz, T12rz, T21rz,T11dz, T12dz, T21dz,...
                Datarz, Datadz, Q,q,n_u,n_e,xmax,xmin,MaxIter);
    end
    
elseif algo == 2
    
    % -------- KELLEY'S CPM -------- %
    
    switch SumOrMax
        case 1
            %         Weighted Minmax
            
            [xk,frz,fdz]=f_KelleyCPM_GenMix_Optimizer...
                (N,NQ,xk,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz,T12rz,...
                T21rz,T11dz, T12dz, T21dz, Datarz, Datadz, Q,q,n_u,...
                n_e,MaxIter,xmax,xmin);
            
%             %--------- Added for HSV IO WITH Tniu ---------%
%             [xk,frz,fdz]=f_KelleyCPM_GenMix_Optimizer_With_Tniu...
%                 (N,NQ,xk,Mobjrz,Mobjdz,Mobjniz,Mconrz,Mcondz,T11rz,...
%                 T12rz, T21rz,T11dz, T12dz, T21dz,T11niz,T12niz,...
%                 T21niz, Datarz, Datadz,Dataniz, Q,q,n_u,n_e,MaxIter,...
%                 xmax,xmin);
            
        case 2
            %         Weighted sum
            [xk,fo]=f_KelleyCPM_GenMix_Optimizer_Sum...
                (N,NQ,xk,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz,...
                T12rz, T21rz,T11dz, T12dz, T21dz, Datarz,Datadz,...
                Q,q,n_u,n_e,MaxIter,xmax,xmin);
    end
    
elseif algo == 3
    % -------- SOLVOPT -------- %
    
    opts(1) = -1;       % negative => minimization
    opts(2) = 1e-4;
    opts(3) = 1e-4;
    opts(4) = MaxIter;  % default num iter 15000
    opts(5) = 0;        % 1->verbose, 0->silent
    NQ = N/(n_u*n_e);
    [xk_solvopt,fx_solvopt,opts_solvopt] = solvopt(x0,...
        @(x)solvopt_fval(x,NQ,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz,...
        T12rz, T21rz,T11dz, T12dz, T21dz, Datarz, Datadz, q,n_u,...
    n_e),@(x)solvopt_sg(x,NQ,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz,...
    T12rz, T21rz,T11dz, T12dz, T21dz, Datarz, Datadz, q,n_u,n_e),opts);
    
else
    disp('Choose a valid algorithm')
end

%% Form Q
Q = f_FormQN(xk, q, n_u, n_e, NQ);
disp(' ')

%%
[forz, Gfo] =feval('f_Hinf', Mobjrz, xk, T11rz, T12rz, T21rz, Q, Datarz.ObjVec);
[fodz, Gfo] =feval('f_Hinf', Mobjdz, xk, T11dz, T12dz, T21dz, Q, Datadz.ObjVec);
fx=max([forz,fodz]);

%% Form K
if YoulaOrZames==1
    K=f_FormK(P_ss,Q,F,L); % youla parameterization
else
    K=f_FormK_ZamesParam(P_ss,Q,F,L); % Zames parameterization
end

T_copr.T11rz = T11rz;
T_copr.T12rz = T12rz;
T_copr.T21rz = T21rz;

T_copr.T11dz = T11dz;
T_copr.T12dz = T12dz;
T_copr.T21dz = T21dz;

% --------- Added for Hypersonic inner-outer WITH Tniu ---------%
% T_copr.T11niz = T11niz;
% T_copr.T12niz = T12niz;
% T_copr.T21niz = T21niz;

%% *************** Inverse Bilinear Transformations ***************
if  Bilinear==1
    [Acp1,Bcp1,Ccp1,Dcp1] = ssdata(K);
    K_BeforeInvBilin=K; % Backup K before inverse bilin transformation
    [Atk1,Btk1,Ctk1,Dtk1]=bilin(Acp1,Bcp1,Ccp1,Dcp1,-1,'Sft_jw',[p2 p1]);
    K=ss(Atk1,Btk1,Ctk1,Dtk1);
    P_ss=P_ss_BeforeBilin;
end

%% Analysis of Control Design

% Frequency vector
wvec2=logspace(-3,3,1000); tvec=linspace(0,10,100);

% Form open and closed loop maps
if strcmp(PlntLabel,'hsv_io')
    K_Design=K; % Backup design K (without integ aug)
    K = minreal(K);
    
    if AugTwoChannel==1
        %         % Augment integrator at input, in all channels
        %         K=series(K,1/s);
        % Augment integrator at output, in first two channels
        Kouter=series(1/s,K(:,1:2));
        Kinner=K(:,3:end);
    else
        Kouter=K(:,1:2);
        Kinner=K(:,3:end);
    end
    
    if AugTwoChannel==1
        K=series(blkdiag(1/s, 1/s, 1),K);
    end
    
    %  Add ROLL-OFF if needed 
    K_NoRolloff=K; % backup
    %     K=series(K_NoRolloff,(58/(s+58))^2);
    
    P_ss=P_ss_TwoChannel; % Plant 2-outputs, without integ
    
    [Lo,Li,So,Si,To,Ti,KS,PS,Tniy,Tniu]=f_CLMapInnerOuter_BigK...
        (P_ss_TwoChannel,K,M);
    
    % Modify weights to remove near-zero dummy values
    W1=W1(1:size(So,1),1:size(So,1));
    W3=W3(1:size(To,1),1:size(To,1));
    Wd2=Wd2(1:size(KS,1),1:size(KS,1));
    n_e=2; n_u=2;
    
else
    %     K = minreal(K);
    % Standard feedback (no inner loop)
    [Lo,Li,So,Si,To,Ti,KS,PS] = f_CLTFM(P_ss,K);
end

K_gms=K;


%% Display max_xk and isstable(To)
max_xk=max(abs(xk))

isstab=isstable(To)


%% CL Performance and Robustness

NormInf = mag2db([hinfnorm(So), hinfnorm(Si), hinfnorm(KS), ...
    hinfnorm(PS), hinfnorm(To), hinfnorm(Ti)])
PerformMeasOutOrigWts=norm([W1*So; W2*KS; W3*To],inf);
PerformMeasInOrigWts=norm([Wd1*Si; Wd2*PS; Wd3*Ti],inf);

% % Bandwidth/crossovers
% BW_20_So = min(getGainCrossover(So,0.1))
% BW_20_Si = min(getGainCrossover(Si,0.1))
% BW_20_To = max(getGainCrossover(To,0.1))
% BW_20_Ti = max(getGainCrossover(Ti,0.1))
% BW_0_KS = max(getGainCrossover(KS,1))
% BW_0_PS = max(getGainCrossover(PS,1))
% % BW_0_Tniu = max(getGainCrossover(Tniu,1))
% BW_0_Lo = max(getGainCrossover(Lo,1))
% BW_0_Li = max(getGainCrossover(Li,1))

% % Time domain properties
% To_stepinfo = stepinfo(To);
% % v_ts = To_stepinfo(1,1).SettlingTime
% % gamma_ts = To_stepinfo(2,2).SettlingTime
% ts1 = To_stepinfo(1,1).SettlingTime
% ts2 = To_stepinfo(2,2).SettlingTime
% KS_stepinfo = stepinfo(KS);
% % peak_FER = KS_stepinfo(1,1).Peak
% % peak_elev = KS_stepinfo(2,2).Peak
% u_peak1 = KS_stepinfo(1,1).Peak
% u_peak2 = KS_stepinfo(2,2).Peak
