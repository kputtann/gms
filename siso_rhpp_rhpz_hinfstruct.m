% MATLAB Code to design controller for SISO unstable and non-minimum 
% phase plant using hinfstruct
% A classic P-K structure OR a hierarchical inner-outer structure can 
% be selected

clear; 
close all;

% Classic P-K structure: set flag=1, or
% Hierarchical inner-outer structure: set flag=2.
flag = 1;

% Transfer functions variable
s = tf('s');

%% Plant
% Plant
p = 1;
z = 10;

P = ((z-s)/(s-p))*((p)/(z));
[n_e,n_u] = size(P);

%% Weighting Functions
Eps=0.01;
Ms=2; wb=0.02; 
% W1 = tf([1/Ms wb], [1 wb*Eps]);
W1 = ss(1);
% Mu=0.1/3; wbu=7500;
% W2 = [tf([1 wbu*Mu],[Eps wbu])];
W2 = ss(0);
My=2; wbc=0.02; 
% W3 = tf([1 wbc/My], [Eps wbc]);
W3 = ss(1);

% Number of exogenous inputs and outputs
n_w = size([W1; W2; W3],2);
n_z = size([W1; W2; W3],1);

%% Define tunable controller structure

% rolloff at a
a = 100;

if flag == 1
    % % Classic P-K structure
    % ---------------------------------------------------
    K_norolloff = ltiblock.pid('K_norolloff','p');
    % K_norolloff = tunablePID('K_norolloff','pi');
    
    % K = gk/s*(s+zk)*(a/(s+a))^2;
    K = K_norolloff;
    % K = series(K_norolloff,(a/(s+a))^2);
    
elseif flag == 2
    % % Hierarchical inner-outer structure
    % ---------------------------------------------------
    Ko_norolloff = ltiblock.pid('Ko_norolloff','p');
    Ko = Ko_norolloff;
%     Ko = series(Ko_norolloff,(a/(s+a))^2);
    
    Ki_norolloff = ltiblock.pid('Ki_norolloff','p');
    Ki = Ki_norolloff;
%     Ki = series(Ki_norolloff,(a/(s+a))^2);
end

%% Define closed-loop interconnection

if flag == 1
    % Classic P-K structure
    % ---------------------------------------------------
    % Using feeedback command
    CL0 = blkdiag(W1*feedback(1,P*K),W3*feedback(P*K,1));
    % CL0 = blkdiag(feedback(1,P*K));
    %
    % % % % Generalized Plant - MSO
    % % systemnames='P W1 W2 W3';
    % % inputvar=['[r{' int2str(n_e) '};u{' int2str(n_u) '}]'];
    % % outputvar=['[W1; W2; W3; r-P]'];
    % % input_to_P='[u]';
    % % input_to_W1='[r-P]';
    % % input_to_W2='[u]';
    % % input_to_W3='[P]';
    % % cleanupsysic='yes';
    % % GenP_mso=sysic;
    % %
    % % [Ag,Bg,Cg,Dg]=ssdata(GenP_mso);
    % % % [Ag,Bg,Cg,Dg]=ssdata(GenP_mso_io);
    % %
    % % % % Matrix blocks of Generalized Plants for LMI
    % % % A=Ag;
    % % % B1=Bg(:,1:n_w);
    % % % B2=Bg(:,n_w+1:end);
    % % % C1=Cg(1:n_z,:);
    % % % C2=Cg(n_z+1:end,:);
    % % % D11=Dg(1:n_z,1:n_w);
    % % % D12=Dg(1:n_z,n_w+1:end);
    % % % D21=Dg(n_z+1:end,1:n_w);
    % % % D22=Dg(n_z+1:end,n_w+1:end);
    % %
    % % % Num. of states of GenP
    % % nx_genp = size(Ag,1);
    % % CL0 = lft(GenP_mso,K);
    
elseif flag == 2
    % Hierarchical inner-outer structure
    % ---------------------------------------------------
    P.InputName = 'u';
    P.OutputName = 'y';
    Ko.InputName = 'e';
    Ko.OutputName = 'uo';
    Ki.InputName = 'y';
    Ki.OutputName = 'ui';
    W1.InputName = 'e';
    W1.OutputName = 'z1';
    W3.InputName = 'y';
    W3.OutputName = 'z3';
    sum_outer=sumblk('e=r-y',1);
    sum_inner=sumblk('u=uo-ui',1);
    WS = connect(P,Ki,Ko,W1,W3,sum_outer,sum_inner,'r','z1');
    WT = connect(P,Ki,Ko,W1,W3,sum_outer,sum_inner,'r','z3');
    CL0 = blkdiag(WS,WT);
    
end
%% Solve $H_{\infty}$ problem with hinfstruct
opts = hinfstructOptions('Display','final','MaxIter',100,...
    'RandomStart',50);%,'TolGain',1e-7);
[CL,gam1] = hinfstruct(CL0,opts);  % CL is tuned version of CL0

if flag == 1
    % % Classic P-K structure
    % % ---------------------------------------------------
    % Get proportional and integral gains
    kp = CL.Blocks.K_norolloff.Kp.Value;
    ki = CL.Blocks.K_norolloff.Ki.Value;
    kd = CL.Blocks.K_norolloff.Kd.Value;
    tau = CL.Blocks.K_norolloff.Tf.Value;
    % Form the final controller
    K_norolloff = kp + kd*(s/(tau*s+1)) + ki/s;
    K = K_norolloff;
    % series(K_norolloff,(a/(s+a))^2);
    
elseif flag == 2
    % Hierarchical inner-outer structure
    % ---------------------------------------------------
    % Ko:
    % Get proportional and integral gains
    kop = CL.Blocks.Ko_norolloff.Kp.Value;
    koi = CL.Blocks.Ko_norolloff.Ki.Value;
    % Form the final controller
    Ko_norolloff = kop + koi/s;
    Ko = Ko_norolloff;
%     Ko = series(Ko_norolloff,(a/(s+a))^2);
    
    % Ki:
    % Get proportional and differential gains
    kip = CL.Blocks.Ki_norolloff.Kp.Value;
    kid = CL.Blocks.Ki_norolloff.Kd.Value;
    tau = CL.Blocks.Ki_norolloff.Tf.Value;
    kii = CL.Blocks.Ki_norolloff.Ki.Value;
    % Form the final controller
    Ki_norolloff = kip + kid*(s/(tau*s+1)) + kii/s;
    Ki = Ki_norolloff;
%     Ki = series(Ki_norolloff,(a/(s+a))^2);
    
end

%% Analyze OL and CL maps

if flag == 1
    % % Classic P-K structure
    % % ---------------------------------------------------
    [Lo,Li,So,Si,To,Ti,KS,PS] = f_CLTFM(P,K);
    S = So;
    T = To;
    zpk(T)
    
    % Plot S and T
    wvec = logspace(-4,3,1000);
    figure; sigma(So,wvec); grid on; hold on; sigma(To,wvec);
    plot_axis;
    [hL,hObj]=legend('S','T');
    plot_legend(hL,hObj);
    
elseif flag == 2
    % % Hierarchical inner-outer structure
    % % ---------------------------------------------------
    P.InputName = 'u';
    P.OutputName = 'y';
    Ko.InputName = 'e';
    Ko.OutputName = 'uo';
    Ki.InputName = 'y';
    Ki.OutputName = 'ui';
    sum_outer=sumblk('e=r-y',1);
    sum_inner=sumblk('u=uo-ui',1);
    
    Pmod = connect(P,Ki,sum_inner,'uo','y');
    zpk(Pmod)
    
    S = connect(P,Ki,Ko,sum_outer,sum_inner,'r','e');
    T = connect(P,Ki,Ko,sum_outer,sum_inner,'r','y');
    zpk(T)
    
    % Plot S and T
    wvec = logspace(-4,2,1000);
    figure; sigma(S,wvec); grid on; hold on; 
    sigma(T,wvec);
    sigma(inv(W1),wvec);
    sigma(inv(W3),wvec);
    title('Sensitivity and Complementary Sensitivity')
    plot_axis;
    [hL,hObj]=legend('S','T','W1^{-1}','W3^{-1}');
    plot_legend(hL,hObj);
end