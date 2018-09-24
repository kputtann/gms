function [condnum_store,L_opt_store,R_opt_store] = f_MinCondNum(...
    ss_Plant,wvec,gamma_vec)

% Function to obtain the minimized condition number of a dynamical system
% Min cond num at each desired frequency value is obtained
% The value of optimization variable gamma is guessed each time before
% successively solving an LMI. The gamma is looped over several values. 
% The minimum condition number obtained over all the gamma values is chosen
% Set wvec, the frequency points where condnum need to be mimimized
% Set gamma_vec, the values of guess for optimization variable gamma at
% each frequency
% 
% Inputs: 
% 1) System whose minimized condition number is to be found
% 2) wvec: the frequency points where condnum need to be mimimized
% 3) gamma_vec, the values of guess for optimization variable gamma at each
% frequency
% Example: 
% wvec=logspace(-3,1,25);
% gamma_vec = logspace(-3,1,50);
% % load FlexNom.mat
% % load Linr-Bolender_NewEng_OldPlm1.mat
% % ss_Plant = HSV_Trim_Data{1,1}.FER_Scaled;
% load('NENP_Rigid');
% ss_Plant=plantRigid(1:2,:);
% % ss_Plant=[1/(s+1) 0; 0 1/(s+2)]*[9 -10; -8 9];
% % ss_Plant=[(10-s)/10 0; 0 1/(s+1)];
% 
% Outputs:
% 1) Cell of length same as that of wvec, each with diagonal matrix matrix 
% L of size equal to number of system outputs
% 2) Cell of length same as that of wvec, each with diagonal matrix matrix 
% R of size equal to number of system inputs
% 3) Vector of Minimized Condition Number at every frequency point given by
% wvec

% addpath(genpath('path to yalmip'))
% addpath(genpath('path to SeDuMi'))

% clear;
% warning off;
% s=tf('s');

[n_y,n_u] = size(ss_Plant); 

% wvec=1e-1;
% Preallocate:
condnum_store = NaN*ones(1,length(wvec));
L_opt_store{length(wvec)} = [];
R_opt_store{length(wvec)} = [];

for j=1:length(wvec)
    
    M = bode(ss_Plant,wvec(j));
    eps = 1e-4;
    
    for jj = 1:length(gamma_vec) 
        g = gamma_vec(jj); % this is the optimization objective
        tic;
        if cond(M) > 1
            P = diag(sdpvar(n_y,1));
            Q = diag(sdpvar(n_u,1));
            MID = M'*P*M;
            C1 = [ Q <= MID, MID <= g^2*Q, Q>eps*eye(2), P>eps*eye(2)];
            % C1 = [Q <= MID, MID <= gamma^2*Q, Q>eps*eye(2)];
            
            % sdpsettings('bmibnb.roottight',[0|1])
            sdpset = sdpsettings('solver','sedumi','verbose',1); ...
                %,'showprogress',1);
            
            % solvesdp(C1,g^2,sdpset);
            diagnostics = optimize(C1,g^2,sdpset);
            
            P = double(P);
            Q = double(Q);
            
            PT = isnan(P);
            QT = isnan(Q);
            PT1 = find(PT == 1);
            QT1 = find(QT == 1);
            PT2 = sum(PT1);
            QT2 = sum(QT1);
            if PT2 ==0 && QT2 ==0
                L = sqrtm(P);
                R = Q^(-0.5);
                condnum(jj) = cond(L*M*R);
                L_cur_gam{jj} = L; 
                R_cur_gam{jj} = R; 
            else
                condnum(jj) = Inf;
                L_cur_gam{jj} = [];
                R_cur_gam{jj} = [];
            end
            clear PT QT PT1 QT1 PT2 QT2
        else
            condnum = 1;
        end
        
        toc;
    end
    
    [condnum_store(j),ind] = min(condnum);
    clear condnum;
    L_opt_store{j}=L_cur_gam{ind}; 
    R_opt_store{j}=R_cur_gam{ind}; 
end

