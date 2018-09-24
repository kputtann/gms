% Illustrative Example: Weighted Sensitivities vs Tradeoff Parameter
% Using nonlinear optimization solver
% % Assumptions:
% P is dynamic 
% K is static
% 2X2 system
% K is of the form [k11 k12; k21 k22]

clear all;
close all;

% Transfer function variable
s = tf('s'); 

% Plant definition
P = 1/s*[10 9; 9 8];

% Weighting functions
% % W1 and Wd1
Eps=0.01; Ms=2; wb1=1; wb2=1;
W1 = [tf([1/Ms wb1], [1 wb1*Eps]) 0; 0 tf([1/Ms wb2], [1 wb2*Eps])];
Eps=0.01; Ms=2; wb1=2; wb2=2;
Wd1 = [tf([1/Ms wb1], [1 wb1*Eps]) 0; 0 tf([1/Ms wb2], [1 wb2*Eps])];
% W2 and Wd2
W2 = ss([0.5 0; 0 2]);
Wd2 = ss([1 0; 0 1]);

% Constraint values
cvalvec = [10; 1000; 0.0001; 10000];

% rho vector
rho_vec = 0:0.05:1; 

% Bounds on optimization variable
k_ub = 20; 
k_lb = -20; 
kvec_ub = k_ub*ones(4,1);
kvec_lb = k_lb*ones(4,1);

% Initial point for the first iteration
kvec_init = [-8; 9; 9; -10];

% Loop for all the different values of rho
for ind = 1:numel(rho_vec)
    rho = rho_vec(ind);
    
    % Call the fmincon solver function
    [kvec_soln,fval,flg] = fmincon(@(x) obj_func(x,P,W1,Wd1,rho),...
        kvec_init,[],[],[],[],kvec_lb,kvec_ub,...
        @(c) constr_func(c,P,W2,Wd2,cvalvec));
    
    % Controller parameters obtained
    k11 = kvec_soln(1);
    k12 = kvec_soln(2);
    k21 = kvec_soln(3);
    k22 = kvec_soln(4);
    K = [k11 k12; k21 k22]
    
    % Store the relevant control properties
    K_store{ind} = K;
    [Lo,Li,So,Si,To,Ti,KS,PS] = f_CLTFM(ss(P),ss(K));
    NormInf = mag2db([hinfnorm(So), hinfnorm(Si), hinfnorm(KS), ...
        hinfnorm(PS), hinfnorm(To), hinfnorm(Ti)])
    norm1_store(ind) = hinfnorm(W1*So)
    norm2_store(ind) = hinfnorm(Wd1*Si)
    So_BW_20 = max(getGainCrossover(So,0.1));
    To_BW_20 = max(getGainCrossover(To,0.1));
    constr_store(:,ind) = [hinfnorm(W2*KS); hinfnorm(Wd2*PS); ...
        So_BW_20; To_BW_20]
    
    % Initial point for next iteration
    kvec_init = kvec_soln;
    
end

figure; plot(rho_vec,mag2db(norm1_store),'-*'); grid on; 
hold on;
plot(rho_vec,mag2db(norm2_store),'-*')
plot_axis;
title('||W*S||_{\infty} vs \rho');
[hL,hObj]=legend('||W_1*S_e||_{\infty}','||W_4*S_c||_{\infty}');
plot_legend(hL,hObj)
ylabel('||W*S||_{\infty} (dB)');
xlabel('\rho');

figure; plot(mag2db(norm1_store),mag2db(norm2_store),'-*'); grid on; 
plot_axis;
ylabel('||W_4*S_c||_{\infty} (dB)');
xlabel('||W_1*S_e||_{\infty} (dB)');
title('||W_4*S_c||_{\infty} vs ||W_1*S_e||_{\infty}');

%% 
function fval = obj_func(kvec,P,W1,Wd1,rho)
% % Objective Function
% Inputs: kvec P W1 Wd1 rho
% Output: fval

% Form the controller
k11 = kvec(1); 
k12 = kvec(2); 
k21 = kvec(3); 
k22 = kvec(4); 
K = [k11 k12; k21 k22];

% Compute the closed loop maps
[~,~,So,Si,~,~,~,~] = f_CLTFM(ss(P),ss(K));

% Objective Function
% fval = rho*hinfnorm(W1*So) + (1-rho)*hinfnorm(Wd1*Si);
fval = max(rho*hinfnorm(W1*So),(1-rho)*hinfnorm(Wd1*Si));
end

function [cineq,ceq] = constr_func(kvec,P,W2,Wd2,cvalvec)
% % Constraint Function
% Inputs: kvec P W2 Wd2
% Output: cineq ceq

% Form the controller
k11 = kvec(1); 
k12 = kvec(2); 
k21 = kvec(3); 
k22 = kvec(4); 
K = [k11 k12; k21 k22];

% Compute the closed loop maps
[~,~,So,Si,To,Ti,KS,PS] = f_CLTFM(ss(P),ss(K));

% Check the bandwidths 
if ~isempty(getGainCrossover(So,0.1))
    So_BW_20 = max(getGainCrossover(So,0.1)); 
else
    So_BW_20 = Inf;
end;
if ~isempty(getGainCrossover(To,0.1))
    To_BW_20 = max(getGainCrossover(To,0.1)); 
else
    To_BW_20 = Inf;
end;

% Inequality Constraints
cineq = [hinfnorm(W2*KS) - cvalvec(1); 
    hinfnorm(Wd2*PS) - cvalvec(2); 
    -So_BW_20 + cvalvec(3); 
    To_BW_20 - cvalvec(4)
    ];

% No Equality Constraint
ceq = [];
end
