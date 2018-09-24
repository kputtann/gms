% Bode Sensitivity Integral
% SISO LTI plant with P-K classic feedback structure
% Open loop TF is rational and has at least 2-pole roll-off
clear; 
close all;
% Sensitivity Integral Relations for 
% Generic upper bound on sensitivity

% Available Bandwidth
wp = 10; 

% Performance Bandwidth Vector
ws_vec = logspace(-1,1,1000);
% Peak Sensitivity Value Vector
M_vec = logspace(0,2,1000); 

% RHP pole 
p = 0;  % p=0 => Stable Plant

% Order of transfer function in a given frequency range
% First slope:
k1 = 1; 
% Second slope:
k2 = 1;
% Third slope:
k3 = 1; 

% Epsilon (magnitude of sensitivity at low freq.)
% Assume that epsilon is 0 (or 0+ to be precise)

% --------------------------------------------------
% % Solve for M for different value of ws
% Expression involving M:
% k2*nthroot(M,k2)*nthroot(M,k3)*ws ...
% - nthroot(M,k3)*(-pi*p - k1*ws + k2*ws + k3*wp) + k3*wp <= 0

% Preallocate
Msoln_vec = NaN*ones(numel(ws_vec),1); 
validM_vec = NaN*ones(numel(ws_vec),1); 

for ws_ind = 1:numel(ws_vec)
    ws = ws_vec(ws_ind); 
    
    % Parameterized function:
    M_param_func = @(M,ws,wp,p,k1,k2,k3) ...
        k2*nthroot(M,k2)*nthroot(M,k3)*ws ...
        - nthroot(M,k3)*(-pi*p - k1*ws + k2*ws + k3*wp) + k3*wp;
    % "Single" variable function:
    M_singlevar_fun = @(M) M_param_func(M,ws,wp,p,k1,k2,k3);
    % Initial Point:
    M_init = 1;
    % Solve for M:
    try 
        M_soln = fzero(M_singlevar_fun,M_init);
    catch
        M_soln = NaN; 
    end
    % Store the solution:
    Msoln_vec(ws_ind) = M_soln; 
    
    % For the assumed upper bound on sensitivity, the relations are valid 
    % under certain assumption. This assumption is based on relation 
    % between ws and wp.
    validM_vec(ws_ind) = (wp/ws)^((k2*k3)/(k2+k3));  
end
% % Plot
figure; 
semilogx(ws_vec,mag2db(Msoln_vec),'-b');
grid on; hold on; 
semilogx(ws_vec,mag2db(validM_vec),'-r');
title('LB on M vs \omega_s');
ylabel('LB on M (dB)');
xlabel('\omega_s (rad/s)');
plot_axis;
ylim([0 30])

% --------------------------------------------------
% Solve for ws for different value of M
% Expression involving ws:
% ws <= (pi*p + k3/nthroot(M,k3)*wp - k3*wp)/(-k1 - k2*nthroot(M,k2) + k2)

% Preallocate
ws_vec = NaN*ones(numel(ws_vec),1); 
validws_vec = NaN*ones(numel(ws_vec),1); 

for M_ind = 1:numel(M_vec)
    M = M_vec(M_ind); 
    
    % Solve for ws
    ws_soln = (pi*p + k3/nthroot(M,k3)*wp - k3*wp)...
        /(-k1 - k2*nthroot(M,k2) + k2);
    % Store ws
    ws_vec(M_ind) = ws_soln;
    
    % For the assumed upper bound on sensitivity, the relations are valid 
    % under certain assumption. This assumption is based on relation 
    % between ws and wp.
    validws_vec(M_ind) = wp/(nthroot(M,k2)*nthroot(M,k3));  
end
% % Plot
figure; 
semilogy(mag2db(M_vec),ws_vec,'-b');
grid on; hold on; 
semilogy(mag2db(M_vec),validws_vec,'-r');
title('UP on \omega_s vs M');
xlabel('M (dB)');
ylabel('UB on Perf. Bandwidth \omega_s (rad/s)');
plot_axis;
ylim([1e-1 1e1])
