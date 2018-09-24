% Sensitivity Limits Imposed by RHP Zero: 
% Sensitivity bounds based on RHPZ and Weighting function
% SISO LTI plant with P-K classic feedback structure
% Stable plant

% -----------------------------------------------------------
clear; 
s = tf('s'); 

% Available Bandwidth
wp = 1e5; 

% RHP zero 
z = 10; 

% % % Generic Weighting FUnction: 
% W = ((s+ws*nthroot(M,k2))^k2 * (s+wp/nthroot(M,k3))^k3)/ ...
%     ((s+ws*nthroot(epsilon,k1))^k1 * (s+ws)^(k2-k1) * (s+wp)^k3);
    
% Order of transfer function in a given frequency range
% First slope:
k1 = 1; 
% Second slope:
k2 = 1;
% Third slope:
k3 = 1; 

% Performance Bandwidth Vector
ws_vec = logspace(-1,1,1000);
% Peak Sensitivity Value Vector
M_vec = logspace(0,2,1000);

% Epsilon (magnitude of sensitivity at low freq.)
% Assume that epsilon is 0 (or 0+ to be precise)

% --------------------------------------------------
% % Solve for M for different value of ws
% Expression involving M and ws
% z^k1 * (z+ws)^(k2-k1) * (z+wp)^k3 ...
%     - (z+ws*nthroot(M,k2))^k2 * (z+wp/nthroot(M,k3))^k3 == 0; 

% Preallocate
Msoln_vec = NaN*ones(numel(ws_vec),1); 
validM_vec = NaN*ones(numel(ws_vec),1); 

for ws_ind = 1:numel(ws_vec)
    ws = ws_vec(ws_ind); 
    
    
    M_param_func = @(M,ws,wp,z,k1,k2,k3) ...
        (z+ws*nthroot(M,k2))^k2 * (z+wp/nthroot(M,k3))^k3 ...
        - z^k1 * (z+ws)^(k2-k1) * (z+wp)^k3;
    % "Single" variable function:
    M_singlevar_fun = @(M) M_param_func(M,ws,wp,z,k1,k2,k3);
    % Initial Point:
    M_init = 1;
    % Solve for M (root of nonlinear expression):
    try
        M_soln = fzero(M_singlevar_fun,M_init);
    catch
        M_soln = NaN; %
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
ylim([-10 60])


% -------------------------------------------------------------
% Solve for ws for different value of M
% Expression involving ws:
% ws <= (pi*p + k3/nthroot(M,k3)*wp - k3*wp)/(-k1 - k2*nthroot(M,k2) + k2)

% Preallocate
ws_vec = NaN*ones(numel(ws_vec),1); 
validws_vec = NaN*ones(numel(ws_vec),1); 

for M_ind = 1:numel(M_vec)
    M = M_vec(M_ind); 
    
    ws_param_func = @(ws,M,wp,z,k1,k2,k3) ...
        (z+ws*nthroot(M,k2))^k2 * (z+wp/nthroot(M,k3))^k3 ...
        - z^k1 * (z+ws)^(k2-k1) * (z+wp)^k3;
    % "Single" variable function:
    ws_singlevar_fun = @(ws) ws_param_func(ws,M,wp,z,k1,k2,k3);
    % Initial Point:
    ws_init = 1;
    % Solve for ws (root of nonlinear expression):
    try 
        ws_soln = fzero(ws_singlevar_fun,ws_init);
    catch
        ws_soln = NaN; %
    end
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
ylim([1e-2 1e1])
