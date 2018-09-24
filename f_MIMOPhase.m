function [phase_vec] = f_MIMOPhase(Plnt,wval_vec)
% Phase of SISO/MIMO Square Plant
%
% Inputs:
%   Plnt: Dynamical system (usually tf/zpk/ss)
%   wvec: [Optional] Frequency vector for obtaining phase at each
%         of those points
%         [Default] wvec = logspace(-3,3,1e3)
% Outputs:
%   phase_vec: Computed phase of system (in deg)
%
% For SISO systems, for getting the phase same as usual Bode phase
% plot, pick wi=-ui
% See Jie Chen, Multivar Gain-Phase ... , 1998
% % Also see Freudenberg Book
% Here, it is assumed that the reference vector wi=ui.

% Frequency vector: Assign Default if not provided
if nargin < 2
    wval_vec = logspace(-3,3,1e3);
end

% System
% s = tf('s');
%
% Plnt = [1/(s+1) 0; 0 1/(s+2)] * [9 -10; -8 9];
% % Plnt = [1/(s+1) 0; 0 1/(s+1)];

% sval = 1i*0.5;

for sval_ind = 1:length(wval_vec)
    
    sval = wval_vec(sval_ind);
    Plnt_jw = evalfr(Plnt,1i*sval);
    
    [u,sv,v] = svd(Plnt_jw);
    
    for singval_ind = 1:size(Plnt,1)
        
        phase_sv(sval_ind,singval_ind) = angle(u(:,singval_ind)'...
            *v(:,singval_ind))*180/pi;
        % fprintf('Angle (deg) corresponding to sv1')
        % angle(u(:,1)'*v(:,1))*180/pi
        
        % phase_sv2(sval_ind) = angle(u(:,2)'*v(:,2))*180/pi;
        % % fprintf('Angle (deg) corresponding to sv2')
        % % angle(u(:,2)'*v(:,2))*180/pi
    end
        
end

figure; 
for singval_ind = 1:size(Plnt,1)
    semilogx(wval_vec,phase_sv(:,singval_ind));
    hold on; 
end
title('\angle u_i^H v_i');
ylabel('(deg)');
xlabel('Frequency (rad/s)');
% legend('sv1','sv2');
plot_axis;
% figure;
% semilogx(wval_vec,phase_sv1,wval_vec,phase_sv2);
% title('\angle u_i^H v_i');
% ylabel('(deg)');
% xlabel('Frequency (rad/s)');
% legend('sv1','sv2');
% plot_axis;
