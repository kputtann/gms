function cond_num=f_CondNum(G,wvec)
% Code to find condition number of a dynamic system
% Inputs: 
%   G: Dynamic system, 
%   wvec: frequency sampling points
% Outputs: 
%   Condition number vector 
sing_val=sigma(G,wvec);
cond_num=sing_val(1,:)./sing_val(end,:);

% % Plot
% figure; loglog(wvec,cond_num); grid on;