function [value sg varargout] = f_Hinf(M, x, T11, T12, T21, Q, vec, varargin)
% Compute H-infinity norm and subgradients
% of Parameterized TFMs for given Youla et al. parameter Q

if nargin == 8
 conval = varargin{1};
 varargout{1} = conval;
end
n = length(x);

[n_u, n_e, n_s] = size(Q);
if isempty(T11)
    Twz = ss([]); 
else
    Twz = parallel(T11,series(series(T21,Q),T12));
end
%Twz = minreal(Twz);
Twz = Twz(vec,:);

% [ninf, fpeak] = norm(Twz, inf, 1e-8);
[ninf, fpeak] = hinfnorm(Twz, 1e-8);
value = ninf;

if fpeak < 1e-5
    fpeak = 1e-5; 
end

if fpeak<1e-5
    wmin=1e-8; wmax=1e-2;
elseif fpeak<1e-2
    wmin=1e-4; wmax=1e0;
elseif fpeak<1e0
    wmin=1e-1; wmax=1e1;
elseif fpeak<10
    wmin=1e-1; wmax=1e2;
elseif fpeak<1e2
    wmin=1e0; wmax=1e4;
elseif fpeak<1e5
    wmin=1e3; wmax=1e7;    
elseif fpeak>=1e5
    wmin=1e4; wmax=1e10;
else
    wmin=max([1, fpeak-10]); wmax=fpeak+10;
end
TwzScaled=prescale(Twz,{wmin,wmax});

Hjwo = freqresp(TwzScaled,fpeak);
% Hjwo = evalfr(TwzScaled,fpeak);

[U,S,V] = svd(Hjwo);	% SVD at W0
if isempty(U)
    uo = []; 
    vo = []; 
else
    uo 		= U(:,1); 		% Maximum Left Singular Vector
    vo 		= V(:,1); 		% Maximum Right Singular Vector 
end
subgradient = [];
for i = 1:n
    Hjwo = freqresp(M{i},fpeak);
    magHjwo = abs(Hjwo);
    subgradient = [subgradient; real(uo'*Hjwo*vo)];
end
sg = subgradient;
