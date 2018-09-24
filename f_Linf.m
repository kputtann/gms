function [value,sg,ConValVec,varargout] = f_Linf(M,x,T11,T12,...
    T21, Q, vec, varargin)
% Compute peak value and subgradients
% of Parameterized TFMs for given Youla et al. parameter Q

tvec = 0:0.001:10;
if nargin == 8
    conval = varargin{1};
end
n = length(x);

Twz = parallel(T11,series(series(T21,Q),T12));
Twz = Twz(vec,:);
% [n_output,n_input] = size(Twz); % n_row = n_output

% subgradient = NaN*zeros(n,length(conval));
Counter=0;
for ii = 1:size(conval,1)
    for jj = 1:size(conval,2)
        %         kk=(ii-1)*size(conval,2)+jj;
        if conval(ii,jj)==Inf
            disp('');
            
        else
            Counter=Counter+1;
            ConValVec(Counter)=conval(ii,jj);
            
            [y,tvec] = step(Twz(ii,jj), tvec);
            
            [ypeak,I] = max(y);
            tpeak = tvec(I);
            value(Counter,1) = ypeak;
            
            for i = 1:n
                %             if n_output == n_input
                %                 [y,tvec] = step(M{i}(ii,jj), tvec);
                %             else
                %                 [y,tvec] = step(M{i}(ii,1), tvec);
                %             end
                [y,tvec] = step(M{i}(ii,jj), tvec);
                subgradient(i,Counter) = y(I);
            end
            
            
            if nargin == 8
                varargout{1} = conval(ii);
            end
            %     if value > conval(ii) % See why this is required
            %         return
            %     end
        end
    end
end
sg = subgradient;