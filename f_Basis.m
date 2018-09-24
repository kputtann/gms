function q = f_Basis(N, p, z, basis_type)
% Form the basis TFs for given basis parameters and type

q{1} = zpk([],[],1);
% q{1} = tf(1,1);
if basis_type == 1  % fixed pole low pass
    for k=2:N
        q{k} = zpk([],-p,p)^(k-1);
    end
elseif basis_type == 2 % fixed pole all pass
    for k=2:N
        q{k} = zpk(p,-p,-1)^(k-1);
    end
elseif basis_type == 3 % variable pole low pass
    for k=2:N
        q{k} = zpk([],-p*(k-1),p*(k-1));
    end
elseif basis_type == 4 % variable pole all pass
    for k=2:N
        q{k} = zpk(p*(k-1),-p*(k-1),-1);
    end
elseif basis_type == 5 % pole and zero
    for k=2:N 
        q{k} = zpk(z,-p,-1)^(k-1);
    end
elseif basis_type == 5 % Laguerre 
    for k=2:N
        q{k} = zpk([],-p,sqrt(2*p))*zpk(p,-p,1)^(k-1); 
    end
end