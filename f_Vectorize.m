function [M Mobj Mcon] = f_Vectorize(T11, T12, T21, qk, N, n_u, n_e, ProblemData)
% Vectorize Problem
% Forms M_{l} = M_{k}^{ij}
% l = (k-1)*nu*ne+(j-1)*nu+i;
% M_{k}^{ij} = T_{12}*B^{ij}*T_{21}*q_k
% T_wz = M_o + sum_{l=1}^{nu*ne*N} M_l x_l
Mobj = {};
Mcon = {};
Bij = zeros(n_u,n_e);
for k = 1:N
    for j = 1:n_e
        for i = 1:n_u
            l = (k-1)*n_u*n_e+(j-1)*n_u+i;
            Bij = zeros(n_u,n_e);
            Bij(i,j) = 1;
            [size_t21 temp] = size(T21.a);
            [temp size_t12] = size(T12.a);
            if isempty(T12)
                M{l} = ss([]); 
            else
                a = [T21.a zeros(size_t21,size_t12);
                    T12.b*Bij*T21.c T12.a];
                b = [T21.b; T12.b*Bij*T21.d];
                c = [T12.d*Bij*T21.c T12.c];
                d = T12.d*Bij*T21.d;
                M{l} = ss(a,b,c,d)*qk{k};
            end
        end
    end
end
for k = 1:N*n_e*n_u
    Mobj{k} = M{k}(ProblemData.ObjVec,:);
end
for i = 1:ProblemData.ConNum
    for k = 1:N*n_e*n_u
        Mcon{i,k} = M{k}(ProblemData.ConVec{i},:);
    end
end
