function QN = f_FormQN(x, qk, n_u, n_e, N)
% From Q_N at given point x using the bases qk

xtemp = reshape(x,n_u*n_e,N);
QN = zeros(n_u, n_e);
for i = 1:N
    X{i} = reshape(xtemp(:,i),n_u,n_e);
    
    % % Find temp = QN + X{i} * qk{i}
    % % Straight forward way.
    % temp = QN + X{i} * qk{i};
    % Alternative way. minreal in later works better when this is used,
    % i.e., order of QN found will be as expected.
    QNss = ss(QN);
    xqss = ss(series(qk{i},X{i}));
    A = blkdiag(QNss.a,xqss.a);
    B = [QNss.b;xqss.b];
    C = [QNss.c,xqss.c];
    D = QNss.d+xqss.d;
    temp = zpk(ss(A,B,C,D));
    
    QN = minreal(temp,1e-6); 
    
end
QN = ss(QN);
