function [xk,frz,fdz]=f_KelleyCPM_GenMix_Optimizer(N,NQ,x0,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz, T12rz, T21rz,T11dz, T12dz, T21dz, Datarz, Datadz, Q,q,n_u,n_e,MaxIter,xmax,xmin)
% Kelley's CPM

warning off
tol_obj = 1e-4;
tol_feas = 1e-4;


% INITIALIZE
% fx = 0;               % Set output to zero
iter = 0;               % Iteration count
xk  = x0;               % Initial query point
xkStore=NaN*ones(length(xk),MaxIter);
ExitFlagStore=NaN*ones(1,MaxIter);
foStore=NaN*ones(2,MaxIter);

nConrz = Datarz.ConNum; nCondz = Datadz.ConNum; % Number of constraints
% Below matrices are used in solving the LP: min c'x   s.t. Aw<b
Ao = [];                % A matrix associated with objective function
bo = [];                % b vector associated with objective function

c = [zeros(N,1); 1];    % cvector associated with the variable x
UkminLkrz=1000; UkminLkdz=1000; constraint_flagrz=1; constraint_flagdz=1;
w=zeros(N+1,1);

options = optimset('Display','off'); 
% options=optimset('MaxIter',500,'TolFun',1e-9,'Display','final')

% START
while UkminLkrz > tol_obj || UkminLkdz > tol_obj || (constraint_flagrz>0) || (constraint_flagdz>0)
    
    Ac=[]; bc=[];
    [forz, Gfo] =feval('f_Hinf', Mobjrz, xk, T11rz, T12rz, T21rz, Q, Datarz.ObjVec);
    if UkminLkrz > tol_obj
        Ao = [Ao; Gfo' -1];
        bo = [bo; Gfo'*xk-forz];
        %     UkminLkrz3=for3-c'*w;
    end
    
%         Constraints rz:
    % Compute fi(x), Gfi(x) and Form Ac, bc
    frz{1}=[];
    for ii = 1:nConrz
        Mrz = Mconrz(ii,:);
        [frz{ii}, Gf{ii}, ConValVec] = ...
            feval(Datarz.ConNam{ii}, Mrz, xk, T11rz, T12rz, T21rz, Q, Datarz.ConVec{ii}, Datarz.ConVal{ii});
        frz{ii} = frz{ii} - ConValVec'; 
        if constraint_flagrz>0
            Ac = [Ac; Gf{ii}' zeros(size(Gf{ii}',1),1)];
            bc = [bc; Gf{ii}'*xk-frz{ii}];
        end
    end

    [fodz, Gfo] =feval('f_Hinf', Mobjdz, xk, T11dz, T12dz, T21dz, Q, Datadz.ObjVec);
    if UkminLkdz > tol_obj
        Ao = [Ao; Gfo' -1];
        bo = [bo; Gfo'*xk-fodz];
        %     UkminLkdiz3=fo3-c'*w;
    end
    
    %     Constraints dz:
    % Compute fi(x), Gfi(x) and Form Ac, bc
    fdz{1}=[];
    for ii = 1:nCondz
        Mdz = Mcondz(ii,:);
        [fdz{ii}, Gf{ii}, ConValVec] = ...
            feval(Datadz.ConNam{ii}, Mdz, xk, T11dz, T12dz, T21dz, Q, Datadz.ConVec{ii}, Datadz.ConVal{ii});
        fdz{ii} = fdz{ii} - ConValVec';  % In AllStep this is changed. Originally it was frz{ii} - Datarz.ConVal{ii}
        if constraint_flagdz>0
            Ac = [Ac; Gf{ii}' zeros(size(Gf{ii}',1),1)];
            bc = [bc; Gf{ii}'*xk-fdz{ii}];
        end
    end

    Ao=[Ao;Ac]; bo=[bo;bc];
    
    % Solve LP (used optimization toolbox function: linprog)
    [w,fval,exitflag] = linprog(c,Ao,bo,[],[],xmin,xmax,xk,options);
    % Check if problem is giving empty w. Try using other algorithms
    optionstemp=options; % temporary option
    if exitflag == -4
        optionstemp.Algorithm='dual-simplex';
        [w,fval,exitflag] = linprog(c,Ao,bo,[],[],xmin,xmax,xk,optionstemp);
    end
    if exitflag == -4
        optionstemp.Algorithm='active-set';
        [w,fval,exitflag] = linprog(c,Ao,bo,[],[],xmin,xmax,xk,optionstemp);
    end
    
    UkminLkrz=forz-c'*w; fprintf('\n%d %1.6f  %1.6f  ', iter,forz, UkminLkrz);

    UkminLkdz=fodz-c'*w; fprintf('%1.6f  %1.6f  ', fodz, UkminLkdz);

    foStore(:,iter+1)=[forz; fodz];
    
    % Update xk
    xk = w(1:N); xkStore(:,iter+1)=xk;ExitFlagStore(1,iter+1)=exitflag;
    
    % Increment iter
    iter = iter + 1;
    % Check if fi(xk) < epsilon for all i
    constraint_flagrz = 0;
    for ii = 1:nConrz
        if frz{ii} > tol_feas
            constraint_flagrz = 1;
        end
    end
    constraint_flagdz = 0;
    for ii = 1:nCondz
        if fdz{ii} > tol_feas
            constraint_flagdz = 1;
        end
    end
    
    if iter == MaxIter
        fprintf('\n');
        fprintf('Max Num of Iter exceeded \n')
        break;
    end
    Q = f_FormQN(xk, q, n_u, n_e, NQ); 
end