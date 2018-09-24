function [x, fx, iter_cnt, perf_meas] = f_ACCPM_GenMixSens_Optimizer(N,NQ,x0,Mobjrz,Mobjdz,Mconrz,Mcondz,T11rz, T12rz, T21rz,T11dz, T12dz, T21dz, Datarz, Datadz, Q,q,n_u,n_e,xmax,xmin,MaxIter)
% ACCPM main file

% Inputs
% xk is initial x

%% User initialisation of the ACCPM parameters
% x0 = x01*ones(N,1);
[problemS,methodS] = UserInitACCPM(N, x0,xmax,xmin,MaxIter);

%% Call the initialization routine
[accpm, accpm2Oracle] = C_InitProxAccpm(problemS, methodS);
clear problemS methodS;

%% Optimization process 
Flag = 1;
% Store the num of iterations and objective func value at each iter
iter_cnt = 0;
perf_meas = NaN*ones(MaxIter,1);
while(Flag)
  % Function evaluate at current point (y) using the Oracle
  [oracleS] = UserOracle(accpm2Oracle, Mobjrz,Mobjdz, Mconrz,Mcondz, T11rz, T12rz, T21rz,T11dz, T12dz, T21dz, Q, Datarz,Datadz);
  % Call the query point generator to get the next point
  [accpm, accpm2Oracle] = ...
  C_ProxAccpmGen(oracleS, accpm, accpm2Oracle);
  x = get(accpm2Oracle,'y');
  Q = f_FormQN(x, q, n_u, n_e, NQ); 
  clear oracleS;
  % Possibly artificial stop
  condition = (get(accpm2Oracle,'ExitCode') ~= 0);
  if (condition)
    Flag = 0;
  end
  iter_cnt = iter_cnt+1;
  perf_meas(iter_cnt) = accpm.ParamS.OptTypeFact * accpm.ManagerS.D.ObjBounds(2);
end
% Minimizer and Objective
x  = get(accpm2Oracle,'y');
fx = accpm.ParamS.OptTypeFact * accpm.ManagerS.D.ObjBounds(2);
clear accpm2Oracle;
% Display Results
C_Display_ProxAccpmResults(accpm);
clear accpm;
return

function [problemS,methodS] = UserInitACCPM(n, x0,xmax,xmin,MaxIter)
% In this function, the user initializes the ACCPM parameters 

% problemS object created
problemS = ProblemS;
problemS = set(problemS,'OptType','min');   % min or max            
problemS = set(problemS,'NbVariables',n);   % Number of dual variables
problemS = set(problemS,'NbSubProblems',1); % Number of subproblems
problemS = set(problemS,'StartingPoint',x0);  % starting point 
problemS = set(problemS,'VarLowerBounds',xmin); % Lower bounds on variables
problemS = set(problemS,'VarUpperBounds',xmax); % Upper bounds on variables

% methodS object created
methodS = MethodS;
methodS = set(methodS,'Proximal',0); % Use the proximal term (0/1)
methodS = set(methodS,'Rho',1);      % Value of the rho in the proximal term
methodS = set(methodS,'Verbose',3);  % Display results (0/1/2/3)
methodS = set(methodS,'MaxOuter',MaxIter); % Maximum number of ACCPM iterations
methodS = set(methodS,'MaxInner',100); % Maximum number of Newton iterations in the computation of the analytic center
methodS = set(methodS,'Tolerance',1e-4); % Relative optimality gap
methodS = set(methodS,'WeightEpigraphCutInit',10); % Initial weight on the epigraph cut
methodS = set(methodS,'WeightEpigraphCutInc',0); % Increment on the epigraph cut

return


function  [oracleS] = UserOracle(accpm2Oracle, Mobjrz,Mobjdz, Mconrz,Mcondz, T11rz, T12rz, T21rz,T11dz, T12dz, T21dz, Q, Datarz,Datadz)
% User oracle

% Current point
x = get(accpm2Oracle,'y');

oracleS = OracleS; 

% Constraint
% Reference to output
for ii = 1:Datarz.ConNum
  Mrz = Mconrz(ii,:);  
  [fc,Gfc,val] = ...
      feval(Datarz.ConNam{ii}, Mrz, x, T11rz, T12rz, T21rz, Q,...
                     Datarz.ConVec{ii}, Datarz.ConVal{ii});
  if (fc > val)
    oracleS = set(oracleS, 'FunctionValues', fc-val); % Value of the objective
    oracleS = set(oracleS, 'SubGradients', Gfc);  % Subgradient
    oracleS = set(oracleS, 'SubProblemIndex', 0); % Nature of the cut (Optimality -> 1, Feasibility -> 0)
    return  
  end  
end
% d_i to output
for ii = 1:Datadz.ConNum
  Mdz = Mcondz(ii,:);  
  [fc,Gfc,val] = ...
      feval(Datadz.ConNam{ii}, Mdz, x, T11dz, T12dz, T21dz, Q,...
                     Datadz.ConVec{ii}, Datadz.ConVal{ii});
  if (fc > val)
    oracleS = set(oracleS, 'FunctionValues', fc-val); 
    oracleS = set(oracleS, 'SubGradients', Gfc); 
    oracleS = set(oracleS, 'SubProblemIndex', 0);
    return  
  end  
end

% Optimality
[forz,Gforz] = f_Hinf(Mobjrz, x, T11rz, T12rz, T21rz, Q, Datarz.ObjVec);  
[fodz,Gfodz] = f_Hinf(Mobjdz, x, T11dz, T12dz, T21dz, Q, Datadz.ObjVec);  
if forz>=fodz
    oracleS = set(oracleS, 'FunctionValues', forz); % Value of the objective
    oracleS = set(oracleS, 'SubGradients', Gforz);  % Subgradient
    oracleS = set(oracleS, 'SubProblemIndex', 1); % Nature of the cut (Optimality -> 1, Feasibility -> 0)
else
    oracleS = set(oracleS, 'FunctionValues', fodz); 
    oracleS = set(oracleS, 'SubGradients', Gfodz);  
    oracleS = set(oracleS, 'SubProblemIndex', 1); 
end
return