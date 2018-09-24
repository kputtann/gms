function [n_e, n_u, DATArz,DATAdz] = f_GenData(P, weights)
% Extract data from problem setup for the GMS methodology

W1 = weights.W1; 
W2 = weights.W2; 
W3 = weights.W3; 
Wd1 = weights.Wd1; 
Wd2 = weights.Wd2; 
Wd3 = weights.Wd3; 

W1c = weights.W1c; 
W2c = weights.W2c; 
W3c = weights.W3c; 
Wd1c = weights.Wd1c; 
Wd2c = weights.Wd2c; 
Wd3c = weights.Wd3c; 

if (isfield(weights,'Wni1'))
    Wni1 = weights.Wni1; 
end
if (isfield(weights,'Wni3'))
    Wni3 = weights.Wni3; 
end

[n_e, n_u, n_s] = size(P);

nObj = 0;
%% Check W1
if ~isempty(W1)
    [noutput, ninput, nstate] = size(W1);
    if noutput ~= ninput
        disp('Error: W1 is not square')
        return
    end
    if noutput ~= n_e
        disp('Error: Dimansion mismatch in W1')
        return
    end    
    nObj = nObj+n_e;
end

%% Check W2
if ~isempty(W2)
    [noutput, ninput, nstate] = size(W2);
    if noutput ~= ninput
        disp('Error: W2 is not square')
        return
    end
    if noutput ~= n_u
        disp('Error: Dimansion mismatch in W2')
        return
    end    
    nObj = nObj+n_u;
end

%% Check W3
if ~isempty(W3)
    [noutput, ninput, nstate] = size(W3);
    if noutput ~= ninput
        disp('Error: W3 is not square')
        return
    end
    if noutput ~= n_e
        disp('Error: Dimansion mismatch in W3')
        return
    end    
    nObj = nObj+n_e;
end

DATArz.ObjVec = 1:nObj;
TotalRows = nObj;
%% rz
ConstraintCounter = 0;
[nRow nCol]=size(W1c);
for i=1:nCol
    W1 = W1c{i}.tfm;
    if ~isempty(W1)
        [noutput, ninput, nstate] = size(W1);
        if noutput ~= ninput
            disp(['Error: W1c{' num2str(i) '} is not square'])
            return
        end
        if noutput ~= n_e
            disp(['Error: Dimansion mismatch in W1c{' num2str(i) '}'])
             return
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATArz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_e;
        DATArz.ConNam{ConstraintCounter} = W1c{i}.Fun;
        DATArz.ConVal{ConstraintCounter} = W1c{i}.Val;
        TotalRows = TotalRows + n_e;
    end
end

%%
[nRow nCol]=size(W2c);
for i=1:nCol
    W2 = W2c{i}.tfm;
    if ~isempty(W2)
        [noutput, ninput, nstate] = size(W2);
        if noutput ~= ninput
            disp(['Error: W2c{' num2str(i) '} is not square'])
            return
        end
        if noutput ~= n_u
            disp(['Error: Dimansion mismatch in W2c{' num2str(i) '}'])
             return
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATArz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_u;
        DATArz.ConNam{ConstraintCounter} = W2c{i}.Fun;
        DATArz.ConVal{ConstraintCounter} = W2c{i}.Val;
        TotalRows = TotalRows + n_u;
    end
end

%%
[nRow nCol]=size(W3c);
for i=1:nCol
    W3 = W3c{i}.tfm;
    if ~isempty(W3)
        [noutput, ninput, nstate] = size(W3);
        if noutput ~= ninput
            disp(['Error: W3c{' num2str(i) '} is not square'])
            return
        end
        if noutput ~= n_e
            disp(['Error: Dimansion mismatch in W3c{' num2str(i) '}'])
             return
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATArz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_e;
        DATArz.ConNam{ConstraintCounter} = W3c{i}.Fun;
        DATArz.ConVal{ConstraintCounter} = W3c{i}.Val;
        TotalRows = TotalRows + n_e;
    end
end
DATArz.ConNum = ConstraintCounter;


%% dz
nObj = 0;
%% Check Wd1
if ~isempty(Wd1)
    [noutput, ninput, nstate] = size(Wd1);
    if noutput ~= ninput
        disp('Error: Wd1 is not square')
        return
    end
    if noutput ~= n_u
        disp('Error: Dimansion mismatch in Wd1')
        return
    end    
    nObj = nObj+n_u;
end

%% Check Wd2
if ~isempty(Wd2)
    [noutput, ninput, nstate] = size(Wd2);
    if noutput ~= ninput
        disp('Error: Wd2 is not square')
        return
    end
    if noutput ~= n_e
        disp('Error: Dimansion mismatch in Wd2')
        return
    end    
    nObj = nObj+n_e;
end

%% Check Wd3
if ~isempty(Wd3)
    [noutput, ninput, nstate] = size(Wd3);
    if noutput ~= ninput
        disp('Error: Wd3 is not square')
        return
    end
    if noutput ~= n_u
        disp('Error: Dimansion mismatch in Wd3')
        return
    end    
    nObj = nObj+n_u;
end

DATAdz.ObjVec = 1:nObj;
TotalRows = nObj;
%%
ConstraintCounter = 0;

[nRow nCol]=size(Wd1c);
for i=1:nCol
    Wd1 = Wd1c{i}.tfm;
    if ~isempty(Wd1)
        [noutput, ninput, nstate] = size(Wd1);
        if noutput ~= ninput
            disp(['Error: Wd1c{' num2str(i) '} is not square'])
            return
        end
        if noutput ~= n_e
            disp(['Error: Dimansion mismatch in Wd1c{' num2str(i) '}'])
             return
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAdz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_e;
        DATAdz.ConNam{ConstraintCounter} = Wd1c{i}.Fun;
        DATAdz.ConVal{ConstraintCounter} = Wd1c{i}.Val;
        TotalRows = TotalRows + n_e;
    end
end

%%
[nRow nCol]=size(Wd2c);
for i=1:nCol
    Wd2 = Wd2c{i}.tfm;
    if ~isempty(Wd2)
        [noutput, ninput, nstate] = size(Wd2);
        if noutput ~= ninput
            disp(['Error: Wd2c{' num2str(i) '} is not square'])
            return
        end
        if noutput ~= n_u
            disp(['Error: Dimansion mismatch in Wd2c{' num2str(i) '}'])
             return
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAdz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_u;
        DATAdz.ConNam{ConstraintCounter} = Wd2c{i}.Fun;
        DATAdz.ConVal{ConstraintCounter} = Wd2c{i}.Val;
        TotalRows = TotalRows + n_u;
    end
end

%%
[nRow nCol]=size(Wd3c);
for i=1:nCol
    Wd3 = Wd3c{i}.tfm;
    if ~isempty(Wd3)
        [noutput, ninput, nstate] = size(Wd3);
        if noutput ~= ninput
            disp(['Error: Wd3c{' num2str(i) '} is not square'])
            return
        end
        if noutput ~= n_e
            disp(['Error: Dimansion mismatch in Wd3c{' num2str(i) '}'])
             return
        end
        ConstraintCounter = ConstraintCounter + 1;
        DATAdz.ConVec{ConstraintCounter} = TotalRows+1:TotalRows+n_e;
        DATAdz.ConNam{ConstraintCounter} = Wd3c{i}.Fun;
        DATAdz.ConVal{ConstraintCounter} = Wd3c{i}.Val;
        TotalRows = TotalRows + n_e;
    end
end


DATAdz.ConNum = ConstraintCounter;
