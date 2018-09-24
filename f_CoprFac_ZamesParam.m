function [T11rz, T12rz, T21rz,T11dz, T12dz, T21dz]=f_CoprFac_ZamesParam(P_ss,F,L,weights)
% Zames Coprime Parameterization
% Works stable plants
% Assumes zero initial controller Ko

[Ap, Bp, Cp, Dp] = ssdata(P_ss);

[n_e, n_u] = size(P_ss);
NumP=P_ss;
DenP=ss(eye(n_u));
NumK=ss(zeros(n_u,n_e));
DenK=ss(eye(n_e));

NumPt=P_ss;
DenPt=ss(eye(n_e));
NumKt=ss(zeros(n_u,n_e));
DenKt=ss(eye(n_u));

% Feedback transfer function matrices
SOut11=DenK*DenPt; SOut12=-NumP; SOut21=DenPt;
KSOut11=NumK*DenPt; KSOut12=DenP; KSOut21=DenPt;
TOut11=NumP*NumKt; TOut12=NumP; TOut21=DenPt;
SensIn11=DenP*DenKt; SensIn12=-DenP; SensIn21=NumPt;
SInP11=NumP*DenKt; SInP12=-NumP; SInP21=NumPt;
% TIn11=NumK; TIn12=DenP; TIn21=NumPt;
TIn11=DenP*NumKt*inv(DenPt)*NumPt; TIn12=DenP; TIn21=NumPt;

% Weights
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

% Parameterization
T11rz1=W1*SOut11; T12rz1=W1*SOut12; T21rz1=SOut21;
T11rz2=W2*KSOut11; T12rz2=W2*KSOut12; T21rz2=KSOut21;
T11rz3=W3*TOut11; T12rz3=W3*TOut12; T21rz3=TOut21;
T11dz1=Wd1*SensIn11; T12dz1=Wd1*SensIn12; T21dz1=SensIn21;
T11dz2=Wd2*SInP11; T12dz2=Wd2*SInP12; T21dz2=SInP21;
T11dz3=Wd3*TIn11; T12dz3=Wd3*TIn12; T21dz3=TIn21;

% Constraint tf parameterization
T11rz1c=[]; T12rz1c=[]; T21rz1c=[];
for ii=1:length(W1c)
    T11rz1c=W1c{ii}.tfm*SOut11; T12rz1c=W1c{ii}.tfm*SOut12; T21rz1c=SOut21;
end
% if isempty(W2c)
%     T11rz2c=[]; T12rz2c=[]; T21rz2c=[];
% else
%     T11rz2c=W2c{1}.tfm*KSOut11; T12rz2c=W2c{1}.tfm*KSOut12; T21rz2c=KSOut21;
%     T11rz2c=[T11rz2c; W2c{2}.tfm*KSOut11]; T12rz2c=[T12rz2c; W2c{2}.tfm*KSOut12];
% end
T11rz2c=[]; T12rz2c=[]; T21rz2c=[];
for ii=1:length(W2c)
    T11rz2c=[T11rz2c; W2c{ii}.tfm*KSOut11]; T12rz2c=[T12rz2c; W2c{ii}.tfm*KSOut12];
end
T11rz3c=[]; T12rz3c=[]; T21rz3c=[];
for ii=1:length(W3c)
    T11rz3c=W3c{ii}.tfm*TOut11; T12rz3c=W3c{ii}.tfm*TOut12; T21rz3c=TOut21;
end
T11dz1c=[]; T12dz1c=[]; T21dz1c=[];
for ii=1:length(Wd1c)
    T11dz1c=Wd1c{ii}.tfm*SensIn11; T12dz1c=Wd1c{ii}.tfm*SensIn12; T21dz1c=SensIn21;
end
T11dz2c=[]; T12dz2c=[]; T21dz2c=[];
for ii=1:length(Wd2c)
    T11dz2c=Wd2c{ii}.tfm*SInP11; T12dz2c=Wd2c{ii}.tfm*SInP12; T21dz2c=SInP21;
end
T11dz3c=[]; T12dz3c=[]; T21dz3c=[];
for ii=1:length(Wd3c)
    T11dz3c=Wd3c{ii}.tfm*TIn11; T12dz3c=Wd3c{ii}.tfm*TIn12; T21dz3c=TIn21;
end

% For Trz1 and Tdiz2
T11rz=[T11rz1; T11rz2; T11rz3; T11rz1c; T11rz2c; T11rz3c]; T12rz=[T12rz1; T12rz2; T12rz3; T12rz1c; T12rz2c; T12rz3c]; T21rz=T21rz1;
T11dz=[T11dz1; T11dz2; T11dz3; T11dz1c; T11dz2c; T11dz3c]; T12dz=[T12dz1; T12dz2; T12dz3; T12dz1c; T12dz2c; T12dz3c]; T21dz=T21dz1;