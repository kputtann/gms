function [T11rz, T12rz, T21rz, T11dz, T12dz, T21dz,T11niz,T12niz,T21niz]=f_CoprFac_hsvio_Tniu(P_ss,F,L,weights)
% Youla Coprime factorization
% For HSV inner-outer with three loop-breaking points

W1 = weights.W1; 
W2 = weights.W2; 
W3 = weights.W3; 
Wd1 = weights.Wd1; 
Wd2 = weights.Wd2; 
Wd3 = weights.Wd3; 
Wni1 = weights.Wni1; 
Wni3 = weights.Wni3; 

W1c = weights.W1c; 
W2c = weights.W2c; 
W3c = weights.W3c; 
Wd1c = weights.Wd1c; 
Wd2c = weights.Wd2c; 
Wd3c = weights.Wd3c; 

[Ap, Bp, Cp, Dp] = ssdata(P_ss);

% Right coprime factorization
NumP.a=Ap-Bp*F; NumP.b=Bp; NumP.c=Cp-Dp*F; NumP.d=Dp; NumP=ss(NumP.a,NumP.b,NumP.c,NumP.d);
DenP.a=Ap-Bp*F; DenP.b=Bp; DenP.c=-F; DenP.d=eye(size(DenP.c,1)); DenP=ss(DenP.a,DenP.b,DenP.c,DenP.d);
% Controller
NumK.a=Ap-Bp*F; NumK.b=-L; NumK.c=-F; NumK.d=zeros(size(NumK.c,1),size(NumK.b,2)); NumK=ss(NumK.a,NumK.b,NumK.c,NumK.d);
DenK.a=Ap-Bp*F; DenK.b=L; DenK.c=Cp-Dp*F; DenK.d=eye(size(DenK.c,1)); DenK=ss(DenK.a,DenK.b,DenK.c,DenK.d);
% Left coprime factorization
NumPt.a=Ap-L*Cp; NumPt.b=Bp-L*Dp; NumPt.c=Cp; NumPt.d=Dp; NumPt=ss(NumPt.a,NumPt.b,NumPt.c,NumPt.d);
DenPt.a=Ap-L*Cp; DenPt.b=-L; DenPt.c=Cp; DenPt.d=eye(size(DenPt.c,1)); DenPt=ss(DenPt.a,DenPt.b,DenPt.c,DenPt.d);
% Controller
NumKt.a=Ap-L*Cp; NumKt.b=-L; NumKt.c=-F; NumKt.d=zeros(size(NumKt.c,1),size(NumKt.b,2)); NumKt=ss(NumKt.a,NumKt.b,NumKt.c,NumKt.d);
DenKt.a=Ap-L*Cp; DenKt.b=-(Bp-L*Dp); DenKt.c=-F; DenKt.d=eye(size(DenKt.c,1)); DenKt=ss(DenKt.a,DenKt.b,DenKt.c,DenKt.d);

% Feedback transfer function matrices
SOut11=DenK*DenPt; SOut12=-NumP; SOut21=DenPt;
KSOut11=NumK*DenPt; KSOut12=DenP; KSOut21=DenPt;
TOut11=NumP*NumKt; TOut12=NumP; TOut21=DenPt;
SensIn11=DenP*DenKt; SensIn12=-DenP; SensIn21=NumPt;
SInP11=NumP*DenKt; SInP12=-NumP; SInP21=NumPt;
% TIn11=NumK*NumPt; TIn12=DenP; TIn21=NumPt; % Final!!!
% TIn11=DenP*NumKt*inv(DenPt)*NumPt; TIn12=DenP; TIn21=NumPt;
TIn11=eye(size(DenP,1))-DenP*DenKt; TIn12=DenP; TIn21=NumPt;

% Tniu11=NumK*DenPt; Tniu12=DenP; Tniu21=DenPt; % map from sensor noise in the inner-loop to control signal 
% Tniy11=NumK*DenPt; Tniy12=DenP; Tniy21=DenPt; % map from sensor noise in the inner-loop to output
% Tniei11=NumK*DenPt; Tniei12=DenP; Tniei21=DenPt; % map from sensor noise in the inner-loop to output

%%%%%%%%%%%%%%%%%%%%%%%
% Inner-Outer loop - Select required tf's
Tniu11=KSOut11(1:2,3); Tniu12=KSOut12(1:2,:); Tniu21=KSOut21(:,3);
% Tniy11=TOut11(1:2,3); Tniy12=TOut12(1:2,3); Tniy21=TOut21(:,3);
% Tnixp11=TOut11(1:2,3); Tnixp12=TOut12(1:2,3); Tnixp21=TOut21(:,3);
Tniei11=SOut11(3,3); Tniei12=SOut12(3,:); Tniei21=SOut21(:,3);

SOut11=SOut11(1:2,1:2); SOut12=SOut12(1:2,1:2); SOut21=SOut21(:,1:2);
KSOut11=KSOut11(1:2,1:2); KSOut12=KSOut12(1:2,1:2); KSOut21=KSOut21(:,1:2);
TOut11=TOut11(1:2,1:2); TOut12=TOut12(1:2,1:2); TOut21=TOut21(:,1:2);
SInP11=SInP11(1:2,1:2); SInP12=SInP12(1:2,1:2); SInP21=SInP21(:,1:2);
% SensIn21=NumPt(1:2,:); 
% TIn21=NumPt(1:2,:); 
%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added for Integ Augment fix
s=tf('s');
% T11rz1=series(W1*SOut11,ss(s+1e-6)); T12rz1=series(W1*SOut12,ss(s+1e-6)); T21rz1=SOut21;
T11rz2=series(W2*KSOut11,ss(1/(s+1e-1))); T12rz2=series(W2*KSOut12,ss(1/(s+1e-1))); T21rz2=KSOut21;
% % The eps in the integrator (1/(s+eps)) is picked to relatively high value
% % ~0.1 because of a bad effect caused by bilinear transformation. When
% % bilinear transformation is done, the weighted KS or
% % (T11rz2+T12rz2**T21rz2) was going to a high value at low frequencies
% % (near DC), even though the acutal weighted KS (i.e., without bilin) was
% % low at those frequencies. 
% T11dz1=series(Wd1*SensIn11,ss(s+1e-6)); T12dz1=series(Wd1*SensIn12,ss(s+1e-6)); T21dz1=SensIn21;
T11dz2=series(Wd2*SInP11,ss(s+1e-6)); T12dz2=series(Wd2*SInP12,ss(s+1e-6)); T21dz2=SInP21;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameterization

T11rz1=W1*SOut11; T12rz1=W1*SOut12; T21rz1=SOut21;
% T11rz2=W2*KSOut11; T12rz2=W2*KSOut12; T21rz2=KSOut21;
T11rz3=W3*TOut11; T12rz3=W3*TOut12; T21rz3=TOut21;
T11dz1=Wd1*SensIn11; T12dz1=Wd1*SensIn12; T21dz1=SensIn21;
% T11dz2=Wd2*SInP11; T12dz2=Wd2*SInP12; T21dz2=SInP21;
T11dz3=Wd3*TIn11; T12dz3=Wd3*TIn12; T21dz3=TIn21;

T11niz1=Wni1*Tniu11; T12niz1=Wni1*Tniu12; T21niz1=Tniu21;
% T11niz2=Wni2*Tnixp11; T12niz2=Wni1*Tnixp12; T21niz2=Tnixp21;
T11niz3=Wni3*Tniei11; T12niz3=Wni3*Tniei12; T21niz3=Tniei21;

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
T11niz=[T11niz1; T11niz3]; T12niz=[T12niz1; T12niz3]; T21niz=T21niz1;