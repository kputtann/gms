function [RGAMat,RGASumNorm]=f_RGADynSys(G,wvec)
% Code to find rga for a dynamic system
% Inputs: 
% G: Dynamic system, wvec: frequency sampling points
% Outputs:
% RGAMat: RGA matrix components for all specified frequencies, RGASumNorm

% % % Method-1
% RGASumNorm=zeros(length(wvec),1);
% RGAMat=zeros(size(G,1),size(G,2),length(wvec));
% for ii=1:length(wvec)
%     freq=wvec(ii);
%     G_mat=evalfr(G,freq);
%     R=rga(G_mat);
%     RGAMat(:,:,ii)=R;
%     RGASumNorm(ii)=sum(sum(abs(R)));
% end

% % Method-2
% wvec=logspace(1,1,100);
Gss = ss(G);
Gpck = pck(Gss.a,Gss.b,Gss.c,Gss.d);
Gw=frsp(Gpck,wvec);
RGAMat=veval('.*',Gw,vpinv(vtp(Gw)));
RGASumNorm = sum(abs(RGAMat(1:length(wvec),1:2)),2);