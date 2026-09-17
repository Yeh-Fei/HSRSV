% Aurthor: Ye Fei
% E-mail: feiye@nau.edu.cn
% Create Date: May 29, 2026
% MATLAB R2021a

clear;  clc;
warning('off');
addpath(genpath(pwd));

if isempty(gcp('nocreate'))        
    parpool(maxNumCompThreads); 
end

SNRm=30; SNRh=25;
kernel_type={'Gaussian','Uniform'};   
kernel_type=kernel_type{1};
sf=5;
Tab=table;
Tab.name={'PSNR';'RMSE';'ERGAS';'SAM';'UIQI';'SSIM';'DD';'CC';'Time'};
rng(2,'twister');
s=rng;

load('Synthesized.mat');
data.F=Zh_th;  data.Fm=Zm_th;
clear Zh_th Zm_th;
szF=size(data.F);

%%  LR-HSI
kernel_length=9;  shift=1;
if strcmp(kernel_type,'Gaussian')
    sig=(1/(2*(2.7725887)/sf^2))^0.5;
    par.psf=fspecial('gaussian',kernel_length,sig);
    disp('Gaussian kernel generated.') ;
elseif strcmp(kernel_type,'Uniform')
    par.psf=fspecial('average',sf);
    disp('Average kernel generated.') ;
else
    error('Wrong kernel type assigned !!!');
end
par.otf=psf2otf(par.psf,szF(1:2));
H=ConvC(data.F,par.otf);
data.H=H(1+shift:sf:end,1+shift:sf:end,:);
if SNRh~=0
    sigmah=sqrt(sum(data.H(:).^2)/(10^(SNRh/10))/numel(data.H));
else
    sigmah=0;
end
data.H=data.H+sigmah*randn(size(data.H));
data.H(data.H<0)=0;
disp('LR-HSI generated.');
data.upH=imresize(data.H,sf);

%%  HR-MSI
par.R=create_R(R);
data.M=tmprod(data.Fm,par.R,3);
if SNRm~=0
    sigmam=sqrt(sum(data.M(:).^2)/(10^(SNRm/10))/numel(data.M));
else
    sigmam=0;
end
rng(s);
data.M=data.M+sigmam*randn(size(data.M));
data.M(data.M<0)=0;
disp('HR-MSI generated.');
par.sf=sf;  par.shift=shift;  

%%  HSRSV
clearvars -except par data flag Tab SNRm SNRh kernel_type sigmah sigmam szF;
fprintf('============================================================= HSRSV ============================================================= \n');
par.K=5;         par.iter=200;     par.lamS=1e-3;   par.lamM=1e0;  
par.lamT=1e-2;   par.lamZ=1e-3;    par.eta=1e-1;    par.lam=1;     
par.gam=10;      par.Scale=0.98;   p={1/2,2/3};     SparseX={'1','mcp'};
flag.Sum2One=0;  flag.SparseX=SparseX{2};           flag.p=p{2};
Time.HSRSV=tic;
    HSRSV=HSRSV_main(data,par,flag);  
Time.HSRSV=toc(Time.HSRSV);

[Eval.PSNR,Eval.RMSE,Eval.ERGAS,Eval.SAM,Eval.UIQI,Eval.SSIM,Eval.DD,Eval.CC] ...
    =quality_assessment(double(im2uint8(data.F)),double(im2uint8(HSRSV)),0,1.0/par.sf);
eval(['Tab.', 'Proposed',  ' =[Eval.PSNR; Eval.RMSE; Eval.ERGAS; Eval.SAM; Eval.UIQI; Eval.SSIM; Eval.DD; Eval.CC; Time.HSRSV];']);
disp(Tab);
clear flag p SparseX

%%  Variable
F=data.F;  H=data.H;   M=data.M;     upH=data.upH;  
R=par.R;   sf=par.sf;  psf=par.psf;  otf=par.otf;  
shift=par.shift;       p=sqrt(diag(psf));
P1=conv2mat(p,szF(1:2));   P1=P1(1+shift:sf:end,:);   
P2=conv2mat(p',szF(1:2));  P2=P2(1+shift:sf:end,:);  
clear data par

%%  BGS-GTF 
fprintf('============================================================= BGS-GTF ============================================================= \n');
BGS_GTF.Pr=cell(1,3);    BGS_GTF.Pr{1,1}=P1;          BGS_GTF.Pr{1,2}=P2;    
BGS_GTF.Pr{1,3}=R;       BGS_GTF.nu=1.2;              BGS_GTF.lambda=1e0;    
BGS_GTF.addatoms=[0,0];  BGS_GTF.turank=[100,100,30]; BGS_GTF.shape_re=[4,25,4,25,2,15];
BGS_GTF.MaxLoop=500;     BGS_GTF.IniLoop=500;         BGS_GTF.ConvTol=1e-3;
BGS_GTF.rho=1e-7;        BGS_GTF.mu=1e0;             
Time.BGS=tic;
    BGS_GTF.HSR=ReTF_mod(H,M,BGS_GTF.Pr,BGS_GTF.lambda,BGS_GTF.mu,BGS_GTF.rho,BGS_GTF.nu,BGS_GTF.turank,BGS_GTF.addatoms,BGS_GTF.shape_re,BGS_GTF.IniLoop,BGS_GTF.MaxLoop,BGS_GTF.ConvTol);
Time.BGS=toc(Time.BGS);

[Eval.PSNR,Eval.RMSE,Eval.ERGAS,Eval.SAM,Eval.UIQI,Eval.SSIM,Eval.DD,Eval.CC] ...
    =quality_assessment(double(im2uint8(F)),double(im2uint8(BGS_GTF.HSR)),0,1.0/sf);
eval(['Tab.', 'BGS_GTF',  ' =[Eval.PSNR; Eval.RMSE; Eval.ERGAS; Eval.SAM; Eval.UIQI; Eval.SSIM; Eval.DD; Eval.CC; Time.BGS];']);
disp(Tab);

%%  FuVar 
fprintf('============================================================= FuVar ============================================================= \n');
p=5;  psfR=rot90(psf,2);
N1=size(H,1);  N2=size(H,2);  M1=size(M,1);  
M2=size(M,2);  L=size(H,3);   data_r=(reshape(H,N1*N2,L)')';
M0=vca_FuVar(data_r','Endmembers',p,'verbose','off');
A_FCLSU=FCLSU(data_r',M0)';   A_init=reshape(A_FCLSU',N1,N2,p);
A_init=imresize(A_init,sf);   A_init=reshape(A_init,M1*M2,p)';
FuVar.lambda_m=1;    FuVar.lambda_a=1e-2;   FuVar.lambda_1=1e0;  
FuVar.lambda_2=1e3;  Psi_init = ones(L,p);
Time.FuVar=tic;
    [FuVar.HSRh,FuVar.HSRm,~,~]=FuVar_main(H,M,A_init,M0,Psi_init,R,sf,psfR,FuVar.lambda_m,FuVar.lambda_a,FuVar.lambda_1,FuVar.lambda_2);
Time.FuVar=toc(Time.FuVar);
FuVar.HSRh=reshape(FuVar.HSRh',M1,M2,L);

[Eval.PSNR,Eval.RMSE,Eval.ERGAS,Eval.SAM,Eval.UIQI,Eval.SSIM,Eval.DD,Eval.CC] ...
    =quality_assessment(double(im2uint8(F)),double(im2uint8(FuVar.HSRh)),0,1.0/sf);
eval(['Tab.', 'FuVar',  ' =[Eval.PSNR; Eval.RMSE; Eval.ERGAS; Eval.SAM; Eval.UIQI; Eval.SSIM; Eval.DD; Eval.CC; Time.FuVar];']);
disp(Tab);
clear p psfR N1 N2 M1 M2 L data_r M0 A_FCLSU A_init Psi_init

%%  GSFus 
fprintf('============================================================= GSFus ============================================================= \n');
GSFus.subspace=3;      GSFus.lambda_1=1e-1; 
GSFus.lambda_2=1e-5;   GSFus.mu=1e-1;       
Time.GSFus=tic;
    GSFus.HSR=GSFus_main(H,M,R,otf,sf,F,GSFus.subspace,GSFus.lambda_1,GSFus.lambda_2,GSFus.mu);
Time.GSFus=toc(Time.GSFus);

[Eval.PSNR,Eval.RMSE,Eval.ERGAS,Eval.SAM,Eval.UIQI,Eval.SSIM,Eval.DD,Eval.CC] ...
    =quality_assessment(double(im2uint8(F)),double(im2uint8(GSFus.HSR)),0,1.0/sf);
eval(['Tab.', 'GSFus',  ' =[Eval.PSNR; Eval.RMSE; Eval.ERGAS; Eval.SAM; Eval.UIQI; Eval.SSIM; Eval.DD; Eval.CC; Time.GSFus];']);
disp(Tab);

%%  BTD-Var 
fprintf('============================================================= BTD-Var ============================================================= \n');
BTD_Var.R=4;        BTD_Var.L=20;       BTD_Var.nIter=150;  
BTD_Var.lamda=1e4;  BTD_Var.r=BTD_Var.L*ones(1,BTD_Var.R);
Time.BTD_Var=tic;
    fprintf('Initializing ...\n');
    BTD_Var.C0=VCA(tens2mat(upH,3),'Endmembers',BTD_Var.R,'SNR',0);
    A=BTD_Var.C0\tens2mat(upH,3);  A=mat2tens(A,[szF(1:2),BTD_Var.R],3);
    [BTD_Var.A0,BTD_Var.B0]=deal(zeros(0));
    for i=1:size(A,3)
        [u,s,v]=svd(A(:,:,i),'econ');
        u=u(:,1:BTD_Var.L);  s=s(1:BTD_Var.L,1:BTD_Var.L);  v=v(:,1:BTD_Var.L);
        BTD_Var.A0=[BTD_Var.A0,deal(u*s^(0.5))];  
        BTD_Var.B0=[BTD_Var.B0,deal(v*s^(0.5))]; 
    end
    [A_hat,B_hat,~,C_hat,~,~,~]=BTD_Var_main(F,H,M,P1,P2,R,BTD_Var.R,BTD_Var.B0,BTD_Var.C0,R*BTD_Var.C0,BTD_Var.nIter,BTD_Var.lamda);
    BTD_Var.HSR=ll1gen({A_hat,B_hat,C_hat},BTD_Var.r);
Time.BTD_Var=toc(Time.BTD_Var);

[Eval.PSNR,Eval.RMSE,Eval.ERGAS,Eval.SAM,Eval.UIQI,Eval.SSIM,Eval.DD,Eval.CC] ...
    =quality_assessment(double(im2uint8(F)),double(im2uint8(BTD_Var.HSR)),0,1.0/sf);
eval(['Tab.', 'BTD_Var',  ' =[Eval.PSNR; Eval.RMSE; Eval.ERGAS; Eval.SAM; Eval.UIQI; Eval.SSIM; Eval.DD; Eval.CC; Time.BTD_Var];']);
disp(Tab);
clear A u s v A_hat B_hat C_hat i
