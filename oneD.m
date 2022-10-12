%AUTHOR:JIANG BAICHUAN
%2022/10/3
%

%% Normal Epi,M,Endo initialization
clc;clear;close all;
Volt = -86.2; Cai = 0.00007; CaSR = 1.3; CaSS = 0.00007; Nai = 7.67;
Ki = 138.3; INa_m = 0.; INa_h = 0.75; INa_j = 0.75; IKr_xr1 = 0.;
IKr_xr2 = 1.; IKs_xs = 0.; Ito_r = 0.; Ito_s = 1.; ICaL_d = 0.;
ICaL_f = 1.; ICaL_f2 = 1.; ICaL_fCaSS = 1.; RR = 1.;
CellNumber = 50;
X0_ini = [Volt Cai CaSR CaSS Nai Ki INa_m INa_h INa_j IKr_xr1 IKr_xr2 IKs_xs ...
    Ito_r Ito_s ICaL_d ICaL_f ICaL_f2 ICaL_fCaSS RR];%第一个细胞的初始状态

%% VenCon and VenDif set up
% D,G set up, D:cm^2/ms, Cm:pF
deltax = 0.025;
Dmm = 0.00154; 
Gmm = Dmm./(deltax); 
%for w=1:1:300
%% 外部电刺激变化
    StimNum=4;
    timegap=350;
    %timegap=250+w;
    %FileName=num2str(timegap);
    IexMatrix=zeros(2,StimNum);
for i=1:StimNum
    IexMatrix(1,i)= 100 + (i-1)*timegap;
    IexMatrix(2,i) = IexMatrix(1,i)+1;
end

%% 6 M(EAD)-F coupling no legend
dur = 1800;
celltype = 'M';
NumF = 4;
Ggap = 3; Gf = 0.2; Ef = -50;
Cf = 6.3; Cm = 185; Gmf = Ggap./Cm; Gfm = Ggap./Cf;

%是否EAD化
PGkr = 1; PGks = 1; PGCaL = 1; PGNa = 1; 
%PGkr = 0.5; PGks = 0.5; PGCaL = 2.7; PGNa = 1;
dT = 0.02; dT_rev = 1/dT;

%储存单元
X0=zeros(CellNumber,19);
col=CellNumber;
parfor q = 1:col
    X0(q,:)=X0_ini;
end
APRS = zeros(dur,CellNumber+1); 

%Istim_mm=zeros(1,CellNumber);
Iex=zeros(CellNumber,1);
count=1;
t = 0; 
VF = Ef*ones(NumF,1);
for step = 1:dur/dT
    t = t + dT; Istim_fm = 0;
    %count=count+1;
    % Update FtM
    if NumF > 0
        for i = 1:NumF
            Istim_fm = Istim_fm + (VF(i) - X0(1))*Gmf;
        end
    end
    % Update fibrosis
    for i = 1:NumF
        If = Gf*(VF(i) - Ef);
        Istim_mf = (X0(1) - VF(i))*Gfm;
        dv = (Istim_mf - If)/Cf;
        VF(i) = dv*dT + VF(i);
    end
    % Update MtM
    parfor j=2:CellNumber
        Iex(j)=(X0(j)-X0(j-1))*Gmm;
    end
    
    % update Istim_ex
    amp = -52.0; Istim_ex = 0;
    for k = 1:size(IexMatrix,2)
        if t > IexMatrix(1,k)- 1e-6 && t < IexMatrix(2,k) + 1e-6
            Istim_ex = amp;
        end
    end
    %外部刺激相加
    Iex(1) = Istim_ex-Istim_fm;
   
    % Update myoctle Input:t,X0,flag,PGkr,PGks,PGCaL
    
    X0(1,:)= Tnnp06_model(t,X0(1,:),1,dT,PGkr,PGks,PGCaL,PGNa,celltype,Iex(1));
    for z=2:CellNumber
    X0(z,:) = Tnnp06_model(t,X0(z,:),1,dT,PGkr,PGks,PGCaL,PGNa,celltype,Iex(z));
    end
    if rem(step,dT_rev) == 0
        APRS(step*dT,1) = t; 
        for c=1:CellNumber
        APRS(step*dT,c+1) = X0(c);
        end
    end
end

%figure(w)
plot(APRS(:,1),APRS(:,2),APRS(:,1),APRS(:,50),'linewidth',2);
title('drug-induced heterogeneity')
%set(gca,'looseInset',[0 0 0 0]);
%set(gcf,'position',[50 100 500 500],'color',[1 1 1]);
axis([-inf inf -90 50]); 
%legend('No fib','1 fib','2 fib','3 fib','4 fib')
%saveas(w,FileName,'jpg')
%end
