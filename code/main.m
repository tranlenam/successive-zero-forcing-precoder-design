rng(1)
clear 
clc
N=8;
Nr =[2 2 2 2];
K=length(Nr);
P=10; % dBW
P=10.^(P./10);
% randomly genreate channels
%{ 
H=cell(K,1);
for k=1:K
    H{k}=1/sqrt(2)*(randn(Nr(k),N)+1i*randn(Nr(k),N));
end
%}
load('../data/channel.mat');
%cvx_solver mosek
Nr_cum_sum=cumsum([0,Nr]);
Nr_bar=N-Nr_cum_sum;

MAX_ITER=15;
InitialPrecoder_AllZero=cell(K,1);
for k=1:K
    InitialPrecoder_AllZero{k}=0;
end
SumRate_Zero=Algorithm1(H,P,MAX_ITER,InitialPrecoder_AllZero);
% Random initialization
InitialPrecoder_Random1=cell(K,1);
for k=1:K
    T=rand(Nr_bar(k),Nr(k));
    InitialPrecoder_Random1{k}=T*T';
end
SumRate_Random1=Algorithm1(H,P,MAX_ITER,InitialPrecoder_Random1);

% Another random initialization
InitialPrecoder_Random2 = cell(K,1);
for k=1:K
    T=rand(Nr_bar(k),Nr(k));
    InitialPrecoder_Random2{k}=T*T';
end
SumRate_Random2=Algorithm1(H,P,MAX_ITER,InitialPrecoder_Random2);

iter=1:length(SumRate_Zero);
linewidth=0.8;
h1=figure('Units','inches','Pos',[2 2  6 3]);

plot(iter,SumRate_Zero,'-ro','LineWidth',linewidth,'MarkerEdgeColor','r',...
                'MarkerFaceColor','r','MarkerSize',4.5);
hold on
plot(iter,SumRate_Random1,'-bs','LineWidth',linewidth,'MarkerEdgeColor','b',...
                'MarkerFaceColor','k','MarkerSize',4.5);
plot(iter,SumRate_Random2,'-md','LineWidth',linewidth,'MarkerEdgeColor','m',...
                'MarkerFaceColor','m','MarkerSize',4.5);

xlim([1-0.1 MAX_ITER]);
ymin=min([SumRate_Zero(1) SumRate_Random1(1) SumRate_Random2(1)]);
ymax=max([SumRate_Zero(end) SumRate_Random1(end) SumRate_Random2(end)]);
ylim([ymin-0.2 ymax+0.1]);
xlabel('Iteration','fontsize',8);
ylabel('Sum rate (b/z/Hz)','fontsize',8);
h=legend('Algorithm 1, zero-matrix initialization','Algorithm 1, random initialization',...
    'Algorithm 1, random initialization','Location','southeast');
legend('boxoff')
set(h,'fontsize',8)
set(gca,'fontsize',7);
saveas(gcf,'../results/convergence.png')
