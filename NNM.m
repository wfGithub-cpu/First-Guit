clear all;
close all;
clc
%% signal model
K = 21;%信源数 
source = deg2rad(acosd(linspace(0.9,-0.9,K)));
f0 = 100e6;   % 频率
lc = 3e8;   % 波速
lambda = lc/f0;  % 波长
dl = lambda/2;  % 阵元间距
m = 3;
n = 5;
Omega1 = [0:m:m*(n-1)];
Omega2 = [0:n:n*(2*m-1)];
Omega = unique([Omega1,Omega2]);% 互质阵列
M = length(Omega);  % 阵元数
geoa = Omega.'*dl;  % 阵列几何形状
L = 400;
SNR = 10;
sigmas = 10^(SNR/10);
s = sqrt(sigmas/2)*(randn(K,L) + 1i*randn(K,L));
A = exp(1i*2*pi/lambda*geoa*cos(source));
noise = sqrt(1/2)*(randn(M,L) + 1i*randn(M,L)); % 噪声功率为1
Y = A*s+noise;

R = Y *Y'/L;   % 接收数据协方差矩阵
vR = R(:);  %向量化协方差矩阵
ind = reshape(repmat(Omega.',1,M)-repmat(Omega,M,1),M^2,1);
[C,ia,~] = unique(ind); %找出第一次出现的数据，并排序
vx = vR(ia);
SM = n*(2*m-1)+1;
OnOmega = C((length(C)+1)/2:end)+1;

cvx_quiet true
cvx_precision default
cvx_solver sdpt3
cvx_begin sdp 

  variable Rv(SM,SM) hermitian toeplitz,

  minimize norm_nuc(Rv);
  subject to 
  Rv(OnOmega,1) == vx(C>=0) 

cvx_end
%% MUSIC
[V,D] = eig(Rv); 
[~, order] = sort(abs(diag(D)));
Vp = V(:, order(1:end-K));%噪声空间
Pp = Vp*Vp';
step = 0.01;
thetaptall = [0:step:180];
pmusic = zeros(1,length(thetaptall));
angle = zeros(1,length(thetaptall));
 for iang = 1:length(thetaptall)
   angle(iang) = thetaptall(iang);
   phim = deg2rad(angle(iang));
   atheta = exp(1i*2*pi/lambda*dl*(0:SM-1)'.*cos(phim));
   pmusic(iang) = 10*(log10(abs((atheta'*atheta)/(atheta'*Pp*atheta))));
 end
angle_NNM = PeakFinding(pmusic,K,acosd(1),acosd(-1),step);
mse = sqrt(mean(abs(angle_NNM-rad2deg(source)).^2));
%% plot result
figure(1);
set(gcf,'unit','centimeters','position',[12 8 8.8 5.5]);
set(gca,'Position',[.15 .2 .8 .7]);
plot(angle,pmusic/max(pmusic),'b');
grid on;
hold on;
xlabel('θ(deg)','FontName','Arial','FontSize',8);
ylabel('Spatial Spectrum','FontName','Arial','FontSize',8);
set(gca,'FontName','Arial','FontSize',8);
x = rad2deg(source)';
yy = [0:1/(K-1):1];
for ii = 1:K
    xx(ii,:) = x(ii)*ones(1,K);
    plot(xx(ii,:),yy,'r-.');
    hold on;
end
axis([1 180 0 1])