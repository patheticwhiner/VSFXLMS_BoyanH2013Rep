clear; close all; clc;
%% 次级通路建模，免去辨识
Wc = 0.4*pi;
M = 24;
[numz,denz] = butter(M,Wc/pi);
S = [numz;denz];
Sest = S; 
%% 仿真
freqN = 3; 
freqW = 0.1*pi*[1,2,3];
Nstr = 16000;
n = 1:Nstr;
xa = cos(n'*freqW);
xb = sin(n'*freqW);
X = [xa,xb];
coeffsA = [2 1 0.5]';
coeffsB = [-1 -0.5 0.1]';
Coeffs = [coeffsA; coeffsB];
X1 = X(1:Nstr/2,:); X2 = X(Nstr/2+1:end,:);
d1 = X1*Coeffs + mvnrnd(0,0.2,Nstr/2);
d2 = X2*(-Coeffs) + mvnrnd(0,0.2,Nstr/2);
%% 最优权重计算
% IIR滤波器的最优权重暂时还未经推导
%% FULMS滤波器设计
FilterParams.Length = 2*freqN;
FilterParams.StepSizeConst = 0.00025*[1 1 2 2 1 1];
FilterParams.SecondaryPath = S;
FilterParams.SecondaryPathEst = Sest;
%% 滤波
[y1,e1,Wgt1] = FULMS_NANC(X1,d1,FilterParams);
[y2,e2,Wgt2] = FULMS_NANC(X2,d2,FilterParams);
%% 滤波图示
figure;
subplot(2,1,1); 
plot(n,[d1;d2]); hold on; 
plot(n,[e1,e2]);
xlabel('Iteration number n'); legend('Primary noise','Output','Error');
subplot(2,1,2);
plot(n,[y1,y2]); hold on; 
plot(n,[e1,e2]);
% subplot(2,1,2); plot(Wopt,'-o'); hold on; plot(Wgt1(:,end),'--s'); plot(-Wgt2(:,end),':*');
% xlabel('DFC sequence number'); ylabel('Weight'); legend('Optimal Weight','1st Half','2nd Half');
tightfig;