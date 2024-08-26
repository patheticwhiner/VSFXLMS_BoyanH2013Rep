clear; close all; clc;
%% 次级通路建模
wc = 0.4*pi;
M = 64;
window = blackman(M + 1).';
S = fir1(M,wc/pi,window);
%% 次级通路辨识
Nitr = 25000; 
stepLMS = 0.0025;
randomSignal = mvnrnd(0,1,Nitr) + 0.01*mvnrnd(0,0.1,Nitr);
bandpassSignal = filter(S,1,randomSignal);
W0 = zeros(1,68+1);
[y,e,wgt] = simLMS(randomSignal,bandpassSignal,stepLMS,W0);
Sest = mean(wgt(:,end-999:end),2);
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
d1 = X1*Coeffs + mvnrnd(0,0.1,Nstr/2);
d2 = X2*(-Coeffs) + mvnrnd(0,0.1,Nstr/2);
%% 最优权重计算
a_opt = zeros(1,freqN);
b_opt = zeros(1,freqN);
for k = 1:freqN
    nn = 0:M;
    xa = cos(nn'*freqW(k));
    xb = sin(nn'*freqW(k));
    a = S*xa;
    b = S*xb;
    Y = inv([a,-b;b,a])*[coeffsA(k);coeffsB(k)];
    a_opt(k) = Y(1); b_opt(k) = Y(2);
end
Wopt = [a_opt,b_opt]';

%% FXLMS/FXRLS/VSS FXLMS滤波器设计
FilterParams.Length = 2*freqN;
FilterParams.SecondaryPath = S;
FilterParams.SecondaryPathEst = Sest.';
FilterParams.StepSizeConst = 0.0035*ones(1,2*freqN);
FilterParams.ForgetConst = 0.999;
varU = 0.002*ones(1,2*freqN);
epl = 0.999*ones(1,2*freqN);
eta = 0.0001*ones(1,2*freqN);
FilterParams.VariStepParams = [varU;epl;eta];
%% 滤波
for filtChoice = 1:3
    switch filtChoice
        case 1            
            % FXLMS 滤波
            [y1,e1,Wgt1] = FXLMS_NANC(X1,d1,FilterParams);
            [y2,e2,Wgt2] = FXLMS_NANC(X2,d2,FilterParams);
            text = 'FXLMS';
        case 2
            % FXRLS 滤波
            [y1,e1,Wgt1] = FXRLS_NANC(X1,d1,FilterParams);
            [y2,e2,Wgt2] = FXRLS_NANC(X2,d2,FilterParams);
            text = 'FXRLS';
        case 3 
            % VSS FXLMS 滤波
            [y1,e1,varStep1,Wgt1] = VSFXLMS_SIMP1(X1,d1,FilterParams);
            [y2,e2,varStep2,Wgt2] = VSFXLMS_SIMP1(X2,d2,FilterParams);
            varStep = [varStep1, varStep2];
            text = 'VSS FXLMS';
    end
    % 滤波图示
%     figure(filtChoice);
%     subplot(2,1,1); plot(n,[d1;d2]); hold on; plot(n,[y1,y2]); plot(n,[e1,e2]);
%     xlabel('Iteration number n'); legend('Primary noise','Output','Error'); title(text);
%     subplot(2,1,2); plot(Wopt,'-o'); hold on; plot(Wgt1(:,end),'--s'); plot(-Wgt2(:,end),':*');
%     xlabel('DFC sequence number'); ylabel('Weight'); legend('Optimal Weight','1st Half','2nd Half'); title(text);
%     tightfig;
    % 绘制Fig
%     Figdemo1(n,Wgt1,Wgt2,a_opt,b_opt,freqN,filtChoice);
    Figdemo2(n,Wgt1,Wgt2,a_opt,filtChoice);
%     figure(3+3+filtChoice); Figdemo3(n,e1,e2); title(text);saveas(gcf,[text,'_errMSE.svg']);
end
% figure(3+3+3+1); Figdemo4(n,varStep,freqN);
% [y3,e3,varStep3,Wgt3] = VSFXLMS_SIMP1(X1,d1,FilterParams,'SIMP1');
% [y4,e4,varStep4,Wgt4] = VSFXLMS_SIMP1(X2,d2,FilterParams,'SIMP1');
% [y5,e5,varStep5,Wgt5] = VSFXLMS_SIMP1(X1,d1,FilterParams,'SIMP2');
% [y6,e6,varStep6,Wgt6] = VSFXLMS_SIMP1(X2,d2,FilterParams,'SIMP2');
% Figdemo5(n,Wgt1,Wgt2,Wgt3,Wgt4,Wgt5,Wgt6,a_opt,b_opt);
%% 
function Figdemo1(n,Wgt1,Wgt2,a_opt,b_opt,freqN,filtChoice)
    figure(3+filtChoice); 
    subplot(1,2,1); hold on;
    for k = 1:freqN
        plot(n,20*log10([(Wgt1(k,:)-a_opt(k)).^2,(Wgt2(k,:)+a_opt(k)).^2])); 
    end
    legend('E[\epsilon^2_{a_1}(n)]','E[\epsilon^2_{a_2}(n)]','E[\epsilon^2_{a_3}(n)]');
    subplot(1,2,2); hold on;
    for k = 1:freqN
        plot(n,20*log10([(Wgt1(k+freqN,:)-b_opt(k)).^2,(Wgt2(k+freqN,:)+b_opt(k)).^2])); 
    end
    legend('E[\epsilon^2_{b_1}(n)]','E[\epsilon^2_{b_2}(n)]','E[\epsilon^2_{b_3}(n)]');
end
%%
function Figdemo2(n,Wgt1,Wgt2,a_opt,filtChoice)
    for k = 1:3
        figure(3+3+k); hold on;
        plot(n,20*log10([(Wgt1(k,:)-a_opt(k)).^2,(Wgt2(k,:)+a_opt(k)).^2])); 
        if filtChoice == 3
            legend('FXLMS','FXRLS','VSS-FXLMS');
            xlabel('Iteration number n'); ylabel('Estimation MSE of DFCs [dB]');
            title(['E[\epsilon^2_',num2str(k),'(n)]']);
            grid on; tightfig; saveas(gcf,['DFC_MSE',num2str(k),'.svg']);
        end
    end
end
%%
function Figdemo3(n,e1,e2)
    plot(n,20*log10([e1.^2,e2.^2]));
    xlabel('Iteration number n');
    ylabel('Mean squared residual noise [dB]');
end
%%
function Figdemo4(n,varStep,freqN)
    subplot(2,1,1); plot(n,varStep(1,:)); hold on; plot(n,varStep(1+freqN,:),'--');
    xlabel('Iteration number n'); ylabel('Variable step size'); legend('\mu_{a_1}(n)','\mu_{a_2}(n)');
    subplot(2,1,2); plot(n,varStep(1,:)); hold on; plot(n,varStep(2,:),'--'); plot(n,varStep(3,:),':');
    legend('\mu_{a_1}(n)','\mu_{a_2}(n)','\mu_{a_3}(n)'); xlabel('Iteration number n'); ylabel('Variable step size');
end
%%
function Figdemo5(n,Wgt1,Wgt2,Wgt3,Wgt4,Wgt5,Wgt6,a_opt,b_opt)
    k = 1;
    figure(3+3+3+1+1); hold on;
    plot(n,20*log10([(Wgt1(k,:)-a_opt(k)).^2,(Wgt2(k,:)+a_opt(k)).^2])); 
    plot(n,20*log10([(Wgt3(k,:)-a_opt(k)).^2,(Wgt4(k,:)+a_opt(k)).^2])); 
    plot(n,20*log10([(Wgt5(k,:)-a_opt(k)).^2,(Wgt6(k,:)+a_opt(k)).^2])); 
    xlabel('Iteration number n'); ylabel('Estimation MSE of DFCs [dB]');
    legend('VSS-FXLMS','VSS-FXLMS-1','VSS-FXLMS-2'); title('E[\epsilon_{a_1}^2(n)]');
    figure(3+3+3+1+2); hold on;
    plot(n,20*log10([(Wgt1(k+3,:)-b_opt(k)).^2,(Wgt2(k+3,:)+b_opt(k)).^2])); 
    plot(n,20*log10([(Wgt3(k+3,:)-b_opt(k)).^2,(Wgt4(k+3,:)+b_opt(k)).^2])); 
    plot(n,20*log10([(Wgt5(k+3,:)-b_opt(k)).^2,(Wgt6(k+3,:)+b_opt(k)).^2])); 
    xlabel('Iteration number n'); ylabel('Estimation MSE of DFCs [dB]');
    legend('VSS-FXLMS','VSS-FXLMS-1','VSS-FXLMS-2'); title('E[\epsilon_{b_1}^2(n)]');
end