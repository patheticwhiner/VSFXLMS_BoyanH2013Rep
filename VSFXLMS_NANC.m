function [y,e,varStep,Wgt] = VSFXLMS_NANC(X,d,FilterParams)
    L = FilterParams.Length;
    step = FilterParams.VariStepParams(1,:);
    epl = FilterParams.VariStepParams(2,:);
    eta = FilterParams.VariStepParams(3,:);
    H = FilterParams.SecondaryPath;
    Hest = FilterParams.SecondaryPathEst;
    % 界限还是必须要设的
    mu_min = 5e-2;
    mu_max = 0.5;
    % 由于论文提出了特殊要求，此处应该增设一部分内容，防止中断影响
    nStr = size(X,1);
    Y = zeros(length(H),1);
    Xrev = zeros(length(Hest),L);
    W = zeros(L,1);
    varStep = zeros(L,nStr);
    Wgt = zeros(L,nStr);
    for i = 1:nStr
        if isnan(W)
            W = W;
        end
        Xrev = [X(i,:); Xrev(1:end-1,:)];
        y(i) = X(i,:)*W; % 权重和原本的x 
        Y = [y(i);Y(1:end-1)];
        e(i) = d(i) - H*Y; % 
        x = Hest*Xrev; % filtered-x
        if i==1
            step = epl.*step + eta.*(e(i).^2*x(1).^2);
        else
            step = epl.*step + eta.*(e(i).*e(i-1).*x(1).*x(2));
        end
        step(step<mu_min) = mu_min;
        step(step>mu_max) = mu_max;
        varStep(:,i) = step;
        W = W + (step.*x)'*e(i);
        if abs(e(i))>1
            W = W;
        end
        Wgt(:,i) = W;
    end
end