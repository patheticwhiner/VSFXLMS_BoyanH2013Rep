function [y,e,Wgt] = FXLMS_NANC(X,d,FilterParams)
    L = FilterParams.Length;
    stepSize = FilterParams.StepSizeConst;
    secPath = FilterParams.SecondaryPath;
    secPathEst = FilterParams.SecondaryPathEst;
    % 由于论文提出了特殊要求，此处应该增设一部分内容，防止中断影响
    Nstr = size(X,1);
    Xrev = zeros(length(secPathEst),L);
    Y = zeros(length(secPath),1);
    W = zeros(L,1);
    Wgt = zeros(L,Nstr);
    y = zeros(1,Nstr);
    e = zeros(1,Nstr);
    
    for i = 1:Nstr
        Xrev = [X(i,:); Xrev(1:end-1,:)];
        y(i) = X(i,:) * W; % 权重和原本的x 
        Y = [y(i);Y(1:end-1)];
        e(i) = d(i) - secPath * Y; % 
        filtX = secPathEst * Xrev; % filtered-x
        W = W + (stepSize.*filtX)'*e(i);
        Wgt(:,i) = W;
    end
end