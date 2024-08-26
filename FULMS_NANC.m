function [y,e,Wgt] = FULMS_NANC(X,d,FilterParams)
    L = FilterParams.Length;
    stepSize = FilterParams.StepSizeConst;
    secPath = FilterParams.SecondaryPath;
    secPathEst = FilterParams.SecondaryPathEst;
    % 由于论文提出了特殊要求，此处应该增设一部分内容，防止中断影响
    secPathN = secPath(1,:);
    secPathD = secPath(2,:);
    secPathEstN = secPathEst(1,:);
    secPathEstD = secPathEst(2,:);
    % 
    Nstr = size(X,1);
    filtX = zeros(1,L);
    Xrev = zeros(size(secPathEst,2),L);
    filtXrev = zeros(size(secPathEst,2),L);
    y = zeros(1,Nstr);
    filty = 0;
    yrev = zeros(size(secPath,2),1);
    filtyrev = zeros(size(secPath,2),1);
    e = zeros(1,Nstr);
    W = zeros(L,1);
    Wgt = zeros(L,Nstr);
    % 
    for i = 1:Nstr
        Xrev = [X(i,:); Xrev(1:end-1,:)];
        filtXrev = [filtX; filtXrev(1:end-1,:)];
        filtX = secPathEstD * Xrev + secPathEstN * filtXrev; % filtered-x
        y(i) = X(i,:) * W;
        yrev = [y(i);yrev(1:end-1)];
        filtyrev = [filty; filtyrev(1:end-1,:)];
        filty = secPathD * yrev + secPathN * filtyrev;
        e(i) = d(i) - filty; % 
        W = W + (stepSize.*filtX)'*e(i);
        Wgt(:,i) = W;
    end
end