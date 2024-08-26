function [y,e,Wgt] = FXRLS_NANC(X,d,FilterParams)
    L = FilterParams.Length;
    lambda = FilterParams.ForgetConst;
    H = FilterParams.SecondaryPath;
    Hest = FilterParams.SecondaryPathEst;
    % 由于论文提出了特殊要求，此处应该增设一部分内容，防止中断影响
    nStr = size(X,1);
    Y = zeros(length(H),1);
    Xrev = zeros(length(Hest),L);
    F = eye(L);
    W = zeros(L,1);
    Wgt = zeros(L,nStr);
    
    for i = 1:nStr
        y(i) = X(i,:)*W; % 权重和原本的x 
        Y = [y(i);Y(1:end-1)];
        e(i) = d(i) - H*Y; % 误差
        Xrev = [X(i,:); Xrev(1:end-1,:)];
        filtX = Xrev'*Hest'; % filtered-x
        for r = 1:2:L-1
            F(r:r+1,r:r+1) = [F(r:r+1,r:r+1)-...
                (F(r:r+1,r:r+1)*filtX(r:r+1)*filtX(r:r+1)'*F(r:r+1,r:r+1))/...
                (lambda+filtX(r:r+1)'*F(r:r+1,r:r+1)*filtX(r:r+1))]./lambda;
            W(r:r+1) = W(r:r+1) + (F(r:r+1,r:r+1)*filtX(r:r+1))*e(i)/...
                (lambda+filtX(r:r+1)'*F(r:r+1,r:r+1)*filtX(r:r+1));
        end
        Wgt(:,i) = W;
    end
end