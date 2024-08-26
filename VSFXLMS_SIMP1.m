function [y,e,varStep,Wgt] = VSFXLMS_SIMP1(X,d,FilterParams,varargin)
    L = FilterParams.Length;
    step = FilterParams.VariStepParams(1,:);
    epl = FilterParams.VariStepParams(2,:);
    eta = FilterParams.VariStepParams(3,:);
    H = FilterParams.SecondaryPath;
    Hest = FilterParams.SecondaryPathEst;
    % 检查step,epl,eta长度是否相同
    if nargin > 3
        switch varargin{1}
            case 'SIMP1'
                mode = 'SIMP1';
                stepL = L/2;
            case 'SIMP2'
                mode = 'SIMP2';
                stepL = 1;
            otherwise
                error('输入参数不合要求');
        end
    else
        mode = 'ORIG';
        stepL = L;
    end
    step = step(1:stepL);
    epl = epl(1:stepL);
    eta = eta(1:stepL);
    % 界限还是必须要设的
    mu_min = 5e-4;
    mu_max = 0.5;
    % 由于论文提出了特殊要求，此处应该增设一部分内容，防止中断影响
    nStr = size(X,1);
    Y = zeros(length(H),1);
    Xrev = zeros(length(Hest),L);
    W = zeros(L,1);
    varStep = zeros(stepL,nStr);
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
        if i == 1 
            xpri = x;
            epri = e(i);
        end
        switch mode
            case 'ORIG'
                step = epl.*step + eta*e(i)*epri.*(x.*xpri);
                stp = sat(step,mu_min,mu_max);
            case 'SIMP1'
                step = epl.*step + eta*e(i)*epri.*(x(1:L/2).*xpri(1:L/2));
                stp = sat([step,step],mu_min,mu_max);
            case 'SIMP2'
                step = epl.*step + eta*e(i)*epri.*(x(1).*xpri(1));
                stp = repmat(sat(step,mu_min,mu_max),[1,L]);
        end
        xpri = x;
        epri = e(i);
        varStep(:,i) = step;
        W = W + (stp.*x)'*e(i);
        if abs(e(i))>1
            W = W;
        end
        Wgt(:,i) = W;
    end
end