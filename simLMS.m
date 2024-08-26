% simLMS.m --- Designed by linmh0130@stu.xjtu.edu.cn
%
% Description
%   function [y,e,W] = simLMS(x,d,mu,W0)
%   function [y,e,W] = simLMS(x,d,mu,W0,v)
%       A simple LMS Adaptive Filter
% Update equation
%       e(i) = d(i) - y(i);
% or    e(i) = d(i) + v(i) - y(i);
%
%       W(next) = W(now) + (mu*e(i)).*u;
% Parameters
%   x:          Input signal
%   d:          Desired output
%   mu:         Stepsize
%   W0:         Initial value of weights
% Return
%   y:          Output signal
%   e:          Error
%   W:          Weights

function [y,e,wgt,varargout] = simLMS(x,d,mu,W0)
    W = W0;
    L = length(W);
    y = zeros(1,length(x));
    e = zeros(1,length(x));
    u = zeros(1,L); % input_sav vector
    % 小试牛刀
    wgt = zeros(L,length(x));         % 初始化权重记录
    R = zeros(L, L);
    P = zeros(L, 1);
    Rd = 0;
    for i = 1 : length(x)
        u = [x(i) u(1:end-1)];
        y(i) = W * u';
        e(i) = d(i)- y(i);
        W = W + (mu*e(i)).*u;
        
        if i >= L
            u_vec = u(:); % 当前输入向量
            R = R + u_vec * u_vec'; % 输入自相关矩阵估计
            P = P + d(i) * u_vec;   % 输入、目标互相关向量
            Rd = Rd + d(i)*d(i);
        end
        wgt(:, i) = W;
    end
    R = R / (length(x) - L + 1);
    P = P / (length(x) - L + 1);
    Rd = Rd / (length(x) - L + 1) - P.'*inv(R)*P;
    % 根据 nargout 确定是否输出 R 和 P
    if nargout == 5
        varargout{1} = R;
        varargout{2} = P;
    end
end
