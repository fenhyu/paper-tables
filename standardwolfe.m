% Copyright (C) 2021 Guanghui Zhou and Yu Cai
% School of Mathematics and Statistics, Huaibei Normal University
% Huaibei, China
% All rights reserved.
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the 'Software'), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.


% function [alpha] = standardwolfe(fun,dfun,xk,dk)
% rho = 0.01; 
% sigma = 0.1;
% alpha = 1;
% a = 0;
% b = Inf; 
% i = 0; 
% imax = 30;
% while i < imax
%     fold = feval(fun,xk);
%     fnew = feval(fun,xk+alpha*dk);
%     gold = feval(dfun,xk);
%     gnew = feval(dfun,xk+alpha*dk);
%     i = i+1;
%     if ~(fnew<=fold+rho*alpha*gold'*dk)
%         b = alpha;
%         alpha = (alpha+a)/2;
%         continue;
%     end
%     if ~(gnew'*dk >= sigma*gold'*dk)
%         a = alpha;
%         alpha = min([2*alpha, (b+alpha)/2]);
%         continue;
%     end
%     if i>=imax
%         break;
%     end
%     break;
% end
% 
% end



function [alpha, nf, ng, lsflag] = standardwolfe(fun, dfun, xk, dk, nf, ng, tStart, timeLimit)
% Standard Wolfe line search with wall-clock time limit.
% lsflag = 0 (normal), -1 (timeout)

rho   = 0.01;
sigma = 0.1;

alpha = 1;
a = 0;
b = Inf;
i = 0;
imax = 30;

lsflag = 0;

if nargin < 8 || isempty(timeLimit)
    timeLimit = 0;   % 0 means no limit
end

% ===== timeout check (before any eval) =====
if timeLimit > 0 && toc(tStart) > timeLimit
    lsflag = -1;
    alpha  = 0;
    return;
end

% ---- initial evaluations ----
fold = feval(fun, xk);
nf = nf + 1;

gold = feval(dfun, xk);
ng = ng + 1;

while i < imax

    % ===== timeout check (inside line search loop) =====
    if timeLimit > 0 && toc(tStart) > timeLimit
        lsflag = -1;
        alpha  = 0;
        return;
    end

    xnew = xk + alpha * dk;

    fnew = feval(fun, xnew);
    nf = nf + 1;

    gnew = feval(dfun, xnew);
    ng = ng + 1;

    i = i + 1;

    % ---- Armijo condition ----
    if ~(fnew <= fold + rho * alpha * (gold' * dk))
        b = alpha;
        alpha = (alpha + a) / 2;
        continue;
    end

    % ---- Wolfe curvature condition ----
    if ~(gnew' * dk >= sigma * (gold' * dk))
        a = alpha;
        alpha = min([2 * alpha, (b + alpha) / 2]);
        continue;
    end

    break;
end

end