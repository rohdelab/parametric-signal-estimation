function [f,df,xtilde]=CDT(J0,J1,t,sigma,Es)
% The goal is to find the mass preserving transform of a 1D pdf into another one
%                       df(x)*J1(f(x))=J0(x)
% Input: 
%               J0,J1  = PDF functions (Normalized 1D functions of the same
%                        size). Note that J0 is the template signal.
%               t      = Domain of J1
%               sigma  = Standard deviation of the additive noise
%               Es     = Energy of the input signal
% Output:       f      = The mass preserving mapping (CDT)
%               df     = The Jacobian (in 1D is just the gradient) of f. 
%               xtilde = Domain of the CDT
%

if length(J0) ~= length(J1)
    error('Signals must have same length!')
end
if (sum(J0<0) ~= 0 || sum(J1<0) ~= 0) 
    error('Signals must be non-negative!')
end
if (size(J0,2) == 1)
    J0=J0';
end
if (size(J1,2) == 1) 
    J1=J1';
end
if (size(t,2) == 1) 
    t=t';
end

xk = t;

cJ0 = cumsum(J0)/sum(J0);% cJ0 is the CDF of J0 
cJ1 = cumsum(J1)/sum(J1);% cJ1 is the CDF of J1 

if(sigma) % Do the de-noising if noise standard deviation is known -jmn
    cJ1=(cJ1*(Es+sigma^2*(t(end)-t(1)))-sigma^2*(t-t(1)))/Es;
    
    % [CDF restoration] Preserve the non-decreasing property of CDFs
    [~,ind] = max(cJ1); cJ1(ind:end) = max(cJ1);
    [~,ind] = min(cJ1); cJ1(1:ind) = min(cJ1);
    for p=2:length(cJ1)
        if cJ1(p)<=cJ1(p-1); cJ1(p) = cJ1(p-1)+1e-7; end
    end
    
    cJ1 = (cJ1 - min(cJ1))/(max(cJ1)-min(cJ1));  % min-max normalization
end

xtilde = linspace(0,1,length(J0)); % A Grid for CDF with the same size as J0

%%  
%    Vq = interp1(X,V,Xq) interpolates to find Vq, the values of the
%    underlying function V=F(X) at the query points Xq.  This does the job
%    of finding shat=S^{-1}(t).  Basically takes us from uniform shat
%    values at irregular S(t) spacing to irregular shat values at regular t
%    spacing
%
XJ0 = interp1(cJ0,t,xtilde,'pchip');% Find the xs which match the grided
                                     % CDF for cJ0
XJ1 = interp1(cJ1,t,xtilde,'pchip');% Find the xs which match the grided 
                                     % CDF for cJ1
%%
u = (XJ0-XJ1); % Calculate the translation function u(XJ0)=XJ0-f(XJ0)
u = interp1(XJ0,u,t,'pchip');% Regrid u to get u(x)=x-f(x)
du = gradient(u,1);         % Calculate gradient of u

f = xk-u; %get the mapping f=x-u
df = 1-du;%get the gradient of f
