
function mCRLB = fn_CRLB_delay(sig1fil, t, dt, tau, variance, bias)

delay = tau;
dtau = 1;
delt = dtau*dt;

a=length(delay);

if bias ~= 0
    dbdtau=1;
else
    dbdtau=0;
end
%%
maxrun=length(variance);
for run=1:maxrun
    x_w = interp1(t,sig1fil,t-tau,'pchip');%circshift(sig1fil,delay);
    x_w2 = interp1(t,sig1fil,t-tau-delt,'pchip');%circshift(sig1fil,delay+dtau);
    
    I = (1/variance(run))*(sum(((x_w2-x_w)/delt).^2));
    mCRLB(run,:) = (1+dbdtau)^2/I; % CRLB for tau
end



