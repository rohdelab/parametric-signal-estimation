
function mCRLB = fn_CRLB_delay_scale(sig1fil, t, dt, w, tau, variance, bias)

delay = tau;
dtau = 1;
delt = dtau*dt;

disp = w;
dw = 0.001;

a=length(disp);

if bias ~= 0
    dbdtau=1;
else
    dbdtau=0;
end

%%
maxrun=length(variance);
for run=1:maxrun

for b = 1:a   
    
    x_g = interp1(t,sig1fil,t*disp(b) - delay(b),'pchip');
    
    x_g_w = interp1(t,sig1fil,t*(disp(b)+dw) - delay(b),'pchip');
    x_g_t = interp1(t,sig1fil,t*disp(b) - delay(b)-delt,'pchip');

    I11 = sum(((x_g_w-x_g)/dw).^2);
    I22 = sum(((x_g_t-x_g)/delt).^2);
    I12 = sum(((x_g_t-x_g)/delt).*((x_g_w-x_g)/dw));
    I21 = I12;
    
    I = (1/variance(run))*[I11, I12; I21, I22];
    Iinv = (1+dbdtau)^2*pinv(I);
    
    CRLB(b,:) = [Iinv(1,1), Iinv(2,2)]; % CRLB for w and tau respectively
        
end
mCRLB(run,:)=mean(CRLB,1);
end



