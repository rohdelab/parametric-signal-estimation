
function mCRLB = fn_CRLB_chirp_delay(sig1fil,win,tm,f, t, dt, w, tau, variance, bias)

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
    
    %x_g = (2*disp(b)*t).*interp1(t,sig1fil,(t.^2)*disp(b) - delay(b),'pchip');
    gwin=exp(-(disp(b)*(t.^2)-tm-delay(b)).^2/(2*win^2));
    x_g=(2*disp(b)*t).*gwin.*sin(2*pi*f*(disp(b)*(t.^2)-delay(b)));
    
    %x_g_w = (2*(disp(b)+dw)*t).*interp1(t,sig1fil,(t.^2)*(disp(b)+dw) - delay(b),'pchip');
    gwin=exp(-((disp(b)+dw)*(t.^2)-tm-delay(b)).^2/(2*win^2));
    x_g_w=(2*(disp(b)+dw)*t).*gwin.*sin(2*pi*f*((disp(b)+dw)*(t.^2)-delay(b)));
    
    %x_g_t = (2*disp(b)*t).*interp1(t,sig1fil,(t.^2)*disp(b) - delay(b)-delt,'pchip');
    gwin=exp(-(disp(b)*(t.^2)-tm-delay(b)-delt).^2/(2*win^2));
    x_g_t=(2*disp(b)*t).*gwin.*sin(2*pi*f*(disp(b)*(t.^2)-delay(b)-delt));

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


