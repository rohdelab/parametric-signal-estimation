%%
% Generate plots (simulation results) from the paper
%
% Reference: Parametric Signal Estimation Using the Cumulative Distribution Transform
%
%%
clc; clear; close all;
addpath('../')


% Select any of the estimation problems listed below. Uncomment the
% corresponding line and then comment the other three.

param='delay';                 %gp(t) = t - tau;    Fig. 6
%param='delay_dispersion';      %gp(t) = wt - tau;   Fig. 7
%param='delay_quadtratic';      %gp(t) = kt^2 - tau; Fig. 10


% Number of realizations for estimating. Results shown in the paper
% were generated using 1000 realizations
NumRealizations=10;


N=400;                          % number of points in signal
dt=0.025;                       % timestep
Fs = 1/dt;
eps=1e-12;                      % "small" value for use in CDT estimation
clip_CDF=25;                    % Number of points to clip off beginning and end of CDF estimate (where no signal present)

if strcmp(param, 'delay_quadtratic')
    t=0:dt:(N-1)*dt;            % time array for chirp (can't be 0 centered)
    tm=t(end)/2;
else
    t=-N/2*dt:dt:(N/2-1)*dt;    % time array (0 centered)
    tm=0;                       % initial centering of pulse
end


%% Define original signal (before time delay)
f=1;                                % modulation frequency
win=1;                              % Width of the pulse (all models assume this =1, don't change)

gwin=exp(-(t-tm).^2/(2*win^2));     % Apodization function
z=gwin.*sin(2*pi*f*t);              % The clean input signal
s=z.^2/sum(z.^2);                   % Squared, normalized input signal.

if strcmp(param, 'delay_quadtratic')
    gwin=exp(-((t.^2)-tm).^2/(2*win^2));
    sc=gwin.*sin(2*pi*f*(t.^2));
    Es = mean((2*t.*sc).^2)*(t(end)-t(1));
else
    Es = mean(z.^2)*(t(end)-t(1));
end


s0=ones(1,N);                       % Reference signal
[shat,df]=CDT(s0,s+eps,t,0,Es);     % CDT of clean "PDF" signal s


%% Define the true parameter values
if strcmp(param, 'delay_quadtratic')
    tau=50.3*dt;                 % true delay (make this non-integer # of time-steps to make discrete XC a biased estimator)
    omega=0.5;
else
    tau=10.3*dt;
    omega=0.75;
end

switch param
    case 'delay'
        % Define signal after time delay
        gwin=exp(-(t-tm-tau).^2/(2*win^2));                     % Apodization function
        zg=gwin.*sin(2*pi*f*(t-tau));                           % altered input signal
        Esg=sum(zg.^2);                                         % Signal energy - should be the same for both input and propagated signal
    case 'delay_dispersion'
        % Define signal after time delay+linear dispersion
        gwin=exp(-(omega*t-tm-tau).^2/(2*win^2));               % Apodization function
        zg=gwin.*sin(2*pi*f*(omega*t-tau));                     % altered input signal
        Esg=sum(zg.^2);
    case 'delay_quadtratic'
        % Define time delay+quadratic dispersion
        gwin=exp(-(omega*(t.^2)-tm-tau).^2/(2*win^2));          % Apodization function
        zg=(2*omega*t).*gwin.*sin(2*pi*f*(omega*(t.^2)-tau));   % altered input signal
        Esg=sum(zg.^2);
end


%% Define SNR range to look at
SNRdb = 0:3:15;
SNR = 10.^(SNRdb/10);
sigmavals =(pi)^(1/4)*sqrt(win)./(sqrt(2*SNR*N*dt))    % standard deviation


%% CRLB
switch param
    case 'delay'
        cCRLB = fn_CRLB_delay(z, t, dt, tau, sigmavals.^2, 0);
        CRLB_tau = cCRLB;
    case 'delay_dispersion'
        cCRLB = fn_CRLB_delay_scale(z, t, dt, omega, tau, sigmavals.^2, 0);
        CRLB_tau = cCRLB(:,2);
        CRLB_omega = cCRLB(:,1);
    case 'delay_quadtratic'
        cCRLB = fn_CRLB_chirp_delay(z,win,tm,f, t, dt, omega, tau, sigmavals.^2, 0);
        CRLB_tau = cCRLB(:,2);
        CRLB_omega = cCRLB(:,1);
end


%% Loop over different SNR levels
indx=1;
for sigma=sigmavals
    disp(['SNR is ',num2str(Es/sigma^2/(t(end)-t(1))),' which should match prescribed SNR of ',num2str(SNR(indx))]); % This should match the SNR
    
    tau_est = zeros([1,NumRealizations]);
    omega_est = zeros([1,NumRealizations]);
    for ind=1:NumRealizations
        sigma
        rng shuffle
        eta=sigma*randn(1,N);       % The noise
        zgn=zg+eta;                 % The altered clean signal + noise
        
        
        %% Estimation maximizing wide-band ambiguity function (WBAF)
        switch param
            case 'delay'
                x0 = 0;
                % local minimization
                fn_zg = @(x) exp(-(t-x-tm).^2/(2*win^2)).*sin(2*pi*f*(t-x));
                fn = @(x) -abs(sum(zgn.*fn_zg(x))).^2;
                lb = 0;
                ub = t(end);
                
                opts = optimoptions('fmincon','Algorithm','sqp');   % sqp = sequential quadratic programming
                xr = fmincon(fn,x0,[],[],[],[],lb,ub,[],opts);
                tau_AF_local(ind)=xr;
                
                % global search
                xr = maximize_AF(zgn, t, tm, f, win, x0);
                tau_AF(ind)=xr;
            case 'delay_dispersion'
                x0 = [1, 0];
                % local minimization
                fn_zg = @(x) exp(-(x(1)*t-x(2)-tm).^2/(2*win^2)).*sin(2*pi*f*(x(1)*t-x(2)));
                fn = @(x) -abs(sqrt(abs(x(1)))*sum(zgn.*fn_zg(x))).^2;
                lb = [0.001, 0];
                ub = [20, t(end)];
                
                opts = optimoptions('fmincon','Algorithm','sqp');
                xr = fmincon(fn,x0,[],[],[],[],lb,ub,[],opts);
                omega_AF_local(ind)=xr(1);
                tau_AF_local(ind)=xr(2);
                
                % global search
                xr = maximize_AF(zgn, t, tm, f, win, x0);
                omega_AF(ind)=xr(1);
                tau_AF(ind)=xr(2);
        end
        
        
        %% XC estimator
        [XC,lags2]=xcorr(zgn,z,N);          % cross correlation over all possible "2N" lags
        [mvp,iposXC]=max(XC);
        tau_XC1(ind)=dt*(iposXC-(N+1));     % could take this as final estimator, however folks claim better accuracy using Hilbert envelope
        
        
        % XC + envelope detection estimator
        xh=hilbert(XC);
        [mvp,iposXC2]=max(abs(xh));         % envelope detector
        tau_XC2(ind)=dt*(iposXC2-(N+1));
        
        
        %% CDT estimation process
        r=zgn.^2/sum(zgn.^2);                       % Squared, altered signal + noise
        [rhat,df2]=CDT(s0,r+eps,t,sigma,Es);        % Noise-corrected CDT of r=sg+noise
        
        midrange=clip_CDF:length(rhat)-clip_CDF;    % Restrict domain of CDT to use in estimation.  Has pronounced effect on bias but not variance
        
        switch param
            case 'delay'
                tau_est(ind)=mean(rhat(midrange))-mean(shat(midrange)); % Estimator is just mean difference of CDTs
                tau_est_f(ind)=tau_est(ind);
            case 'delay_dispersion'               
                B = [rhat(midrange)', -1*ones(length(midrange),1)];
                estimate = (B'*B)\(B'*shat(midrange)');
                omega_est(ind) = estimate(1);
                tau_est(ind)=estimate(2);
            case 'delay_quadtratic'                
                B = [(rhat(midrange).^2)', -1*ones(length(midrange),1)];
                c = (B'*B)\(B'*shat(midrange)');
                omega_est(ind) = c(1);
                tau_est(ind) = c(2);
        end
        
        %% Subspace methods
        if strcmp(param, 'delay')|| strcmp(param, 'delay_dispersion')
            
            % ESPRIT (only delay estimation)
            delta = 1;
            fz = fft(z);
            S = diag(fz);
            S1 = S(1:N-delta,1:N-delta);
            S2 = S(1+delta:N, 1+delta:N);
            Xw = fft(zgn);
            [E,L,~] = svd(Xw');
            es = E(:,1);
            en = E(:,2:end);
            e1 = es(1:N-delta,1);
            e2 = es(1+delta:N,1);
            
            Z_hat = (S2*e1)\(S1*e2);
            lam_hat = eig(Z_hat);
            tau_ESP(ind) = dt*N*atan(imag(lam_hat)/real(lam_hat))/(2*pi);
        
            
            % MUSIC
            ws = 2*pi*horzcat(-linspace(0,N/2,N/2)*Fs/N,linspace(N/2,0,N/2)*Fs/N)';
            D = diag(gradient(fz));
            Pe = en*en';
            fn_v = @(x) exp((x*ws*1.0i));
            fn_G = @(x) [S*fn_v(x),  -1* D*fn_v(x)];
            
            
            fn = @(x) ([1;x(1)]'*real(fn_G(x(2))'*Pe*fn_G(x(2)))*[1;x(1)]) / ([1;x(1)]'*real(fn_G(x(2))'*fn_G(x(2)))*[1;x(1)]);
            x0 = [0.0, tau_ESP(ind)];
            
            xr = fminsearch(fn, x0);
            omega_MUSIC(ind)= -xr(1)/2;
            tau_MUSIC(ind)= xr(2);
        end        
    end    
    
    %% Calculate MSE
    if strcmp(param, 'delay')|| strcmp(param, 'delay_dispersion')
        mse_tau(indx)=mean((tau_est - tau*ones(size(tau_est))).^2);   
        mse_XC1(indx)=mean((tau_XC1 - tau*ones(size(tau_XC1))).^2);
        mse_XC2(indx)=mean((tau_XC2 - tau*ones(size(tau_XC2))).^2);
        mse_ESP(indx)=mean((tau_ESP - tau*ones(size(tau_ESP))).^2);
        mse_MUSIC(indx)=mean((tau_MUSIC - tau*ones(size(tau_MUSIC))).^2);
    end
    if strcmp(param, 'delay')|| strcmp(param, 'delay_dispersion')
        mse_tau_AF_local(indx)=mean((tau_AF_local - tau*ones(size(tau_AF_local))).^2);
        mse_tau_AF(indx)=mean((tau_AF - tau*ones(size(tau_AF))).^2);
    end
    if strcmp(param,'delay_dispersion')
        mse_omega(indx)=mean((omega_est - omega*ones(size(omega_est))).^2); 
        mse_omega_AF_local(indx)=mean((omega_AF_local - omega*ones(size(omega_AF_local))).^2);
        mse_omega_AF(indx)=mean((omega_AF - omega*ones(size(omega_AF))).^2);
        mse_omega_MUSIC(indx)=mean((omega_MUSIC - omega*ones(size(omega_MUSIC))).^2);
    end
    if strcmp(param,'delay_quadtratic')
        mse_tau(indx)=mean((tau_est - tau*ones(size(tau_est))).^2);   
        mse_XC1(indx)=mean((tau_XC1 - tau*ones(size(tau_XC1))).^2);
        mse_XC2(indx)=mean((tau_XC2 - tau*ones(size(tau_XC2))).^2);
        mse_omega(indx)=mean((omega_est - omega*ones(size(omega_est))).^2); 
    end
    
    indx=indx+1;
end


%% MSE plots

if strcmp(param, 'delay')||strcmp(param,'delay_dispersion')    
    % MSE of tau estimate
    figure(3)
    if strcmp(param, 'delay_dispersion')
        ph=semilogy(SNRdb,mse_tau,'k-+',SNRdb,mse_XC2,'k--^',SNRdb,mse_ESP,'k-.d',SNRdb,mse_MUSIC,'k-.X',SNRdb,mse_tau_AF_local,'k-.s',SNRdb,mse_tau_AF,'k--*',SNRdb,CRLB_tau,'k--'); 
        axis([SNRdb(1) SNRdb(end) 10^-8 10])
        lh=legend('CDT','XC','ESPRIT','MUSIC','WBAF (local)','WBAF (global)','CRLB'); set(ph,'LineWidth',2.0, 'MarkerSize',10.0); 
        title('Estimate Time Delay ($$\tau$$), when $$g(t) = \omega t -\tau$$','FontSize',22,'Interpreter','Latex'); 
        xlabel('SNR (dB)','FontSize',22); ylabel('$$mse(\hat{\tau})$$','FontSize',28,'Interpreter','Latex')
    else
        ph=semilogy(SNRdb,mse_tau,'k-+',SNRdb,mse_XC2,'k--^',SNRdb,mse_ESP,'k-.d',SNRdb,mse_tau_AF_local,'k-.s',SNRdb,mse_tau_AF,'k--*',SNRdb,CRLB_tau,'k--'); 
        axis([SNRdb(1) SNRdb(end) 10^-8 10])
        lh=legend('CDT','XC','ESPRIT','MLE (local)','MLE (global)','CRLB'); set(ph,'LineWidth',2.0, 'MarkerSize',10.0); 
        title('Estimate Time Delay ($$\tau$$), when $$g(t) = t -\tau$$','FontSize',22,'Interpreter','Latex'); 
        xlabel('SNR (dB)','FontSize',22); ylabel('$$mse(\hat{\tau})$$','FontSize',28,'Interpreter','Latex')
    end
    set(gca,'FontSize',22,'LineWidth',2.0)
end

if strcmp(param,'delay_dispersion')  
    % MSE of omega estimate
    figure(6)
    ph=semilogy(SNRdb,mse_omega,'k-+',SNRdb,mse_omega_AF_local,'k-.s',SNRdb,mse_omega_AF,'k--*',SNRdb,CRLB_omega,'k--'); 
    axis([SNRdb(1) SNRdb(end) 10^-8 10])
    lh=legend('CDT','WBAF (local)','WBAF (global)','CRLB'); set(ph,'LineWidth',2.0, 'MarkerSize',10.0);
    title('Estimate Dispersion ($$\omega$$), when $$g(t) = \omega t -\tau$$','FontSize',22,'Interpreter','Latex');
    xlabel('SNR (dB)','FontSize',22); ylabel('$$mse(\hat{\omega})$$','FontSize',28,'Interpreter','Latex')
    set(gca,'FontSize',22,'LineWidth',2.0)
end
if strcmp(param,'delay_quadtratic')    
    % MSE of tau estimate
    figure(3)
    ph=semilogy(SNRdb,mse_tau,'k-+',SNRdb,mse_XC2,'k--^',SNRdb,CRLB_tau,'k--'); 
    axis([SNRdb(1) SNRdb(end) 10^-8 10^2])
    lh=legend('CDT','XC','CRLB'); set(ph,'LineWidth',2.0, 'MarkerSize',10.0); 
    set(gca,'FontSize',22,'LineWidth',2.0)
    title('Estimate Time Delay ($$\tau$$), when $$g(t) = \kappa t^2 - \tau$$','FontSize',22,'Interpreter','Latex'); 
    xlabel('SNR','FontSize',20); ylabel('$$mse(\hat{\tau})$$','FontSize',28,'Interpreter','Latex')

    % MSE of omega estimate
    figure(6)
    ph=semilogy(SNRdb,mse_omega,'k-+', SNRdb,CRLB_omega,'k--'); 
    axis([SNRdb(1) SNRdb(end) 10^-8 10^2])
    lh=legend('CDT','CRLB'); set(ph,'LineWidth',2.0, 'MarkerSize',10.0); 
    set(gca,'FontSize',22,'LineWidth',2.0)
    title('Estimate Quadratic Dispersion ($$\kappa$$), when $$g(t) = \kappa t^2 - \tau$$','FontSize',22,'Interpreter','Latex'); 
    xlabel('SNR','FontSize',22); ylabel('$$mse(\hat{\kappa})$$','FontSize',28,'Interpreter','Latex')
end



