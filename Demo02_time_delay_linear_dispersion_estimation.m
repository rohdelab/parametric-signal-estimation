%%
% Demo02::Time delay and linear dispersion estimation using the CDT
% g_p(t) = wt - tau
%
% Reference: Parametric Signal Estimation Using the Cumulative Distribution Transform
%
%%
clc; clear; close all;

N=400;                       % number of points in signal
dt=0.025;                    % timestep
Fs = 1/dt;
eps=1e-12;                   % "small" value for use in CDT estimation

t=-N/2*dt:dt:(N/2-1)*dt;     % time array (0 centered)
tm=0;                        % initial centering of pulse


%% Define original signal (before time delay)
f=1;                            % modulation frequency
win=1;                          % Width of the pulse (all models assume this =1, don't change)

gwin=exp(-(t-tm).^2/(2*win^2)); % Apodization function
z=gwin.*sin(2*pi*f*t);          % The clean input signal
s=z.^2/sum(z.^2);               % Squared, normalized input signal.
Es = mean(z.^2)*(t(end)-t(1));  % Energy of the signal


%% Define signal after time delay+linear dispersion
tau=25.3*dt;                                % True time delay in seconds
omega=0.75;                                 % True linear dispersion parameter
gwin=exp(-(omega*t-tau-tm).^2/(2*win^2));   % Apodization function
zg=gwin.*sin(2*pi*f*(omega*t-tau));         % Altered input signal


%% Define SNR
SNRdb=10;
SNR=10.^(SNRdb/10);
sigma=(pi)^(1/4)*sqrt(win)./(sqrt(2*SNR*N*dt));    % Standard deviation
disp(['SNR: ' num2str(SNRdb) 'dB'])


%% Add noise to the altered (delay+dispersion) signal
noise=sigma*randn(1,N);       % The noise
zgn=zg+noise;                 % The altered clean signal + noise
        

%% CDT estimation process
s0=ones(1,N);                                   % Reference signal 
[shat,df1,xtilde]=CDT(s0,s+eps,t,0,Es);         % CDT of clean "PDF" signal s

r=zgn.^2/sum(zgn.^2);                           % Squared, altered signal + noise
[rhat,df2,xtilde]=CDT(s0,r+eps,t,sigma,Es);     % Noise-corrected CDT of r=sg+noise

clip_CDF=25;                                    % Number of points to clip off beginning and end of CDF estimate (where no signal present)
midrange=clip_CDF:length(rhat)-clip_CDF;        % Restrict domain of CDT to use in estimation.  Has pronounced effect on bias but not variance
   
% Calculate the time delay and linear dispersion using CDTs
B = [rhat(midrange)', -1*ones(length(midrange),1)];
estimate = (B'*B)\(B'*shat(midrange)');
omega_est = estimate(1);
tau_est=estimate(2);


%% Plots
figure;
plot(t,z,'b', 'Linewidth',2.0), hold on
plot(t,zgn,'r', 'Linewidth',2.0)
xlabel('$$t$$','interpreter','latex')
ylabel('PDF')
legend({'$$z(t)$$','$$z_g(t)+\eta(t)$$'},'interpreter','latex')
set(gca,'FontSize',16,'LineWidth',2.0)
title('Input and measured signals')

figure;
plot(t,s,'b', 'Linewidth',2.0), hold on
plot(t,r,'r', 'Linewidth',2.0)
xlabel('$$t$$','interpreter','latex')
ylabel('PDF')
legend({'$$s(t)$$','$$r(t)$$'},'interpreter','latex')
set(gca,'FontSize',16,'LineWidth',2.0)
title('Normalized input PDFs')

figure;
plot(xtilde,shat,'b', 'Linewidth',2.0), hold on
plot(xtilde,rhat,'r', 'Linewidth',2.0)
xlabel('$$t$$','interpreter','latex')
ylabel('CDT')
legend({'$$\hat{s}(t)$$','$$\hat{r}(t)$$'},'interpreter','latex')
set(gca,'FontSize',16,'LineWidth',2.0)
title('CDTs of the signals')


%% Display the results
disp(['True time delay: ' num2str(tau) ' seconds'])
disp(['Estimated time delay: ' num2str(tau_est) ' seconds'])

disp(['True linear dispersion: ' num2str(omega) ' seconds'])
disp(['Estimated linear dispersion: ' num2str(omega_est) ' seconds'])

