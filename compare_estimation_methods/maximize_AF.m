function xr = maximize_AF(zgn, t, tm, f, win, x0)

% rng default
if length(x0)>1
    fn_zg = @(x) exp(-(x(1)*t-x(2)-tm).^2/(2*win^2)).*sin(2*pi*f*(x(1)*t-x(2)));
%     fn = @(x) -sum(zgn.*fn_zg(x));    % negative of WBAF

%     fn_zg = @(x) exp(-(x(1)*(t-x(2))-tm).^2/(2*win^2)).*sin(2*pi*f*(x(1)*(t-x(2))));
%     fn = @(x) -sqrt(abs(x(1)))*sum(zgn.*fn_zg(x));    % negative of WBAF
    fn = @(x) -abs(sqrt(abs(x(1)))*sum(zgn.*fn_zg(x))).^2;

    lb = [0.001, 0];
    ub = [20, t(end)];
else
    fn_zg = @(x) exp(-(t-x-tm).^2/(2*win^2)).*sin(2*pi*f*(t-x));
    %fn = @(x) -abs(sum(zgn.*fn_zg(x))).^2;
    fn = @(x) -(sum(zgn.*fn_zg(x)));
    lb = 0;
    ub = t(end);
end

% rng default % For reproducibility
opts = optimoptions(@fmincon,'Algorithm','sqp');
%opts = optimoptions(@fmincon);
problem = createOptimProblem('fmincon','objective',...
    fn,'x0',x0,'lb',lb,'ub',ub,'options',opts);   % lb-> lower bounds, ub-> upper bounds
                                                         % bounds: w>=0
gs = GlobalSearch;
[xr,f] = run(gs,problem);

% if length(xr)>1
%     xr(2) = xr(1)*xr(2);
% end
