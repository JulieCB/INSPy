
% connectivity matrix (number of nodes is determined below, 
% from shape of C)
wizard = @(r) (1-abs(r)).*exp(-abs(r))/4;
N = 20;

% feed forward excitation 1/0
%C = diag(ones(N-1, 1), -1);

% symmetric excitation 1/0
%C = diag(ones(N-1, 1), 1) + diag(ones(N-1, 1), -1);

% mexican hat (w/ inhibition)
C = wizard(abs(meshgrid(0:N) - meshgrid(0:N)')/20);
C = C - diag(diag(C));

% integration parameters
tf = 6000.0;
dt = 0.05;
ds = 10; % downsample factor
Q = ones(2*size(C, 1), 1)*25e-2;

% other model parameters;
parameter_regimes.sparse_spiking = [10 0.65 3 0.4 0.5*0.37];  % narrow spikes
% parameter_regimes.subcritical_hopf = [1.5 0.7 2 0.55 0.37];  % damped, noise driven oscillations
parameter_regimes.subcritical_hopf = [1.5 0.7 2 0.55 0.36];  % damped, noise driven oscillations

par = parameter_regimes.subcritical_hopf;

tau = par(1);    % time scale of cubic nullcline (how fast is it)
b = par(2);      % horizontal position of quadratic nullcline
c = par(3);      % ohrizontal scale of quadratic nullcline
a0 = par(4);     % vertical position of quadratic nullcline (excitability)
g = par(5);      % scaling of connectivity matrix coupling nodes of network

% coupling & frequency scaling
% 0.2 -> ~ 45 Hz, 0.4 -> 90 Hz
f0 = 0.1;

% allocate memory for simulation output
ys = zeros(2*size(C, 1), fix(tf/dt/ds));
ts = zeros(fix(tf/dt/ds), 1);

% initial conditions
y = repmat([-1, 0]', size(C, 1), 1);
y = rand(2*size(C, 1), 1);

% parameter to vary with time
as = zeros(size(ys, 2), 1);

% 1 if x > 0 else 0
H = @(x) (x > 0);

% show animated node amplitudes? (runs a bit slower)
animate = 1;

if animate
    figure(2)
    clf
    h = plot(ys(1:2:end, 1), 'ko-');
    ylim( [-1.5, 1])
    grid
end

r = 0;

% integrate equations
for i=2:fix(tf/dt);
    t = i*dt;
    
    % differential equations
    xt = y(1:2:end);
    yt = y(2:2:end);
    
    dx = tau*(xt - xt.^3 + yt);
    dy = a0 - r - yt - c*(xt + b).^2  + g*C*H(xt);
    
    %r = r + dt*(-r  + mean(H(xt+0.1)))/1000;
    
    % stimulation on first node to put the node into oscillatory regime
    ai =  H(sin(t/1000 * 2 * pi * 1.0) - 0.8)*0.1;  % pulse train
    dy(1:3) = dy(1:3) + ai;
    
    % interleaved variables
    dys = reshape([dx dy]', [], 1);
    
    % Euler update w/ white additive noise
    y = y + f0*dt*(dys + Q.*randn(size(ys, 1), 1));
    
    % update plot
    if animate && mod(i, 350) == 0
        set(h, 'YData', y(1:2:end));
        drawnow
    end
    
    % output solution
    if mod(i, ds) == 0
        ys(:, fix(i/ds)) = y;
        ts(fix(i/ds)) = t;
        as(fix(i/ds)) = ai;
    end
    
end

% filter & downsample 


%% plots

% ms -> s
ts = ts/1000; 

figure(1), clf
subplot (531)
plot(ys(1:2:end, :)', ys(2:2:end, :)')
hold on
title('phase space')
[X, Y] = meshgrid(-1.5:0.1:1, -1.5:0.1:1.5);
DX = -X.^3 + X + Y;
DY = a0 - Y - c*(X + b).^2;
%quiver(X, Y, DX, DY);
contour(X, Y, DX, [0 0]);
contour(X, Y, DY, [0 0]);
grid on
subplot (532)
imagesc(C)
title('connectivity matrix')
colorbar
subplot (533)
for i=1:size(C, 1)
    %y_ = mean(ys(1:2:end, 100:end));
    y_ = ys((i-1)*2+1, 100:end);
    y_ = y_ - mean(y_(:));
    fs = (1/(ts(2)-ts(1))) / 2 * linspace(0, 1, length(y_));
    Fy = abs(fft(y_));
    semilogy(fs(1:fix(length(fs)/2)), Fy(1:fix(length(fs)/2)))
    hold on
end
xlim([1, 100])
grid on

ax{1} = subplot (512);
for i=1:size(C, 1)60
    plot(ts, ys((i-1)*2+1, :) + i, 'k')
    hold on
end
hold off
title('time courses')
ylabel('x')
grid on
ax{2} = subplot(513);
imagesc(ts, 1:N, ys(1:2:end, :))
set(gca, 'CLim', [-1.5, 1])
ax{3} = subplot (514);
plot(ts, as)
hold on
plot(ts, (mean(ys(1:2:end, :)) - mean(mean(ys(1:2:end, :))))*5 + 1)
xlabel('time (s)')
grid on

subplot (515);
ei = find(as(2:end) > mean(as) & as(1:end-1)<mean(as));
ms100 = fix(0.1 / (ts(2) - ts(1)));
te = (ts(2) - ts(1))*(-ms100 : 8*ms100);
event = zeros(length(ei), length(te));
for i=1:length(ei)
    if ei(i)-ms100 > 0 && ei(i)+ms100*8 < size(ys, 2)
        event(i, :) = mean(ys(1:2:end, ei(i)-ms100 : ei(i) + ms100*8));
        %plot(te, event(i, :));
        hold on
    end
end
hold on
evst = std(event);
evmt = mean(event);
evm = mean(event(:));
evs = std(event(:));
plot(te, (evmt - evm)/evs, 'k', 'linewidth', 3)
plot(te, (evmt - evm + evst)/evs, 'k')
plot(te, (evmt - evm - evst)/evs, 'k')
plot(te, as(ei(2)-ms100 : ei(2) + ms100*8), 'r')
%plot(te, 5*mean(exp(1j*(event - evm))));
hold off


% zoom time together
if exist('linkaxes')
    linkaxes([ax{:}], 'x')
end


