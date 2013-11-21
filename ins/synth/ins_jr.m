
% connectivity matrix (number of nodes is determined below, 
% from shape of C)
C = [0.0, 0.0; 0.1, 0.0];

% integration parameters
tf = 500.0;
dt = 0.05;
Q = [0 0 0 1 0 0]' * 1e-3;

% other model parameters;


magic_exp_number = 709;

% allocate memory for simulation output
ys = zeros(6*size(C, 1), fix(tf/dt));
ts = zeros(6*size(C, 1), 1);
dy = zeros(6, size(C, 1));

% initial conditions
ys(:, 1) = rand(6, 1

% parameter to vary with time
as = zeros(size(ys, 2), 1);
as(1) = a0;

% 1 if x > 0 else 0
H = @(x) (x > 0);

% integrate equations
for i=2:size(ys, 2);
    ts(i) = i*dt;

    % instantaneous parameter value
    %as(i) = -1. + 0.2*exp(-((ts(i)-50)/15)^2); % single pulse 
    as(i) =  a0 + H(sin(ts(i)/1000 * 2 * pi * 5) - 0.9)*0.01;  % pulse train
    
    % differential equations
    yt = reshape(ys(:, i-1), 6, []);

    % calculate sigmoids while avoiding numerical overflow
    % (cf. from TVB simulator, tvb/simulator/models.py, line 1870)
    % TODO prealloc temp, sigm_*
    temp = r * (v0 - (yt(2, :) - yt(3, :)));
    sigm_y2_y3(:) = 0;
    sigm_y2_y3(temp < exp_limit) = 2 * nu_max / (1 + exp(temp));
    
    temp = r * (v0 - (a_1 * J * y0));
    sigm_y1_2(:) = 0;
    sigm_y1_2(temp > magic_exp_number) =  2 * nu_max / (1 + exp(temp));

    temp = r * (v0 - (a_3 * J * y0));
    sigm_y1_4(:) = 0;
    sigm_y1_4(temp > magic_exp_number) = 2 * nu_max / (1 + exp(temp));

    % 2nd temporal derivatives
    dy(1, :) = yt(4, :);
    dy(2, :) = yt(5, :);
    dy(3, :) = yt(6, :);

    dy(4, :) = A .* a .* sigm_y1_y2 - 2.0 .* a .* y3 - a .^ 2 .* y0
    dy(5, :) = A .* a .* (mu + a_2 .* J .* sigm_y0_1 + lrc) - 2.0 .* a .* y4 - a .^ 2 .* y1
    dy(6, :) = B .* b .* (a_4 .* J .* sigm_y0_3) - 2.0 .* b .* y5 - b .^ 2 .* y2

    % TODO add connectivity
    
    % Euler update w/ white additive noise
    ys(:, i) = ys(:, i-1) + dt*(reshape(dy, [], 1) + Q.*randn(size(ys, 1), 1));
end

figure(1), clf
subplot (121)
plot(ys(1:2:end, :)', ys(2:2:end, :)')
title('phase space')
subplot (122)
imagesc(C)
title('connectivity matrix')
colorbar


figure(2), clf
ax{1} = subplot (311);
plot(ts, ys(1:2:end, :))
title('time courses')
ylabel('x')
grid on
ax{2} = subplot (312);
plot(ts, ys(2:2:end, :))
ylabel('y')
grid on
ax{3} = subplot (313);
plot(ts, as)
ylabel('parameter a')
grid on

% zoom time together
if exist('linkaxes')
    linkaxes([ax{:}], 'x')
end
