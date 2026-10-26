clear;
clc;
close all;


% ML parameters
gCa = 4.4;
gK  = 8.0;
gL  = 2.0;

VCa = 120;
VK  = -84;
VL  = -60;

phi = 0.02;

V1 = -1.2;
V2 = 18;
V3 = 2;
V4 = 30;
V5 = 2;
V6 = 30;

C = 20;

Iext = 0;

p.gCa = gCa;
p.gK = gK;
p.gL = gL;

p.VCa = VCa;
p.VK = VK;
p.VL = VL;

p.phi = phi;

p.V1 = V1;
p.V2 = V2;
p.V3 = V3;
p.V4 = V4;
p.V5 = V5;
p.V6 = V6;

p.C = C;
p.Iext = Iext;

% Method -1 - Reducing the system to one variable
m_inf = @(V) 0.5*(1+ tanh((V-V1)/V2));
w_inf = @(V) 0.5*(1+ tanh((V-V3)/V4));

F = @(V) Iext - gCa*m_inf(V)*(V-VCa) - gK*w_inf(V)*(V-VK) - gL*(V-VL) ;

% Equilibrium state

V_eq = fzero(F,-60);   %starting voltage is chosen to be leark reversal potential
w_eq = w_inf(V_eq);


fprintf('M-1: Equilibrium:\n');
fprintf('V* = %.6f mV\n',V_eq);
fprintf('w* = %.6f\n', w_eq);

% Method- 2 Solving both eqution silmultaneously

tau_w = @(V) 1 ./ cosh((V-V5)/V6);
F2 = @(x) [(Iext - gCa*m_inf(x(1))*(x(1)-VCa) - gK*x(2)*(x(1)-VK) - gL*(x(1)-VL))/C ; ...
    phi*(w_inf(x(1)) - x(2)) / tau_w(x(1)) ] ;

x0 = [-60; 0.02];

x_eq = fsolve(F2,x0);
fprintf('M-2: Equilibrium:\n');
fprintf('V* = %.6f mV\n', x_eq(1));
fprintf('w* = %.6f\n', x_eq(2));

% plots
figure;

% Grid for phase plane
Vgrid = linspace(-80,50,60);
wgrid = linspace(0,0.5,45);

[V,W] = meshgrid(Vgrid,wgrid);

% Vector field
dV = (Iext -gCa*m_inf(V).*(V-VCa) - gK*W.*(V-VK) - gL*(V-VL))/C;

dW = phi*(w_inf(V)-W)./tau_w(V);

% Convert vector field to plotted coordinates (V,100w)
dV_plot = dV;
dW_plot = 100*dW;

% Normalize
speed = sqrt(dV_plot.^2 + dW_plot.^2);

dVn = dV_plot ./ (speed + eps);
dWn = dW_plot ./ (speed + eps);

% Plot vector field
quiver(V,100*W,dVn,dWn,0.7,'LineWidth',1);

hold on;

% Nullclines
Vnull = linspace(-80,50,1000);

W_w = w_inf(Vnull);

W_v = (Iext - gCa*m_inf(Vnull).*(Vnull-VCa)- gL*(Vnull-VL))./ (gK*(Vnull-VK));

plot(Vnull,100*W_w,'LineWidth',2);
plot(Vnull,100*W_v,'LineWidth',2);

% Equilibrium
plot(V_eq,100*w_eq,'ko','MarkerFaceColor','k','MarkerSize',8);
%{
%Trajectory;
x0 = [V_eq + 50; w_eq];
x01 = [V_eq + 42; w_eq];

[t,x] = ode45(@(t,x) morrisLecar(t,x,p),[0 300],x0);
[t1,x1] = ode45(@(t,x) morrisLecar(t,x,p),[0 300],x01);

plot(x(:,1),100*x(:,2),'k','LineWidth',1.5,'Color','black');
plot(x1(:,1),100*x1(:,2),'k','LineWidth',1.5,'Color','blue');
%}

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Phase Plane');

xlim([-80 50]);
ylim([0 50]);

legend('Vector field','w-nullcline','V-nullcline', 'Equilibrium','Location','northeast');

grid on;

% --------------------------------------------------------
% Q3 - Jacobian and eigenvalues
% --------------------------------------------------------

syms V w

% Symbolic steady-state functions
m_inf_sym = 0.5*(1 + tanh((V - V1)/V2));
w_inf_sym = 0.5*(1 + tanh((V - V3)/V4));
tau_w_sym = 1/cosh((V - V5)/V6);

% MLE equations
f_sym = (Iext ...
    - gCa*m_inf_sym*(V - VCa) ...
    - gK*w*(V - VK) ...
    - gL*(V - VL))/C;

g_sym = phi*(w_inf_sym - w)/tau_w_sym;

% Jacobian
J = jacobian([f_sym; g_sym],[V w]);

% Evaluate at equilibrium
J_eq = double(subs(J,{V,w},{V_eq,w_eq}));

fprintf('\nJacobian at equilibrium:\n');
disp(J_eq);

% Eigenvalues
lambda = eig(J_eq);

fprintf('Eigenvalues:\n');
disp(lambda);


% --------------------------------------------------------
% Q5 - Effect of phi on action potential
% --------------------------------------------------------

phis = [0.01,0.02 0.04];
p.Iext = 0;
figure;
hold on;

for i = 1:length(phis)

    p.phi = phis(i);

    % Trigger action potential by depolarizing V
    x0 = [V_eq + 50; w_eq];

    [t,x] = ode45(@(t,x) morrisLecar(t,x,p), ...
                  [0 300],x0);

    plot(t,x(:,1),'LineWidth',1.5);

end

xlabel('Time (ms)');
ylabel('V (mV)');
title('Morris-Lecar Action Potentials');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;


% Phase-plane plot
figure;
hold on;

for i = 1:length(phis)

    p.phi = phis(i);

    x0 = [V_eq + 50; w_eq];

    [~,x] = ode45(@(t,x) morrisLecar(t,x,p), ...
                  [0 300],x0);

    plot(x(:,1),100*x(:,2),'LineWidth',1.5);

end


plot(Vnull,100*W_w,'LineWidth',2);
plot(Vnull,100*W_v,'LineWidth',2);

% Equilibrium
plot(V_eq,100*w_eq,'ko','MarkerFaceColor','k','MarkerSize',8);
%IC
plot(V_eq+50,100*w_eq,'ko','MarkerFaceColor','k','MarkerSize',7);


xlim([-80 50]);
ylim([0 50]);

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Phase Plane');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04','w-nullcline','V-nullcline', 'Equilibrium','IC');
grid on;
hold off;


% --------------------------------------------------------
% Q6 - Threshold depolarization
% --------------------------------------------------------

dV0 = 0:1:80;
p.Iext = 0;
p.phi = 0.02;
V0_values = V_eq + dV0;

maxV = zeros(size(V0_values));

for i = 1:length(V0_values)

    x0 = [V0_values(i); w_eq];

    [t,x] = ode45(@(t,x) morrisLecar(t,x,p), ...
                  [0 300],x0);

    maxV(i) = max(x(:,1));

end

figure;

plot(dV0,maxV,'LineWidth',1.5);

xlabel('Initial depolarization \DeltaV_0 (mV)');
ylabel('Maximum V (mV)');
title('Maximum Voltage vs Initial Depolarization');
grid on;

low = V_eq + 46; %from plot
high = V_eq + 49;

for k = 1:25

    mid = (low + high)/2;

    x0 = [mid; w_eq];

    [~,x] = ode45(@(t,x) morrisLecar(t,x,p), ...
                  [0 300],x0);

    maxV_mid = max(x(:,1));

    if maxV_mid > 0
        high = mid;
    else
        low = mid;
    end

end

V_threshold = (low + high)/2;
deltaV_threshold = V_threshold - V_eq;

fprintf('\nThreshold:\n');
fprintf('V_threshold = %.10f mV\n',V_threshold);
fprintf('Delta V_threshold = %.10f mV\n',deltaV_threshold);


V0_sub = V_threshold - 0.1;
V0_sup = V_threshold + 0.1;

x0_sub = [V0_sub; w_eq];
x0_sup = [V0_sup; w_eq];

[t_sub,x_sub] = ode45(@(t,x) morrisLecar(t,x,p), ...
                     [0 300],x0_sub);

[t_sup,x_sup] = ode45(@(t,x) morrisLecar(t,x,p), ...
                     [0 300],x0_sup);

figure;

plot(x_sub(:,1),100*x_sub(:,2),'LineWidth',1.5);
hold on;

plot(x_sup(:,1),100*x_sup(:,2),'LineWidth',1.5);

plot(V_eq,100*w_eq,'ko', ...
     'MarkerFaceColor','k');

plot(V0_sub,100*w_eq,'ks', ...
     'MarkerFaceColor','k');

plot(V0_sup,100*w_eq,'ks', ...
     'MarkerFaceColor','k');

plot(Vnull,100*W_w,'LineWidth',2,'Color','black');
plot(Vnull,100*W_v,'LineWidth',2,'Color','black');

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Threshold Trajectories');

xlim([-80 50]);
ylim([0 50]);


legend('Subthreshold','Suprathreshold', ...
       'Equilibrium','Initial conditions','w-nullcline','V-nullcline');

grid on;


% --------------------------------------------------------
% Q7 - MLE at Iext = 86
% --------------------------------------------------------

p.Iext = 86;
p.phi = 0.02;

% Find equilibrium for Iext = 86
F86 = @(x) [ ...
    (p.Iext ...
    - p.gCa*m_inf(x(1))*(x(1)-p.VCa) ...
    - p.gK*x(2)*(x(1)-p.VK) ...
    - p.gL*(x(1)-p.VL))/p.C; ...

    p.phi*(w_inf(x(1))-x(2))/tau_w(x(1))
    ];

x_eq86 = fsolve(F86,[-28;0.12]);

fprintf('\nQ7: Equilibrium at Iext = 86\n');
fprintf('V* = %.6f mV\n',x_eq86(1));
fprintf('w* = %.6f\n',x_eq86(2));


% Trial 1: old equilibrium from Iext = 0
x01 = [V_eq; w_eq];

% Trial 2: equilibrium for Iext = 86
x02 = x_eq86;

% Trial 3: specified off-equilibrium point
x03 = [-27.9;0.17];


% Simulate
[t1,x1] = ode45(@(t,x) morrisLecar(t,x,p), ...
                [0 500],x01);

[t2,x2] = ode45(@(t,x) morrisLecar(t,x,p), ...
                [0 500],x02);

[t3,x3] = ode45(@(t,x) morrisLecar(t,x,p), ...
                [0 500],x03);


% Phase plane
figure;
hold on;

plot(x1(:,1),100*x1(:,2),'LineWidth',1.5);
plot(x2(:,1),100*x2(:,2),'LineWidth',1.5);
plot(x3(:,1),100*x3(:,2),'LineWidth',1.5);

% Mark the I=86 equilibrium
plot(x_eq86(1),100*x_eq86(2),'ko', ...
    'MarkerFaceColor','k','MarkerSize',8);
plot(x01(1),100*x01(2),'kd', ...
    'MarkerFaceColor','k','MarkerSize',8);
plot(x03(1),100*x03(2),'ks', ...
    'MarkerFaceColor','k','MarkerSize',8);

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar at I_{ext} = 86 \muA/cm^2');

legend('Trial 1: I=0 equilibrium', ...
       'Trial 2: I=86 equilibrium', ...
       'Trial 3: (-27.9,0.17)', ...
       'I=86 equilibrium',...
       'Trial 1 IC',...
       'Trial 3 IC');

grid on;
hold off;


syms V w

% Symbolic steady-state functions
m_inf_sym = 0.5*(1 + tanh((V - V1)/V2));
w_inf_sym = 0.5*(1 + tanh((V - V3)/V4));
tau_w_sym = 1/cosh((V - V5)/V6);

% MLE equations
f_sym = (p.Iext ...
    - gCa*m_inf_sym*(V - VCa) ...
    - gK*w*(V - VK) ...
    - gL*(V - VL))/C;

g_sym = p.phi*(w_inf_sym - w)/tau_w_sym;

% Jacobian
J = jacobian([f_sym; g_sym],[V w]);

% Evaluate at equilibrium
J_eq = double(subs(J,{V,w},{x02(1),x02(2)}));

fprintf('\nJacobian at new equilibrium (I=86 mv/cm2):\n');
disp(J_eq);

% Eigenvalues
lambda = eig(J_eq);

fprintf('Eigenvalues:\n');
disp(lambda);


% --------------------------------------------------------
% Q9 - Equilibrium and stability for I = 80, 86, 90
% --------------------------------------------------------

currents = [80 86 90];

eq_values = zeros(length(currents),2);
eig_values = zeros(length(currents),2);

for i = 1:length(currents)

    I = currents(i);

    % Equilibrium equations
    F_I = @(x) [ ...
        (I ...
        - gCa*m_inf(x(1))*(x(1)-VCa) ...
        - gK*x(2)*(x(1)-VK) ...
        - gL*(x(1)-VL))/C;

        phi*(w_inf(x(1))-x(2))/tau_w(x(1))
        ];

    % Find equilibrium
    xeq = fsolve(F_I,[-28;0.12]);

    eq_values(i,:) = xeq(:)';

    % Jacobian numerically using symbolic Jacobian J
    J_I = double(subs(J,{V,w},{xeq(1),xeq(2)}));

    lam = eig(J_I);

    eig_values(i,:) = lam(:)';

    fprintf('\nIext = %.1f\n',I);
    fprintf('V* = %.6f mV\n',xeq(1));
    fprintf('w* = %.6f\n',xeq(2));
    fprintf('Eigenvalues:\n');
    disp(lam);

end

currents = 80:100;
firing_rate = zeros(size(currents));
for i = 1:length(currents)

    I = currents(i);

    p.Iext = I;
    p.phi = 0.02;

    % Find equilibrium
    F_I = @(x) [ ...
        (I ...
        - gCa*m_inf(x(1))*(x(1)-VCa) ...
        - gK*x(2)*(x(1)-VK) ...
        - gL*(x(1)-VL))/C;

        phi*(w_inf(x(1))-x(2))/tau_w(x(1))
        ];

    xeq = fsolve(F_I,[-28;0.12]);

    % Start at equilibrium
    x0 = xeq;

    [t,x] = ode45(@(t,x) morrisLecar(t,x,p), ...
                  [0 3000],x0);

    Vt = x(:,1);

    % Upward zero crossings
    idx = find(Vt(1:end-1)<0 & Vt(2:end)>=0);

    % Ignore initial transient
    valid = t(idx) > 1000;

    spike_times = t(idx(valid));

    if length(spike_times) >= 2

        firing_rate(i) = ...
            (length(spike_times)-1) / ...
            ((spike_times(end)-spike_times(1))/1000);

    else

        firing_rate(i) = 0;

    end

end

figure;

plot(currents,firing_rate,'o-','LineWidth',1.5);

xlabel('Applied current I_{ext} (\muA/cm^2)');
ylabel('Firing rate (Hz)');
title('Morris-Lecar Firing Rate vs Applied Current');

grid on;

% --------------------------------------------------------
% Q10 - Morris-Lecar: 3 Equilibria at Iext = 30
% --------------------------------------------------------

% Q10 parameters
gCa = 4.0;
gK  = 8.0;
gL  = 2.0;

VCa = 120;
VK  = -84;
VL  = -60;

phi = 0.0667;

V1 = -1.2;
V2 = 18;

V3 = 12;
V4 = 17.4;

V5 = 12;
V6 = 17.4;

C = 20;

Iext = 30;


% Update parameter structure used by morrisLecar()
p.gCa = gCa;
p.gK  = gK;
p.gL  = gL;

p.VCa = VCa;
p.VK  = VK;
p.VL  = VL;

p.phi = phi;

p.V1 = V1;
p.V2 = V2;
p.V3 = V3;
p.V4 = V4;
p.V5 = V5;
p.V6 = V6;

p.C = C;
p.Iext = Iext;


% Steady-state functions
m_inf = @(V) 0.5*(1 + tanh((V-V1)/V2));
w_inf = @(V) 0.5*(1 + tanh((V-V3)/V4));

tau_w = @(V) 1 ./ cosh((V-V5)/V6);


% Equilibrium equations
F10 = @(x) [ ...
    (Iext ...
    - gCa*m_inf(x(1))*(x(1)-VCa) ...
    - gK*x(2)*(x(1)-VK) ...
    - gL*(x(1)-VL))/C;

    phi*(w_inf(x(1))-x(2))/tau_w(x(1))
    ];


% --------------------------------------------------------
% Find all 3 equilibrium points
% --------------------------------------------------------

% Different initial guesses are required because
% fsolve may converge to different equilibria.
guesses = [
    -50, 0.01;
    -10, 0.03;
      10, 0.30
];

eq10 = zeros(3,2);

for i = 1:3
    eq10(i,:) = fsolve(F10,guesses(i,:)).';
end


fprintf('\n====================================\n');
fprintf('Q10: Equilibria for Iext = 30\n');
fprintf('====================================\n');

for i = 1:3
    fprintf('Equilibrium %d: V = %.6f mV, w = %.6f\n', ...
        i,eq10(i,1),eq10(i,2));
end



% --------------------------------------------------------
% Jacobian and stability of the 3 equilibria
% --------------------------------------------------------

syms Vs ws

m_inf_sym = 0.5*(1 + tanh((Vs-V1)/V2));
w_inf_sym = 0.5*(1 + tanh((Vs-V3)/V4));
tau_w_sym = 1/cosh((Vs-V5)/V6);

f_sym = (Iext ...
    - gCa*m_inf_sym*(Vs-VCa) ...
    - gK*ws*(Vs-VK) ...
    - gL*(Vs-VL))/C;

g_sym = phi*(w_inf_sym-ws)/tau_w_sym;

J10 = jacobian([f_sym;g_sym],[Vs ws]);


fprintf('\nEigenvalues:\n');

for i = 1:3

    Vstar = eq10(i,1);
    wstar = eq10(i,2);

    Jstar = double(subs(J10,{Vs,ws},{Vstar,wstar}));

    lambda = eig(Jstar);

    fprintf('\nEquilibrium %d:\n',i);
    fprintf('Jacobian:\n');
    disp(Jstar);

    fprintf('Eigenvalues:\n');
    disp(lambda);

end

% --------------------------------------------------------
% Q10 Phase Plane
% --------------------------------------------------------

Vnull = linspace(-80,50,1500);

W_w = w_inf(Vnull);

W_v = (Iext ...
    - gCa*m_inf(Vnull).*(Vnull-VCa) ...
    - gL*(Vnull-VL)) ...
    ./ (gK*(Vnull-VK));


figure;
hold on;

plot(Vnull,100*W_w,'LineWidth',2);
plot(Vnull,100*W_v,'LineWidth',2);

% Plot all three equilibria
plot(eq10(:,1),100*eq10(:,2),'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Phase Plane: I_{ext}=30 \muA/cm^2');

xlim([-80 50]);
ylim([0 50]);

legend('w-nullcline', ...
       'V-nullcline', ...
       'Equilibria', ...
       'Location','best');

grid on;
hold off;


% --------------------------------------------------------
% Stable and unstable manifolds of saddle
% --------------------------------------------------------

V_s = eq10(2,1);
w_s = eq10(2,2);

% Jacobian at saddle
Js = double(subs(J10,{Vs,ws},{V_s,w_s}));

% Eigenvalues and eigenvectors
[vecs,vals] = eig(Js);

lambda_s = diag(vals);

% Identify stable and unstable directions
stable_idx   = find(real(lambda_s) < 0,1);
unstable_idx = find(real(lambda_s) > 0,1);

vs = real(vecs(:,stable_idx));
vu = real(vecs(:,unstable_idx));

% Normalize
vs = vs/norm(vs);
vu = vu/norm(vu);

% Small displacement
epsilon = 1e-6;

saddle = [V_s; w_s];

% Initial points on the four branches
xs1 = saddle + epsilon*vs;
xs2 = saddle - epsilon*vs;

xu1 = saddle + epsilon*vu;
xu2 = saddle - epsilon*vu;

% Stop when trajectory leaves plotting region
options = odeset( ...
    'Events', @(t,x) manifoldEvent(t,x) );

% --------------------------------------------------------
% Stable manifold
% Integrate backward in time
% --------------------------------------------------------

[t_s1,x_s1] = ode45(@(t,x) morrisLecar(t,x,p), ...
                    [0 -500], xs1, options);

[t_s2,x_s2] = ode45(@(t,x) morrisLecar(t,x,p), ...
                    [0 -500], xs2, options);

% --------------------------------------------------------
% Unstable manifold
% Integrate forward in time
% --------------------------------------------------------

[t_u1,x_u1] = ode45(@(t,x) morrisLecar(t,x,p), ...
                    [0 500], xu1, options);

[t_u2,x_u2] = ode45(@(t,x) morrisLecar(t,x,p), ...
                    [0 500], xu2, options);

% --------------------------------------------------------
% Q10 Phase Plane with Manifolds
% --------------------------------------------------------

Vnull = linspace(-80,50,1500);

W_w = w_inf(Vnull);

W_v = (Iext ...
    - gCa*m_inf(Vnull).*(Vnull-VCa) ...
    - gL*(Vnull-VL)) ...
    ./ (gK*(Vnull-VK));

figure;
hold on;

% Nullclines
plot(Vnull,100*W_w,'LineWidth',2);
plot(Vnull,100*W_v,'LineWidth',2);

% Equilibria
plot(eq10(:,1),100*eq10(:,2),'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

% Stable manifold (Blue Dotted)
plot(x_s1(:,1), 100*x_s1(:,2), 'b--', 'LineWidth', 2);
plot(x_s2(:,1), 100*x_s2(:,2), 'b--', 'LineWidth', 2);

% Unstable manifold (Green Dotted)
plot(x_u1(:,1), 100*x_u1(:,2), 'g:', 'LineWidth', 2);
plot(x_u2(:,1), 100*x_u2(:,2), 'g:', 'LineWidth', 2);

% Saddle point highlighted
plot(V_s,100*w_s,'ks', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',10);

xlabel('V (mV)');
ylabel('100w');

title('Morris-Lecar Phase Plane: I_{ext}=30 \muA/cm^2');

xlim([-80 50]);
ylim([0 50]);

legend('w-nullcline', ...
       'V-nullcline', ...
       'Equilibria', ...
       'Stable manifold', ...
       '', ...
       'Unstable manifold', ...
       '', ...
       'Saddle point', ...
       'Location','best');

grid on;
hold off;

% --------------------------------------------------------
% Q11 - Equilibrium continuation from I = 30 to 50
% --------------------------------------------------------

% Q10 parameters remain active
p.Iext = 30;
p.phi = 0.0667;

% Equilibrium current as a function of V
Ieq = @(V) ...
    gCa*m_inf(V).*(V-VCa) ...
    + gK*w_inf(V).*(V-VK) ...
    + gL*(V-VL);

% Current range
currents = 30:0.05:50;

% Storage
eqV = cell(size(currents));
eqW = cell(size(currents));
eqType = cell(size(currents));

Vscan = linspace(-80,50,5000);

for k = 1:length(currents)

    I = currents(k);

    F = @(V) Ieq(V) - I;

    Fscan = F(Vscan);

    roots = [];

    % Find sign changes
    for j = 1:length(Vscan)-1

        if Fscan(j)*Fscan(j+1) < 0

            root = fzero(F,[Vscan(j),Vscan(j+1)]);

            % Avoid duplicate roots
            if isempty(roots) || abs(root-roots(end)) > 1e-4
                roots(end+1) = root;
            end

        end

    end

    eqV{k} = roots;
    eqW{k} = w_inf(roots);

end

% --------------------------------------------------------
% Stability classification
% --------------------------------------------------------

for k = 1:length(currents)

    roots = eqV{k};
    wroots = eqW{k};

    types = cell(size(roots));

    for j = 1:length(roots)

        Vstar = roots(j);
        wstar = wroots(j);

        Jstar = double(subs(J10,{Vs,ws},{Vstar,wstar}));

        lambda = eig(Jstar);

        r = real(lambda);

        if all(r < 0)

            if abs(imag(lambda(1))) > 1e-8
                types{j} = 'Stable focus';
            else
                types{j} = 'Stable node';
            end

        elseif any(r > 0) && any(r < 0)

            types{j} = 'Saddle';

        elseif all(r > 0)

            if abs(imag(lambda(1))) > 1e-8
                types{j} = 'Unstable focus';
            else
                types{j} = 'Unstable node';
            end

        else
            types{j} = 'Critical';

        end

    end

    eqType{k} = types;

end

figure;
hold on;

for k = 1:length(currents)

    roots = eqV{k};
    types = eqType{k};

    for j = 1:length(roots)

        if strcmp(types{j},'Stable node')
            plot(currents(k),roots(j),'bo');

        elseif strcmp(types{j},'Saddle')
            plot(currents(k),roots(j),'ks');

        elseif strcmp(types{j},'Unstable focus')
            plot(currents(k),roots(j),'rx');

        elseif strcmp(types{j},'Stable focus')
            plot(currents(k),roots(j),'g.');

        end

    end

end

xlabel('I_{ext} (\muA/cm^2)');
ylabel('Equilibrium V^* (mV)');
title('Morris-Lecar Equilibrium Branches');
legend('Stable node','Saddle','Unstable focus');
grid on;


figure;
hold on;

for k = 1:length(currents)

    if currents(k) < 39 || currents(k) > 40.5
        continue;
    end

    roots = eqV{k};
    types = eqType{k};

    for j = 1:length(roots)

        if strcmp(types{j},'Stable node')
            plot(currents(k),roots(j),'bo');

        elseif strcmp(types{j},'Saddle')
            plot(currents(k),roots(j),'ks');

        elseif strcmp(types{j},'Unstable focus')
            plot(currents(k),roots(j),'rx');

        end

    end

end

xlabel('I_{ext} (\muA/cm^2)');
ylabel('Equilibrium V^* (mV)');
title('Equilibrium Behavior Near 39-40 \muA/cm^2');

xlim([39 40]);
grid on;

% --------------------------------------------------------
% Firing rate from I = 30 to 45
% --------------------------------------------------------

currents_FR = 30:0.1:45;
firing_rate = zeros(size(currents_FR));

for k = 1:length(currents_FR)

    I = currents_FR(k);

    p.Iext = I;
    p.phi = 0.0667;

    % Find all equilibria
    F = @(V) Ieq(V) - I;

    Fscan = F(Vscan);
    roots = [];

    for j = 1:length(Vscan)-1

        if Fscan(j)*Fscan(j+1) < 0

            root = fzero(F,[Vscan(j),Vscan(j+1)]);

            if isempty(roots) || abs(root-roots(end)) > 1e-4
                roots(end+1) = root;
            end

        end

    end

    % Use the highest-voltage equilibrium
    % (this is the equilibrium that becomes unstable)
    Veq = roots(end);
    weq = w_inf(Veq);

    x0 = [Veq;weq];

    % After the saddle-node, the equilibrium is unstable.
    % Add a tiny perturbation so numerical integration can
    % leave the unstable equilibrium.
    if I > 39.963
        x0 = x0 + [1e-6;0];
    end

    [t,x] = ode45(@(t,x)morrisLecar(t,x,p), ...
                  [0 5000],x0);

    Vt = x(:,1);

    % Remove transient
    idx_time = t > 2500;

    t_ss = t(idx_time);
    V_ss = Vt(idx_time);

    % Detect upward crossings of 0 mV
    idx = find(V_ss(1:end-1)<0 & V_ss(2:end)>=0);

    spike_times = t_ss(idx);

    if length(spike_times)>=2

        period = mean(diff(spike_times));

        firing_rate(k) = 1000/period;

    else

        firing_rate(k) = 0;

    end

end

figure;

plot(currents_FR,firing_rate,'o-','LineWidth',1.5);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('Firing rate (Hz)');
title('Morris-Lecar Firing Rate');
grid on;

clear;

% --------------------------------------------------------
% Q12 - Hodgkin-Huxley Model
% --------------------------------------------------------

% Parameters
p.GK  = 36;
p.EK  = -72;

p.GNa = 120;
p.ENa = 55;

p.GL  = 0.3;
p.EL  = -60;      % Temporary; Q13 will determine this

p.C   = 1;

% T = 6.3 deg C -> phi = 1
p.phi = 1;

p.Iext = 0;


% --------------------------------------------------------
% Initial condition
% --------------------------------------------------------

V0 = -60;

[an,bn,am,bm,ah,bh] = hh_rates(V0);

n0 = an/(an+bn);
m0 = am/(am+bm);
h0 = ah/(ah+bh);

X0 = [V0;n0;m0;h0];


% --------------------------------------------------------
% Simulation
% --------------------------------------------------------

options = odeset('RelTol',1e-6,'AbsTol',1e-9);

[t,X] = ode45(@(t,X) hodgkinHuxley(t,X,p), ...
              [0 100],X0,options);

fprintf('Number of time points = %d\n',length(t));
fprintf('Initial V = %.6f\n',X(1,1));
fprintf('Final V   = %.6f\n',X(end,1));

disp('Size of t:');
disp(size(t));

disp('Size of X:');
disp(size(X));

disp('First few t values:');
disp(t(1:min(5,end)));

disp('First few X values:');
disp(X(1:min(5,end),:));

disp('Any NaN?');
disp(any(isnan(X(:))));

disp('Any Inf?');
disp(any(isinf(X(:))));


% --------------------------------------------------------
% Plot voltage
% --------------------------------------------------------

figure;

plot(t,X(:,1),'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('V (mV)');
title('Hodgkin-Huxley Membrane Potential');

grid on;


% --------------------------------------------------------
% Plot gating variables
% --------------------------------------------------------

figure;

plot(t,X(:,2),'LineWidth',1.2);
hold on;
plot(t,X(:,3),'LineWidth',1.2);
plot(t,X(:,4),'LineWidth',1.2);

xlabel('Time (ms)');
ylabel('Gating variables');

title('Hodgkin-Huxley Gating Variables');

legend('n','m','h');

grid on;

% --------------------------------------------------------
% Q13 - Determine EL and verify resting potential
% --------------------------------------------------------

Vrest = -60;

% Calculate gating variables at rest
[an,bn,am,bm,ah,bh] = hh_rates(Vrest);

n_rest = an/(an + bn);
m_rest = am/(am + bm);
h_rest = ah/(ah + bh);

% Potassium current
IK_rest = p.GK*n_rest^4*(Vrest - p.EK);

% Sodium current
INa_rest = p.GNa*m_rest^3*h_rest*(Vrest - p.ENa);

% Determine leak reversal potential
p.EL = Vrest + (IK_rest + INa_rest)/p.GL;

% No applied current
p.Iext = 0;

fprintf('\n====================================\n');
fprintf('Q13 - HH Resting State\n');
fprintf('====================================\n');

fprintf('n_rest = %.8f\n',n_rest);
fprintf('m_rest = %.8f\n',m_rest);
fprintf('h_rest = %.8f\n',h_rest);

fprintf('IK_rest  = %.8f uA/cm^2\n',IK_rest);
fprintf('INa_rest = %.8f uA/cm^2\n',INa_rest);

fprintf('EL = %.8f mV\n',p.EL);


% Initial condition at rest
X0 = [Vrest;
      n_rest;
      m_rest;
      h_rest];

% Required tolerances
options = odeset('AbsTol',1e-6, ...
                 'RelTol',1e-3);

% Stiff solver
[t,X] = ode15s(@(t,X)hodgkinHuxley(t,X,p), ...
               [0 100],X0,options);


% Plot voltage
figure;

plot(t,X(:,1),'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('V (mV)');
title('HH Resting Potential Verification');

grid on;

% --------------------------------------------------------
% Q14 - Stability of HH resting state
% --------------------------------------------------------

% Make sure Q13 parameters are being used
p.Iext = 0;

% Resting state from Q13
Vrest = -60;

[an,bn,am,bm,ah,bh] = hh_rates(Vrest);

n_rest = an/(an + bn);
m_rest = am/(am + bm);
h_rest = ah/(ah + bh);

Xrest = [Vrest;
         n_rest;
         m_rest;
         h_rest];

% Symbolic state variables
syms V n m h

% Rate constants
ifexpr_an = -0.01*(V+50)/(exp(-(V+50)/10)-1);
ifexpr_am = -0.1*(V+35)/(exp(-(V+35)/10)-1);

bn_sym = 0.125*exp(-(V+60)/80);
bm_sym = 4*exp(-(V+60)/18);
ah_sym = 0.07*exp(-(V+60)/20);
bh_sym = 1/(exp(-(V+30)/10)+1);

% Membrane equation
IK_sym = p.GK*n^4*(V-p.EK);
INa_sym = p.GNa*m^3*h*(V-p.ENa);
IL_sym = p.GL*(V-p.EL);

dV_sym = (p.Iext - IK_sym - INa_sym - IL_sym)/p.C;

dn_sym = ifexpr_an*(1-n) - bn_sym*n;
dm_sym = ifexpr_am*(1-m) - bm_sym*m;
dh_sym = ah_sym*(1-h) - bh_sym*h;

% Jacobian
JHH = jacobian([dV_sym;dn_sym;dm_sym;dh_sym], ...
               [V n m h]);

% Evaluate at resting state
Jrest = double(subs(JHH, ...
    {V,n,m,h}, ...
    {Vrest,n_rest,m_rest,h_rest}));

fprintf('\nQ14 - HH Jacobian at rest:\n');
disp(Jrest);

% Eigenvalues
lambda_rest = eig(Jrest);

fprintf('Eigenvalues at rest:\n');
disp(lambda_rest);

% --------------------------------------------------------
% Threshold investigation
% --------------------------------------------------------

p.Iext = 0;

dV_values = 0:1:30;

maxV_HH = zeros(size(dV_values));

for i = 1:length(dV_values)

    dV0 = dV_values(i);

    X0 = [Vrest + dV0;
          n_rest;
          m_rest;
          h_rest];

    [t,X] = ode15s(@(t,X)hodgkinHuxley(t,X,p), ...
                   [0 300],X0);

    maxV_HH(i) = max(X(:,1));

end

figure;

plot(dV_values,maxV_HH,'o-','LineWidth',1.5);

xlabel('Initial depolarization \DeltaV_0 (mV)');
ylabel('Maximum V (mV)');
title('HH Threshold: Maximum Voltage vs Initial Depolarization');

grid on;

low = Vrest + 6;
high = Vrest + 7;

for k = 1:30

    mid = (low + high)/2;

    X0 = [mid;
          n_rest;
          m_rest;
          h_rest];

    [~,X] = ode15s(@(t,X)hodgkinHuxley(t,X,p), ...
                   [0 300],X0);

    maxV_mid = max(X(:,1));

    if maxV_mid > 0
        high = mid;
    else
        low = mid;
    end

end

V_threshold_HH = (low + high)/2;

deltaV_threshold_HH = V_threshold_HH - Vrest;

fprintf('\nHH threshold:\n');
fprintf('V_threshold = %.10f mV\n',V_threshold_HH);
fprintf('Delta V_threshold = %.10f mV\n', ...
        deltaV_threshold_HH);

V0_sub = V_threshold_HH - 0.01;
V0_sup = V_threshold_HH + 0.01;

X0_sub = [V0_sub;
          n_rest;
          m_rest;
          h_rest];

X0_sup = [V0_sup;
          n_rest;
          m_rest;
          h_rest];

[t_sub,X_sub] = ode15s( ...
    @(t,X)hodgkinHuxley(t,X,p), ...
    [0 30],X0_sub);

[t_sup,X_sup] = ode15s( ...
    @(t,X)hodgkinHuxley(t,X,p), ...
    [0 30],X0_sup);

figure;

plot(t_sub,X_sub(:,1),'LineWidth',1.5);
hold on;

plot(t_sup,X_sup(:,1),'LineWidth',1.5);

xlabel('Time (ms)');
ylabel('V (mV)');
title('HH Threshold: Subthreshold vs Suprathreshold');

legend('Subthreshold','Suprathreshold');

grid on;

% --------------------------------------------------------
% Q15 - HH steady current injection
% --------------------------------------------------------
currents = 8:0.1:12;

eqV = zeros(size(currents));
eqN = zeros(size(currents));
eqM = zeros(size(currents));
eqH = zeros(size(currents));

% Start near the Q13 resting voltage
Vguess = -55;

for k = 1:length(currents)

    I = currents(k);

    % Equilibrium equation
    F = @(V) ...
        p.GK * (hh_inf(V,1))^4 * (V-p.EK) + ...
        p.GNa * (hh_inf(V,2))^3 * hh_inf(V,3) * (V-p.ENa) + ...
        p.GL * (V-p.EL) - I;

    % Find equilibrium
    Vstar = fzero(F,Vguess);

    % Store equilibrium
    gates = [hh_inf(Vstar,1), ...
             hh_inf(Vstar,2), ...
             hh_inf(Vstar,3)];

    eqV(k) = Vstar;
    eqN(k) = gates(1);
    eqM(k) = gates(2);
    eqH(k) = gates(3);

    % Use this equilibrium as the next guess
    Vguess = Vstar;

end

% --------------------------------------------------------
% Stability of equilibrium for each current
% --------------------------------------------------------

maxRealEig = zeros(size(currents));

for k = 1:length(currents)

    Vstar = eqV(k);
    nstar = eqN(k);
    mstar = eqM(k);
    hstar = eqH(k);

    Jstar = double(subs(JHH, ...
        {V,n,m,h}, ...
        {Vstar,nstar,mstar,hstar}));

    lambda = eig(Jstar);

    maxRealEig(k) = max(real(lambda));

end

figure;

plot(currents,maxRealEig,'o-','LineWidth',1.5);

yline(0,'--');

xlabel('I_{ext} (\muA/cm^2)');
ylabel('Maximum Re(\lambda)');
title('HH Equilibrium Stability');

grid on;


figure;

plot(currents,eqV,'o-','LineWidth',1.5);

xlabel('I_{ext} (\muA/cm^2)');
ylabel('Equilibrium V^* (mV)');
title('HH Equilibrium Voltage vs Applied Current');

grid on;

% Starting point = equilibrium for Iext = 0
X0_zero = [Vrest;
           n_rest;
           m_rest;
           h_rest];

% Apply steady current of 9
p.Iext = 9;

[t_a,X_a] = ode15s( ...
    @(t,X)hodgkinHuxley(t,X,p), ...
    [0 500],X0_zero,options);

Xeq9 = [eqV(currents==9);
        eqN(currents==9);
        eqM(currents==9);
        eqH(currents==9)];

[t_b,X_b] = ode15s( ...
    @(t,X)hodgkinHuxley(t,X,p), ...
    [0 500],Xeq9,options);

figure;

plot(X_a(:,1),X_a(:,2),'LineWidth',1.5);
hold on;

plot(X_b(:,1),X_b(:,2),'LineWidth',1.5);
plot(Xeq9(1), Xeq9(2), 'ko', ...
     'MarkerFaceColor','k', 'MarkerSize',8);

xlabel('V (mV)');
ylabel('n');
title('HH Dynamics at I_{ext}=9 \muA/cm^2');

legend('Start at I=0 equilibrium', ...
       'Start at I=9 equilibrium' ,'');

grid on;


%q 16
V0 = -60;

[an,bn,~,~,~,~] = hh_rates(V0);

n0 = an/(an+bn);

X0 = [V0;n0];

p.Iext = 9;

[t,X] = ode15s(@(t,X)HH_reduced(t,X,p), ...
               [0 100],X0,options);

% --------------------------------------------------------
% Full HH vs Reduced HH action potential
% --------------------------------------------------------

p.Iext = 10;

% Same resting initial condition
X0_full = [Vrest; n_rest; m_rest; h_rest];
X0_red  = [Vrest; n_rest];

% Simulate full HH
[t_full,X_full] = ode15s( ...
    @(t,X) hodgkinHuxley(t,X,p), ...
    [0 100],X0_full,options);

% Simulate reduced HH
[t_red,X_red] = ode15s( ...
    @(t,X) HH_reduced(t,X,p), ...
    [0 100],X0_red,options);

% Overlay voltage
figure;
plot(t_full,X_full(:,1),'LineWidth',1.5);
hold on;
plot(t_red,X_red(:,1),'--','LineWidth',1.5);

xlabel('Time (ms)');
ylabel('V (mV)');
title('Full HH vs Reduced HH');
legend('Full HH','Reduced HH');
grid on;

figure;

plot(X(:,1),X(:,2),'LineWidth',1.5);

xlabel('V (mV)');
ylabel('n');
title('Reduced HH Phase Plane');

grid on;

%Q-17

X0 = [Vrest;
      n_rest;
      m_rest;
      h_rest];
% --------------------------------------------------------
% Q17 - Anode-break excitation
% --------------------------------------------------------

p.Iext = 0;   % current is handled inside pulse function
p.EL = -49.40107903;

options = odeset('AbsTol',1e-6, ...
                 'RelTol',1e-3);

[t17,X17] = ode15s( ...
    @(t,X)hodgkinHuxleyPulse(t,X,p), ...
    [0 50],X0,options);

figure;

plot(t17,X17(:,1),'LineWidth',1.5);
hold on;

xline(20,'--');

xlabel('Time (ms)');
ylabel('V (mV)');
title('HH Anode-Break Excitation');

legend('V(t)','End of hyperpolarizing pulse');

grid on;

figure;

plot(t17,X17(:,2),'LineWidth',1.3);
hold on;
plot(t17,X17(:,3),'LineWidth',1.3);
plot(t17,X17(:,4),'LineWidth',1.3);

xline(20,'--');

xlabel('Time (ms)');
ylabel('Gating variables');

title('Gating Variables During Anode-Break Excitation');

legend('n','m','h','Pulse end');

grid on;

% --------------------------------------------------------
% Q18 - Reduced HH fast subsystem: V and m
% n and h are fixed
% --------------------------------------------------------

fprintf('\n================================================\n');
fprintf('Q18 - Reduced HH Fast Subsystem (V,m)\n');
fprintf('================================================\n');


%% -------------------------------------------------------
% 1. Get n and h at the END of the Q17 pulse
% --------------------------------------------------------

% Interpolate exactly at t = 20 ms
V_end = interp1(t17,X17(:,1),20);
n_end = interp1(t17,X17(:,2),20);
m_end = interp1(t17,X17(:,3),20);
h_end = interp1(t17,X17(:,4),20);

fprintf('\nValues at end of hyperpolarizing pulse:\n');
fprintf('V_end = %.8f mV\n',V_end);
fprintf('n_end = %.8f\n',n_end);
fprintf('m_end = %.8f\n',m_end);
fprintf('h_end = %.8f\n',h_end);


%% -------------------------------------------------------
% 2. Two cases
% -------------------------------------------------------

% Case 1: n and h fixed at resting values
n1 = n_rest;
h1 = h_rest;

% Case 2: n and h fixed at end-of-pulse values
n2 = n_end;
h2 = h_end;


%% -------------------------------------------------------
% 3. m-nullcline
% -------------------------------------------------------

Vnull = linspace(-80,60,3000);

m_mnull = hh_inf(Vnull,2);


%% -------------------------------------------------------
% 4. V-nullclines
% -------------------------------------------------------

m_vnull_1 = nan(size(Vnull));
m_vnull_2 = nan(size(Vnull));

for k = 1:length(Vnull)

    V = Vnull(k);

    % -------------------------
    % Case 1
    % -------------------------

    numerator1 = ...
        -p.GK*n1^4*(V-p.EK) ...
        -p.GL*(V-p.EL);

    denominator1 = ...
        p.GNa*h1*(V-p.ENa);

    if abs(denominator1) > 1e-12

        ratio1 = numerator1/denominator1;

        % m^3 = ratio
        if ratio1 >= 0
            mvalue = ratio1^(1/3);

            % Only physical gating values
            if mvalue >= 0 && mvalue <= 1
                m_vnull_1(k) = mvalue;
            end
        end
    end


    % -------------------------
    % Case 2
    % -------------------------

    numerator2 = ...
        -p.GK*n2^4*(V-p.EK) ...
        -p.GL*(V-p.EL);

    denominator2 = ...
        p.GNa*h2*(V-p.ENa);

    if abs(denominator2) > 1e-12

        ratio2 = numerator2/denominator2;

        if ratio2 >= 0
            mvalue = ratio2^(1/3);

            if mvalue >= 0 && mvalue <= 1
                m_vnull_2(k) = mvalue;
            end
        end
    end

end


%% -------------------------------------------------------
% 5. Find ALL equilibria systematically
%
% At equilibrium:
%
%       m = m_inf(V)
%
% and
%
%       dV/dt = 0
%
% So we solve one scalar equation in V.
% -------------------------------------------------------

% Current balance for Case 1
F1 = @(V) ...
    -p.GK*n1^4*(V-p.EK) ...
    -p.GNa*(hh_inf(V,2)^3)*h1*(V-p.ENa) ...
    -p.GL*(V-p.EL);

% Current balance for Case 2
F2 = @(V) ...
    -p.GK*n2^4*(V-p.EK) ...
    -p.GNa*(hh_inf(V,2)^3)*h2*(V-p.ENa) ...
    -p.GL*(V-p.EL);


% Scan voltage range
Vscan = linspace(-80,60,10000);

roots1 = [];
roots2 = [];

Fscan1 = arrayfun(F1,Vscan);
Fscan2 = arrayfun(F2,Vscan);


% -------------------------
% Case 1 roots
% -------------------------

for k = 1:length(Vscan)-1

    fa = Fscan1(k);
    fb = Fscan1(k+1);

    if isfinite(fa) && isfinite(fb) && fa*fb < 0

        root = fzero(F1,[Vscan(k),Vscan(k+1)]);

        if isempty(roots1) || abs(root-roots1(end)) > 1e-5
            roots1(end+1) = root;
        end

    end
end


% -------------------------
% Case 2 roots
% -------------------------

for k = 1:length(Vscan)-1

    fa = Fscan2(k);
    fb = Fscan2(k+1);

    if isfinite(fa) && isfinite(fb) && fa*fb < 0

        root = fzero(F2,[Vscan(k),Vscan(k+1)]);

        if isempty(roots2) || abs(root-roots2(end)) > 1e-5
            roots2(end+1) = root;
        end

    end
end


% Corresponding m values
eq1 = [roots1(:), hh_inf(roots1(:),2)];

eq2 = [roots2(:), hh_inf(roots2(:),2)];


%% -------------------------------------------------------
% 6. Print equilibria
% -------------------------------------------------------

fprintf('\nCase 1: n,h fixed at resting values\n');

for k = 1:size(eq1,1)
    fprintf('Equilibrium %d: V = %.8f mV, m = %.8f\n', ...
        k,eq1(k,1),eq1(k,2));
end


fprintf('\nCase 2: n,h fixed at end-of-pulse values\n');

for k = 1:size(eq2,1)
    fprintf('Equilibrium %d: V = %.8f mV, m = %.8f\n', ...
        k,eq2(k,1),eq2(k,2));
end


%% -------------------------------------------------------
% 7. Simulate Case 1
% -------------------------------------------------------

% Start exactly at resting equilibrium
X0_1 = [Vrest;m_rest];

[t_vm1,X_vm1] = ode15s( ...
    @(t,X)HH_Vm_reduced(t,X,p,n1,h1), ...
    [0 50],X0_1,options);


%% -------------------------------------------------------
% 8. Simulate Case 2
% -------------------------------------------------------

% Start from state immediately after Q17 pulse
X0_2 = [V_end;m_end];

[t_vm2,X_vm2] = ode15s( ...
    @(t,X)HH_Vm_reduced(t,X,p,n2,h2), ...
    [0 50],X0_2,options);


%% =======================================================
% 9. CASE 1 PHASE PLANE
% =======================================================

figure;
hold on;

% m-nullcline
h_m = plot(Vnull,m_mnull,'LineWidth',2);

% V-nullcline
h_v = plot(Vnull,m_vnull_1,'LineWidth',2);

% Equilibria
h_eq = plot(eq1(:,1),eq1(:,2),'ko', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

% Trajectory
h_traj = plot(X_vm1(:,1),X_vm1(:,2), ...
    'LineWidth',1.5);

xlabel('V (mV)');
ylabel('m');

title('Q18: V-m Phase Plane, n and h at Rest');

legend([h_m,h_v,h_eq,h_traj], ...
    'm-nullcline', ...
    'V-nullcline', ...
    'Equilibria', ...
    'Trajectory', ...
    'Location','best');

xlim([-80 60]);
ylim([0 1]);

grid on;
hold off;


%% =======================================================
% 10. CASE 2 PHASE PLANE
% =======================================================

figure;
hold on;

% m-nullcline
h_m = plot(Vnull,m_mnull,'LineWidth',2);

% V-nullcline
h_v = plot(Vnull,m_vnull_2,'LineWidth',2);

% Equilibria, if any
if ~isempty(eq2)
    h_eq = plot(eq2(:,1),eq2(:,2),'ko', ...
        'MarkerFaceColor','k', ...
        'MarkerSize',8);
end

% Trajectory
h_traj = plot(X_vm2(:,1),X_vm2(:,2), ...
    'LineWidth',1.5);

% Initial point
h_ic = plot(X0_2(1),X0_2(2),'ks', ...
    'MarkerFaceColor','k', ...
    'MarkerSize',8);

xlabel('V (mV)');
ylabel('m');

title('Q18: V-m Phase Plane, n and h at End of Pulse');

if isempty(eq2)

    legend([h_m,h_v,h_traj,h_ic], ...
        'm-nullcline', ...
        'V-nullcline', ...
        'Trajectory', ...
        'Initial condition', ...
        'Location','best');

else

    legend([h_m,h_v,h_eq,h_traj,h_ic], ...
        'm-nullcline', ...
        'V-nullcline', ...
        'Equilibria', ...
        'Trajectory', ...
        'Initial condition', ...
        'Location','best');

end

xlim([-80 60]);
ylim([0 1]);

grid on;
hold off;

% ========================================================
% FUNCTIONS
% ========================================================


function dx = HH_Vm_equations(X,p,n_fixed,h_fixed)

    V = X(1);
    m = X(2);

    [~,~,alpha_m,beta_m,~,~] = hh_rates(V);

    IK = p.GK*n_fixed^4*(V-p.EK);

    INa = p.GNa*m^3*h_fixed*(V-p.ENa);

    IL = p.GL*(V-p.EL);

    dV = (-IK - INa - IL)/p.C;

    dm = alpha_m*(1-m) - beta_m*m;

    dx = [dV;dm];

end

function dXdt = HH_Vm_reduced(t,X,p,n_fixed,h_fixed)

    V = X(1);
    m = X(2);

    % Rate constants
    [~,~,alpha_m,beta_m,~,~] = hh_rates(V);

    % Ionic currents
    IK = p.GK*n_fixed^4*(V-p.EK);

    INa = p.GNa*m^3*h_fixed*(V-p.ENa);

    IL = p.GL*(V-p.EL);

    % Current after anode-break pulse
    Iext = 0;

    % Differential equations
    dV = (Iext - IK - INa - IL)/p.C;

    dm = alpha_m*(1-m) - beta_m*m;

    dXdt = [dV;dm];

end

function dXdt = hodgkinHuxleyPulse(t,X,p)

    V = X(1);
    n = X(2);
    m = X(3);
    h = X(4);

    % Rate constants
    [alpha_n,beta_n,alpha_m,beta_m,alpha_h,beta_h] = ...
        hh_rates(V);

    % Hyperpolarizing pulse
    if t < 20
        Iext = -3;
    else
        Iext = 0;
    end

    % Ionic currents
    IK  = p.GK*n^4*(V-p.EK);
    INa = p.GNa*m^3*h*(V-p.ENa);
    IL  = p.GL*(V-p.EL);

    % Voltage equation
    dV = (Iext - IK - INa - IL)/p.C;

    % Gating equations
    dn = p.phi*(alpha_n*(1-n) - beta_n*n);
    dm = p.phi*(alpha_m*(1-m) - beta_m*m);
    dh = p.phi*(alpha_h*(1-h) - beta_h*h);

    dXdt = [dV;dn;dm;dh];

end



function dXdt = HH_reduced(t,X,p)

    V = X(1);
    n = X(2);

    % Fast sodium activation
    [~,~,alpha_m,beta_m,~,~] = hh_rates(V);

    m_inf = alpha_m/(alpha_m + beta_m);

    % Approximate sodium inactivation
    h = 0.8 - n;

    % Ionic currents
    IK = p.GK*n^4*(V-p.EK);

    INa = p.GNa*m_inf^3*h*(V-p.ENa);

    IL = p.GL*(V-p.EL);

    % Voltage equation
    dV = (p.Iext - IK - INa - IL)/p.C;

    % Potassium gating equation
    [alpha_n,beta_n,~,~,~,~] = hh_rates(V);

    dn = alpha_n*(1-n) - beta_n*n;

    dXdt = [dV;dn];

end

function value = hh_inf(V,which)

    [an,bn,am,bm,ah,bh] = hh_rates(V);

    ninf = an./(an+bn);
    minf = am./(am+bm);
    hinf = ah./(ah+bh);

    if which == 1
        value = ninf;
    elseif which == 2
        value = minf;
    elseif which == 3
        value = hinf;
    end

end

function [alpha_n,beta_n,alpha_m,beta_m,alpha_h,beta_h] = hh_rates(V)

    % alpha_n
    if isscalar(V)
        if abs(V + 50) < 1e-8
            alpha_n = 0.1;
        else
            alpha_n = -0.01*(V + 50) / ...
                      (exp(-(V + 50)/10) - 1);
        end
    else
        alpha_n = -0.01*(V + 50) ./ ...
                  (exp(-(V + 50)./10) - 1);
    end

    % beta_n
    beta_n = 0.125*exp(-(V + 60)./80);

    % alpha_m
    if isscalar(V)
        if abs(V + 35) < 1e-8
            alpha_m = 1.0;
        else
            alpha_m = -0.1*(V + 35) / ...
                      (exp(-(V + 35)/10) - 1);
        end
    else
        alpha_m = -0.1*(V + 35) ./ ...
                  (exp(-(V + 35)./10) - 1);
    end

    % beta_m
    beta_m = 4*exp(-(V + 60)./18);

    % alpha_h
    alpha_h = 0.07*exp(-(V + 60)./20);

    % beta_h
    beta_h = 1 ./ ...
        (exp(-(V + 30)./10) + 1);

end


function dXdt = hodgkinHuxley(t,X,p)

    V = X(1);
    n = X(2);
    m = X(3);
    h = X(4);

    % Rate constants
    [alpha_n,beta_n,alpha_m,beta_m,alpha_h,beta_h] = ...
        hh_rates(V);

    % Ionic currents
    IK  = p.GK*n^4*(V-p.EK);
    INa = p.GNa*m^3*h*(V-p.ENa);
    IL  = p.GL*(V-p.EL);

    % Membrane equation
    dV = (p.Iext - IK - INa - IL)/p.C;

    % Gating equations
    dn = p.phi*(alpha_n*(1-n) - beta_n*n);
    dm = p.phi*(alpha_m*(1-m) - beta_m*m);
    dh = p.phi*(alpha_h*(1-h) - beta_h*h);

    dXdt = [dV;dn;dm;dh];

end


%MorrisLEcar fun:
function dxdt = morrisLecar(t, x, p)

    V = x(1);
    w = x(2);

    minf = 0.5*(1 + tanh((V-p.V1)/p.V2));
    winf = 0.5*(1 + tanh((V-p.V3)/p.V4));

    tauw = 1/cosh((V-p.V5)/p.V6);

    dV = (p.Iext ...
        - p.gCa*minf*(V-p.VCa) ...
        - p.gK*w*(V-p.VK) ...
        - p.gL*(V-p.VL))/p.C;

    dw = p.phi*(winf-w)/tauw;

    dxdt = [dV; dw];
end

function [value,isterminal,direction] = manifoldEvent(t,x)

    V = x(1);
    w = x(2);

    % Stop if V leaves the displayed region
    % or w becomes physically irrelevant
    value = min([ ...
        V + 80, ...
        50 - V, ...
        w, ...
        0.5 - w]);

    isterminal = 1;
    direction = -1;

end