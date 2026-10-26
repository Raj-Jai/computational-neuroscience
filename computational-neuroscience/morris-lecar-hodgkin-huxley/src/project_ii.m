%% Computational Neuroscience Project 2
% Refined MATLAB implementation for better formating using AI (AI only used for refinement)


clear;
clc;
close all;

%% ============================================================
%  MORRIS-LECAR PARAMETERS AND COMMON SOLVER SETTINGS
% =============================================================

ml.gCa = 4.4;
ml.gK  = 8.0;
ml.gL  = 2.0;

ml.VCa = 120;
ml.VK  = -84;
ml.VL  = -60;

ml.phi = 0.02;

ml.V1 = -1.2;
ml.V2 = 18;
ml.V3 = 2;
ml.V4 = 30;
ml.V5 = 2;
ml.V6 = 30;

ml.C = 20;
ml.Iext = 0;

% General plotting / solver settings
odeOptsML = odeset('RelTol',1e-6,'AbsTol',1e-8);
odeOptsHH = odeset('RelTol',1e-3,'AbsTol',1e-6);

%% ============================================================
%  Q1-Q2 - MLE EQUILIBRIUM, VECTOR FIELD AND NULLCLINES
% =============================================================

mInf = @(V) 0.5*(1 + tanh((V-ml.V1)/ml.V2));
wInf = @(V) 0.5*(1 + tanh((V-ml.V3)/ml.V4));
tauW = @(V) 1 ./ cosh((V-ml.V5)/ml.V6);

F_ml_eq = @(V) ...
    ml.Iext ...
    - ml.gCa*mInf(V)*(V-ml.VCa) ...
    - ml.gK*wInf(V)*(V-ml.VK) ...
    - ml.gL*(V-ml.VL);

% Equilibrium from reduced one-variable equation
V_eq = fzero(F_ml_eq,-60);
w_eq = wInf(V_eq);

fprintf('\n================ Q1-Q2: Morris-Lecar ================\n');
fprintf('Equilibrium from reduced equation:\n');
fprintf('V* = %.8f mV\n',V_eq);
fprintf('w* = %.8f\n',w_eq);

% Equilibrium from simultaneous equations
F_ml_2 = @(x) ml_rhs(x,ml);
x_eq = fsolve(F_ml_2,[V_eq;w_eq]);

fprintf('Equilibrium from simultaneous equations:\n');
fprintf('V* = %.8f mV\n',x_eq(1));
fprintf('w* = %.8f\n',x_eq(2));

% Phase-plane grid
Vgrid = linspace(-80,50,60);
wgrid = linspace(0,0.5,45);
[V,W] = meshgrid(Vgrid,wgrid);

[dV,dW] = ml_vector_field(V,W,ml);

% Plot in (V,100w) coordinates
dVplot = dV;
dWplot = 100*dW;
speed = hypot(dVplot,dWplot);
dVn = dVplot ./ (speed + eps);
dWn = dWplot ./ (speed + eps);

figure('Name','MLE Phase Plane');
quiver(V,100*W,dVn,dWn,0.7,'LineWidth',1);
hold on;

Vnull = linspace(-80,50,1200);
W_w = wInf(Vnull);
W_v = ( ...
    ml.Iext ...
    - ml.gCa*mInf(Vnull).*(Vnull-ml.VCa) ...
    - ml.gL*(Vnull-ml.VL)) ...
    ./ (ml.gK*(Vnull-ml.VK));

hVF = plot(nan,nan,'LineWidth',1);
hW  = plot(Vnull,100*W_w,'LineWidth',2);
hV  = plot(Vnull,100*W_v,'LineWidth',2);
hEQ = plot(V_eq,100*w_eq,'ko','MarkerFaceColor','k','MarkerSize',8);

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Phase Plane');
xlim([-80 50]);
ylim([0 50]);
legend([hVF,hW,hV,hEQ], ...
    'Vector field','w-nullcline','V-nullcline','Equilibrium', ...
    'Location','northeast');
grid on;
hold off;

%% ============================================================
%  Q3 - MLE JACOBIAN AND EIGENVALUES
% =============================================================

J_ml = numericalJacobian(@(x) ml_rhs(x,ml),[V_eq;w_eq],[1e-5;1e-7]);
lambda_ml = eig(J_ml);

fprintf('\n================ Q3: MLE Stability ================\n');
fprintf('Jacobian at equilibrium:\n');
disp(J_ml);
fprintf('Eigenvalues:\n');
disp(lambda_ml);

%% ============================================================
%  Q5 - EFFECT OF phi ON ACTION POTENTIAL
% =============================================================

phis = [0.01 0.02 0.04];

figure('Name','MLE Action Potentials');
hold on;

for k = 1:numel(phis)
    p5 = ml;
    p5.phi = phis(k);

    x0 = [V_eq + 50; w_eq];

    [t5,x5] = ode45(@(t,x) morrisLecar(t,x,p5), ...
        [0 300],x0,odeOptsML);

    plot(t5,x5(:,1),'LineWidth',1.5);
end

xlabel('Time (ms)');
ylabel('V (mV)');
title('Morris-Lecar Action Potentials');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;
hold off;

% Phase plane comparison
figure('Name','MLE Phase Plane - phi');
hold on;

for k = 1:numel(phis)
    p5 = ml;
    p5.phi = phis(k);

    x0 = [V_eq + 50; w_eq];

    [~,x5] = ode45(@(t,x) morrisLecar(t,x,p5), ...
        [0 300],x0,odeOptsML);

    plot(x5(:,1),100*x5(:,2),'LineWidth',1.5);
end

plot(Vnull,100*W_w,'LineWidth',2);
plot(Vnull,100*W_v,'LineWidth',2);
plot(V_eq,100*w_eq,'ko','MarkerFaceColor','k','MarkerSize',8);
plot(V_eq+50,100*w_eq,'ks','MarkerFaceColor','k','MarkerSize',7);

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Phase Plane');
xlim([-80 50]);
ylim([0 50]);
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04', ...
       'w-nullcline','V-nullcline','Equilibrium','Initial condition');
grid on;
hold off;

%% ============================================================
%  Q6 - MLE THRESHOLD DEPOLARIZATION
% =============================================================

p6 = ml;
p6.phi = 0.02;
p6.Iext = 0;

dV0_scan = 0:1:80;
maxV_scan = nan(size(dV0_scan));

for k = 1:numel(dV0_scan)
    x0 = [V_eq + dV0_scan(k);w_eq];

    [~,x6] = ode45(@(t,x) morrisLecar(t,x,p6), ...
        [0 300],x0,odeOptsML);

    maxV_scan(k) = max(x6(:,1));
end

figure('Name','MLE Threshold Scan');
plot(dV0_scan,maxV_scan,'o-','LineWidth',1.5);
xlabel('Initial depolarization \DeltaV_0 (mV)');
ylabel('Maximum V (mV)');
title('MLE Threshold: Maximum Voltage vs Initial Depolarization');
grid on;

% Automatically find the first subthreshold/suprathreshold bracket
isAP = maxV_scan > 0;
idxAP = find(isAP,1,'first');

if isempty(idxAP) || idxAP == 1
    error('Q6: Could not automatically bracket the threshold.');
end

low = V_eq + dV0_scan(idxAP-1);
high = V_eq + dV0_scan(idxAP);

for k = 1:35
    mid = (low + high)/2;

    x0 = [mid;w_eq];

    [~,x6] = ode45(@(t,x) morrisLecar(t,x,p6), ...
        [0 300],x0,odeOptsML);

    if max(x6(:,1)) > 0
        high = mid;
    else
        low = mid;
    end
end

V_threshold = (low + high)/2;
deltaV_threshold = V_threshold - V_eq;

fprintf('\n================ Q6: MLE Threshold ================\n');
fprintf('V_threshold = %.10f mV\n',V_threshold);
fprintf('Delta V_threshold = %.10f mV\n',deltaV_threshold);

V0_sub = V_threshold - 0.1;
V0_sup = V_threshold + 0.1;

[tSub,xSub] = ode45(@(t,x) morrisLecar(t,x,p6), ...
    [0 300],[V0_sub;w_eq],odeOptsML);
[tSup,xSup] = ode45(@(t,x) morrisLecar(t,x,p6), ...
    [0 300],[V0_sup;w_eq],odeOptsML);

figure('Name','MLE Threshold Trajectories');
hold on;
h1 = plot(xSub(:,1),100*xSub(:,2),'LineWidth',1.5);
h2 = plot(xSup(:,1),100*xSup(:,2),'LineWidth',1.5);
h3 = plot(V_eq,100*w_eq,'ko','MarkerFaceColor','k');
h4 = plot(V0_sub,100*w_eq,'ks','MarkerFaceColor','k');
plot(V0_sup,100*w_eq,'ks','MarkerFaceColor','k');
h5 = plot(Vnull,100*W_w,'LineWidth',2,'Color','k');
h6 = plot(Vnull,100*W_v,'LineWidth',2,'Color','k');

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Threshold Trajectories');
xlim([-80 50]);
ylim([0 50]);
legend([h1,h2,h3,h4,h5,h6], ...
    'Subthreshold','Suprathreshold','Equilibrium', ...
    'Initial conditions','w-nullcline','V-nullcline');
grid on;
hold off;

%% ============================================================
%  Q7 - MLE AT Iext = 86
% =============================================================

p7 = ml;
p7.Iext = 86;
p7.phi = 0.02;

x_eq86 = fsolve(@(x) ml_rhs(x,p7),[-28;0.12]);

fprintf('\n================ Q7: MLE at Iext = 86 ================\n');
fprintf('V* = %.8f mV\n',x_eq86(1));
fprintf('w* = %.8f\n',x_eq86(2));

x01 = [V_eq;w_eq];
x02 = x_eq86;
x03 = [-27.9;0.17];

[t1,x1] = ode45(@(t,x) morrisLecar(t,x,p7),[0 500],x01,odeOptsML);
[t2,x2] = ode45(@(t,x) morrisLecar(t,x,p7),[0 500],x02,odeOptsML);
[t3,x3] = ode45(@(t,x) morrisLecar(t,x,p7),[0 500],x03,odeOptsML);

figure('Name','MLE I=86 Phase Plane');
hold on;
h1 = plot(x1(:,1),100*x1(:,2),'LineWidth',1.5);
h2 = plot(x2(:,1),100*x2(:,2),'LineWidth',1.5);
h3 = plot(x3(:,1),100*x3(:,2),'LineWidth',1.5);
h4 = plot(x_eq86(1),100*x_eq86(2),'ko','MarkerFaceColor','k','MarkerSize',8);
h5 = plot(x01(1),100*x01(2),'kd','MarkerFaceColor','k','MarkerSize',8);
h6 = plot(x03(1),100*x03(2),'ks','MarkerFaceColor','k','MarkerSize',8);

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar at I_{ext}=86 \muA/cm^2');
legend([h1,h2,h3,h4,h5,h6], ...
    'Trial 1: I=0 equilibrium', ...
    'Trial 2: I=86 equilibrium', ...
    'Trial 3: (-27.9,0.17)', ...
    'I=86 equilibrium','Trial 1 IC','Trial 3 IC');
grid on;
hold off;

J86 = numericalJacobian(@(x) ml_rhs(x,p7),x_eq86,[1e-5;1e-7]);
lambda86 = eig(J86);

fprintf('\nJacobian at I=86 equilibrium:\n');
disp(J86);
fprintf('Eigenvalues:\n');
disp(lambda86);

%% ============================================================
%  Q9 - MLE EQUILIBRIUM / STABILITY FOR I = 80,86,90
% =============================================================

currentsQ9 = [80 86 90];
eqQ9 = zeros(numel(currentsQ9),2);
eigQ9 = zeros(numel(currentsQ9),2);

for k = 1:numel(currentsQ9)
    p9 = ml;
    p9.Iext = currentsQ9(k);
    p9.phi = 0.02;

    eqQ9(k,:) = fsolve(@(x) ml_rhs(x,p9),[-28;0.12]).';

    J9 = numericalJacobian(@(x) ml_rhs(x,p9),eqQ9(k,:).',[1e-5;1e-7]);
    eigQ9(k,:) = eig(J9).';

    fprintf('\nQ9: Iext = %.1f\n',currentsQ9(k));
    fprintf('V* = %.8f mV\n',eqQ9(k,1));
    fprintf('w* = %.8f\n',eqQ9(k,2));
    fprintf('Eigenvalues:\n');
    disp(eigQ9(k,:).');
end

%% ============================================================
%  Q9 - MLE FIRING RATE, I = 80 TO 100
% =============================================================

currentsFR9 = 80:1:100;
firingRate9 = zeros(size(currentsFR9));

for k = 1:numel(currentsFR9)
    p9 = ml;
    p9.Iext = currentsFR9(k);
    p9.phi = 0.02;

    xeq = fsolve(@(x) ml_rhs(x,p9),[-28;0.12]);
    J9 = numericalJacobian(@(x) ml_rhs(x,p9),xeq,[1e-5;1e-7]);

    x0 = xeq;
    if max(real(eig(J9))) > 0
        x0 = x0 + [1e-5;0];
    end

    [t9,x9] = ode45(@(t,x) morrisLecar(t,x,p9), ...
        [0 3000],x0,odeOptsML);

    post = t9 > 1000;
    tSS = t9(post);
    VSS = x9(post,1);

    idxCross = find(VSS(1:end-1)<0 & VSS(2:end)>=0);
    spikeTimes = tSS(idxCross);

    if numel(spikeTimes) >= 2
        firingRate9(k) = 1000/mean(diff(spikeTimes));
    end
end

figure('Name','MLE Firing Rate');
plot(currentsFR9,firingRate9,'o-','LineWidth',1.5);
xlabel('Applied current I_{ext} (\muA/cm^2)');
ylabel('Firing rate (Hz)');
title('Morris-Lecar Firing Rate vs Applied Current');
grid on;

%% ============================================================
%  Q10 - MLE: THREE EQUILIBRIA AND MANIFOLDS, Iext = 30
% =============================================================

ml10 = ml;

ml10.gCa = 4.0;
ml10.gK  = 8.0;
ml10.gL  = 2.0;

ml10.VCa = 120;
ml10.VK  = -84;
ml10.VL  = -60;

ml10.phi = 0.0667;

ml10.V1 = -1.2;
ml10.V2 = 18;
ml10.V3 = 12;
ml10.V4 = 17.4;
ml10.V5 = 12;
ml10.V6 = 17.4;

ml10.C = 20;
ml10.Iext = 30;

mInf10 = @(V) 0.5*(1+tanh((V-ml10.V1)/ml10.V2));
wInf10 = @(V) 0.5*(1+tanh((V-ml10.V3)/ml10.V4));

Ieq10 = @(V) ...
    ml10.gCa*mInf10(V).*(V-ml10.VCa) + ...
    ml10.gK*wInf10(V).*(V-ml10.VK) + ...
    ml10.gL*(V-ml10.VL);

Vscan10 = linspace(-80,50,10000);
Fscan10 = arrayfun(@(V) Ieq10(V)-ml10.Iext,Vscan10);
roots10 = findAllRoots(@(V) Ieq10(V)-ml10.Iext,Vscan10,Fscan10);

eq10 = [roots10(:),wInf10(roots10(:))];

fprintf('\n================ Q10: MLE I=30 ================\n');
for k = 1:size(eq10,1)
    fprintf('Equilibrium %d: V = %.8f mV, w = %.8f\n', ...
        k,eq10(k,1),eq10(k,2));
end

% Stability classification
fprintf('\nQ10 eigenvalues:\n');
for k = 1:size(eq10,1)
    J10k = numericalJacobian(@(x) ml_rhs(x,ml10),eq10(k,:).',[1e-5;1e-7]);
    lam10 = eig(J10k);
    fprintf('\nEquilibrium %d:\n',k);
    disp(lam10);
end

% Q10 phase plane
Vnull10 = linspace(-80,50,1600);
Ww10 = wInf10(Vnull10);
Wv10 = ( ...
    ml10.Iext ...
    - ml10.gCa*mInf10(Vnull10).*(Vnull10-ml10.VCa) ...
    - ml10.gL*(Vnull10-ml10.VL)) ...
    ./ (ml10.gK*(Vnull10-ml10.VK));

%{
figure('Name','Q10 MLE Phase Plane');
hold on;
h1 = plot(Vnull10,100*Ww10,'LineWidth',2);
h2 = plot(Vnull10,100*Wv10,'LineWidth',2);
h3 = plot(eq10(:,1),100*eq10(:,2),'ko','MarkerFaceColor','k','MarkerSize',8);
xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Phase Plane: I_{ext}=30 \muA/cm^2');
xlim([-80 50]);
ylim([0 50]);
legend([h1,h2,h3],'w-nullcline','V-nullcline','Equilibria','Location','best');
grid on;
hold off; 
%}

% Identify saddle
midIdx = 2;
Vs = eq10(midIdx,1);
ws = eq10(midIdx,2);

Js = numericalJacobian(@(x) ml_rhs(x,ml10),[Vs;ws],[1e-5;1e-7]);
[Vecs,Vals] = eig(Js);
lamS = diag(Vals);

stableIdx = find(real(lamS)<0,1);
unstableIdx = find(real(lamS)>0,1);

vs = real(Vecs(:,stableIdx));   vs = vs/norm(vs);
vu = real(Vecs(:,unstableIdx)); vu = vu/norm(vu);

fprintf('\nSaddle equilibrium:\n');
fprintf('V = %.8f mV\n',Vs);
fprintf('w = %.8f\n',ws);
fprintf('Stable eigenvalue = %.8f%+.8fi\n',real(lamS(stableIdx)),imag(lamS(stableIdx)));
fprintf('Unstable eigenvalue = %.8f%+.8fi\n',real(lamS(unstableIdx)),imag(lamS(unstableIdx)));

epsilonMan = 1e-6;
saddle = [Vs;ws];

xs1 = saddle + epsilonMan*vs;
xs2 = saddle - epsilonMan*vs;
xu1 = saddle + epsilonMan*vu;
xu2 = saddle - epsilonMan*vu;

manOpts = odeset('Events',@(t,x) manifoldEvent(t,x,[-80 50],[0 0.5]), ...
                 'RelTol',1e-7,'AbsTol',1e-9);

[tS1,xS1] = ode45(@(t,x) morrisLecar(t,x,ml10),[0 -150],xs1,manOpts);
[tS2,xS2] = ode45(@(t,x) morrisLecar(t,x,ml10),[0 -150],xs2,manOpts);

[tU1,xU1] = ode45(@(t,x) morrisLecar(t,x,ml10),[0 150],xu1,manOpts);
[tU2,xU2] = ode45(@(t,x) morrisLecar(t,x,ml10),[0 150],xu2,manOpts);

figure('Name','Q10 MLE Manifolds');
hold on;
plot(Vnull10,100*Ww10,'LineWidth',2);
plot(Vnull10,100*Wv10,'LineWidth',2);
plot(eq10(:,1),100*eq10(:,2),'ko','MarkerFaceColor','k','MarkerSize',8);

hS1 = plot(xS1(:,1),100*xS1(:,2),'b--','LineWidth',2);
plot(xS2(:,1),100*xS2(:,2),'b--','LineWidth',2);

hU1 = plot(xU1(:,1),100*xU1(:,2),'g:','LineWidth',2);
plot(xU2(:,1),100*xU2(:,2),'g:','LineWidth',2);

hSaddle = plot(Vs,100*ws,'ks','MarkerFaceColor','k','MarkerSize',10);

xlabel('V (mV)');
ylabel('100w');
title('Morris-Lecar Phase Plane: I_{ext}=30 \muA/cm^2');
xlim([-80 50]);
ylim([0 50]);
legend({'w-nullcline','V-nullcline','Equilibria', ...
        'Stable manifold','','Unstable manifold','','Saddle point'}, ...
        'Location','best');
grid on;
hold off;

%% ============================================================
%  Q11 - MLE EQUILIBRIUM CONTINUATION, I = 30 TO 50
% =============================================================

currents11 = 30:0.05:50;
eqV11 = cell(size(currents11));
eqW11 = cell(size(currents11));
type11 = cell(size(currents11));

Vscan11 = linspace(-80,50,6000);

for k = 1:numel(currents11)
    I = currents11(k);
    F11 = @(V) ...
        ml10.gCa*mInf10(V).*(V-ml10.VCa) + ...
        ml10.gK*wInf10(V).*(V-ml10.VK) + ...
        ml10.gL*(V-ml10.VL) - I;

    Fscan = arrayfun(F11,Vscan11);
    roots = findAllRoots(F11,Vscan11,Fscan);

    eqV11{k} = roots;
    eqW11{k} = wInf10(roots);

    types = cell(size(roots));

    for j = 1:numel(roots)
        xx = [roots(j);wInf10(roots(j))];
        J11k = numericalJacobian(@(x) ml_rhs(x,ml10, I),xx,[1e-5;1e-7]);
        lam = eig(J11k);
        types{j} = classify2DEquilibrium(lam);
    end

    type11{k} = types;
end

figure('Name','Q11 MLE Equilibrium Branches');
hold on;

hStable = plot(nan,nan,'bo','MarkerSize',5);
hSaddle = plot(nan,nan,'ks','MarkerSize',5);
hUnstable = plot(nan,nan,'rx','MarkerSize',6);

for k = 1:numel(currents11)
    for j = 1:numel(eqV11{k})
        switch type11{k}{j}
            case {'Stable node','Stable focus'}
                plot(currents11(k),eqV11{k}(j),'bo','MarkerSize',5);
            case 'Saddle'
                plot(currents11(k),eqV11{k}(j),'ks','MarkerSize',5);
            case {'Unstable node','Unstable focus'}
                plot(currents11(k),eqV11{k}(j),'rx','MarkerSize',6);
        end
    end
end

xlabel('I_{ext} (\muA/cm^2)');
ylabel('Equilibrium V^* (mV)');
title('Morris-Lecar Equilibrium Branches');
legend([hStable,hSaddle,hUnstable], ...
    'Stable','Saddle','Unstable');
grid on;
hold off;

figure('Name','Q11 MLE Equilibria Near 39-40');
hold on;

for k = 1:numel(currents11)
    if currents11(k) < 39 || currents11(k) > 40
        continue;
    end

    for j = 1:numel(eqV11{k})
        switch type11{k}{j}
            case {'Stable node','Stable focus'}
                plot(currents11(k),eqV11{k}(j),'bo','MarkerSize',5);
            case 'Saddle'
                plot(currents11(k),eqV11{k}(j),'ks','MarkerSize',5);
            case {'Unstable node','Unstable focus'}
                plot(currents11(k),eqV11{k}(j),'rx','MarkerSize',6);
        end
    end
end

xlabel('I_{ext} (\muA/cm^2)');
ylabel('Equilibrium V^* (mV)');
title('Equilibrium Behavior Near 39-40 \muA/cm^2');
xlim([39 40]);
grid on;
hold off;

% Firing rate
currents11FR = 30:0.1:45;
firingRate11 = zeros(size(currents11FR));

for k = 1:numel(currents11FR)
    I = currents11FR(k);

    [roots,~] = mlEquilibria(I,ml10,Vscan11);

    if isempty(roots)
        continue;
    end

    % Highest-voltage equilibrium is the branch used in the original code.
    Veq = roots(end);
    weq = wInf10(Veq);

    xeq = [Veq;weq];
    Jk = numericalJacobian(@(x) ml_rhs(x,ml10,I),xeq,[1e-5;1e-7]);
    unstable = max(real(eig(Jk))) > 0;

    x0 = xeq;
    if unstable
        x0 = x0 + [1e-5;0];
    end

    p11 = ml10;
    p11.Iext = I;

    [t11,x11] = ode45(@(t,x) morrisLecar(t,x,p11), ...
        [0 5000],x0,odeOptsML);

    post = t11 > 2500;
    tSS = t11(post);
    VSS = x11(post,1);

    idxCross = find(VSS(1:end-1)<0 & VSS(2:end)>=0);
    spikeTimes = tSS(idxCross);

    if numel(spikeTimes) >= 2
        firingRate11(k) = 1000/mean(diff(spikeTimes));
    end
end

figure('Name','Q11 MLE Firing Rate');
plot(currents11FR,firingRate11,'o-','LineWidth',1.5);
xlabel('I_{ext} (\muA/cm^2)');
ylabel('Firing rate (Hz)');
title('Morris-Lecar Firing Rate');
grid on;

%% ============================================================
%  HODGKIN-HUXLEY PARAMETERS
% =============================================================

hh.GK  = 36;
hh.EK  = -72;

hh.GNa = 120;
hh.ENa = 55;

hh.GL  = 0.3;
hh.EL  = -60;       % Temporary for Q12
hh.C   = 1;

hh.phi = 1;         % T = 6.3 deg C
hh.Iext = 0;

%% ============================================================
%  Q12 - FULL HH MODEL
% =============================================================

V0 = -60;

[an,bn,am,bm,ah,bh] = hh_rates(V0);

n0 = an/(an+bn);
m0 = am/(am+bm);
h0 = ah/(ah+bh);

X0 = [V0;n0;m0;h0];

% Tighter numerical check used for Q12 implementation
q12Opts = odeset('RelTol',1e-6,'AbsTol',1e-9);

[t12,X12] = ode45(@(t,X) hodgkinHuxley(t,X,hh), ...
    [0 25],X0,q12Opts);

fprintf('\n================ Q12: Full HH ================\n');
fprintf('Initial V = %.8f mV\n',X12(1,1));
fprintf('Final V   = %.8f mV\n',X12(end,1));

figure('Name','Q12 HH Membrane Potential');
plot(t12,X12(:,1),'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('V (mV)');
title('Hodgkin-Huxley Membrane Potential');
grid on;

figure('Name','Q12 HH Gating Variables');
plot(t12,X12(:,2),'LineWidth',1.2);
hold on;
plot(t12,X12(:,3),'LineWidth',1.2);
plot(t12,X12(:,4),'LineWidth',1.2);
xlabel('Time (ms)');
ylabel('Gating variables');
title('Hodgkin-Huxley Gating Variables');
legend('n','m','h');
grid on;
hold off;

%% ============================================================
%  Q13 - DETERMINE E_L FOR Vrest = -60 mV
% =============================================================

Vrest = -60;
[an,bn,am,bm,ah,bh] = hh_rates(Vrest);

nRest = an/(an+bn);
mRest = am/(am+bm);
hRest = ah/(ah+bh);

IKrest = hh.GK*nRest^4*(Vrest-hh.EK);
INarest = hh.GNa*mRest^3*hRest*(Vrest-hh.ENa);

hh.EL = Vrest + (IKrest + INarest)/hh.GL;
hh.Iext = 0;

fprintf('\n================ Q13: HH Resting State ================\n');
fprintf('n_rest = %.8f\n',nRest);
fprintf('m_rest = %.8f\n',mRest);
fprintf('h_rest = %.8f\n',hRest);
fprintf('IK_rest  = %.8f uA/cm^2\n',IKrest);
fprintf('INa_rest = %.8f uA/cm^2\n',INarest);
fprintf('EL = %.8f mV\n',hh.EL);

Xrest = [Vrest;nRest;mRest;hRest];

[t13,X13] = ode15s(@(t,X) hodgkinHuxley(t,X,hh), ...
    [0 100],Xrest,odeOptsHH);

figure('Name','Q13 HH Rest Verification');
plot(t13,X13(:,1),'LineWidth',1.5);
yline(Vrest,'--');
xlabel('Time (ms)');
ylabel('V (mV)');
title('HH Resting Potential Verification');
grid on;

%% ============================================================
%  Q14 - HH RESTING STABILITY AND THRESHOLD
% =============================================================

% Full HH Jacobian at rest
Jrest = numericalJacobian(@(x) hh_rhs(x,hh),Xrest,[1e-4;1e-6;1e-6;1e-6]);
lambdaRest = eig(Jrest);

fprintf('\n================ Q14: HH Stability ================\n');
fprintf('Jacobian at rest:\n');
disp(Jrest);
fprintf('Eigenvalues:\n');
disp(lambdaRest);

% Threshold scan
dV14 = 0:0.5:20;
maxV14 = zeros(size(dV14));

for k = 1:numel(dV14)
    X0 = [Vrest+dV14(k);nRest;mRest;hRest];

    [~,X14] = ode15s(@(t,X) hodgkinHuxley(t,X,hh), ...
        [0 300],X0,odeOptsHH);

    maxV14(k) = max(X14(:,1));
end

figure('Name','Q14 HH Threshold Scan');
plot(dV14,maxV14,'o-','LineWidth',1.5);
xlabel('Initial depolarization \DeltaV_0 (mV)');
ylabel('Maximum V (mV)');
title('HH Threshold: Maximum Voltage vs Initial Depolarization');
grid on;

idxAP14 = find(maxV14 > 0,1);

if isempty(idxAP14) || idxAP14 == 1
    error('Q14: Could not automatically bracket the threshold.');
end

low = Vrest + dV14(idxAP14-1);
high = Vrest + dV14(idxAP14);

for k = 1:35
    mid = (low+high)/2;

    X0 = [mid;nRest;mRest;hRest];

    [~,X14] = ode15s(@(t,X) hodgkinHuxley(t,X,hh), ...
        [0 300],X0,odeOptsHH);

    if max(X14(:,1)) > 0
        high = mid;
    else
        low = mid;
    end
end

VthresholdHH = (low+high)/2;
deltaVthresholdHH = VthresholdHH-Vrest;

fprintf('\nHH threshold:\n');
fprintf('V_threshold = %.10f mV\n',VthresholdHH);
fprintf('Delta V_threshold = %.10f mV\n',deltaVthresholdHH);

VsubHH = VthresholdHH-0.01;
VsupHH = VthresholdHH+0.01;

[tSubHH,XsubHH] = ode15s( ...
    @(t,X)hodgkinHuxley(t,X,hh), ...
    [0 30],[VsubHH;nRest;mRest;hRest],odeOptsHH);

[tSupHH,XsupHH] = ode15s( ...
    @(t,X)hodgkinHuxley(t,X,hh), ...
    [0 30],[VsupHH;nRest;mRest;hRest],odeOptsHH);

figure('Name','Q14 HH Threshold Responses');
plot(tSubHH,XsubHH(:,1),'LineWidth',1.5);
hold on;
plot(tSupHH,XsupHH(:,1),'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('V (mV)');
title('HH Threshold: Subthreshold vs Suprathreshold');
legend('Subthreshold','Suprathreshold');
grid on;
hold off;

%% ============================================================
%  Q15 - HH STEADY CURRENT INJECTION, I = 8 TO 12
% =============================================================

currents15 = 8:0.1:12;

eqV15 = zeros(size(currents15));
eqN15 = zeros(size(currents15));
eqM15 = zeros(size(currents15));
eqH15 = zeros(size(currents15));
maxReal15 = zeros(size(currents15));

Vguess = -55;

for k = 1:numel(currents15)
    I = currents15(k);

    F15 = @(V) ...
        hh.GK*hh_inf(V,1)^4*(V-hh.EK) + ...
        hh.GNa*hh_inf(V,2)^3*hh_inf(V,3)*(V-hh.ENa) + ...
        hh.GL*(V-hh.EL) - I;

    Vstar = fzero(F15,Vguess);
    gates = [hh_inf(Vstar,1),hh_inf(Vstar,2),hh_inf(Vstar,3)];

    eqV15(k) = Vstar;
    eqN15(k) = gates(1);
    eqM15(k) = gates(2);
    eqH15(k) = gates(3);

    Xeq = [Vstar;gates(:)];
    J15 = numericalJacobian(@(x) hh_rhs(x,setfield(hh,'Iext',I)), ...
                            Xeq,[1e-4;1e-6;1e-6;1e-6]);
    maxReal15(k) = max(real(eig(J15)));

    Vguess = Vstar;
end

figure('Name','Q15 HH Equilibrium Voltage');
plot(currents15,eqV15,'o-','LineWidth',1.5);
xlabel('I_{ext} (\muA/cm^2)');
ylabel('Equilibrium V^* (mV)');
title('HH Equilibrium Voltage vs Applied Current');
grid on;

figure('Name','Q15 HH Stability');
plot(currents15,maxReal15,'o-','LineWidth',1.5);
yline(0,'--');
xlabel('I_{ext} (\muA/cm^2)');
ylabel('Maximum Re(\lambda)');
title('HH Equilibrium Stability');
grid on;

% Professor's I=9 example
I9idx = find(abs(currents15-9)<1e-12,1);
hh9 = hh;
hh9.Iext = 9;

Xeq9 = [eqV15(I9idx);
        eqN15(I9idx);
        eqM15(I9idx);
        eqH15(I9idx)];

XzeroEq = Xrest;

[tA,XA] = ode15s(@(t,X)hodgkinHuxley(t,X,hh9), ...
    [0 500],XzeroEq,odeOptsHH);

[tB,XB] = ode15s(@(t,X)hodgkinHuxley(t,X,hh9), ...
    [0 500],Xeq9,odeOptsHH);

figure('Name','Q15 HH Dynamics at I=9');
hA = plot(XA(:,1),XA(:,2),'LineWidth',1.5);
hold on;
hB = plot(XB(:,1),XB(:,2),'LineWidth',1.5);
hEQ = plot(Xeq9(1),Xeq9(2),'ko','MarkerFaceColor','k','MarkerSize',8);
xlabel('V (mV)');
ylabel('n');
title('HH Dynamics at I_{ext}=9 \muA/cm^2');
legend([hA,hB,hEQ], ...
    'Start at I=0 equilibrium','Start at I=9 equilibrium','I=9 equilibrium');
grid on;
hold off;

%% ============================================================
%  Q16 - REDUCED HH (V,n): m = m_inf(V), h = 0.8 - n
% =============================================================

hh16 = hh;
hh16.Iext = 9;

[an,bn,~,~,~,~] = hh_rates(Vrest);
n0red = an/(an+bn);

X0red = [Vrest;n0red];

[t16,X16] = ode15s(@(t,X) HH_reduced(t,X,hh16), ...
    [0 100],X0red,odeOptsHH);
%{
figure('Name','Q16 Reduced HH Action Potential');
plot(t16,X16(:,1),'LineWidth',1.5);
xlabel('Time (ms)');
ylabel('V (mV)');
title('Reduced HH Model');
grid on;
%}

% Full HH vs reduced HH under the same condition
hhCompare = hh;
hhCompare.Iext = 10;

X0full = [Vrest;nRest;mRest;hRest];
X0red  = [Vrest;nRest];

[tFull,XFull] = ode15s(@(t,X)hodgkinHuxley(t,X,hhCompare), ...
    [0 100],X0full,odeOptsHH);

[tRed,XRed] = ode15s(@(t,X)HH_reduced(t,X,hhCompare), ...
    [0 100],X0red,odeOptsHH);

figure('Name','Q16 Full vs Reduced HH');
plot(tFull,XFull(:,1),'LineWidth',1.5);
hold on;
plot(tRed,XRed(:,1),'--','LineWidth',1.5);
xlabel('Time (ms)');
ylabel('V (mV)');
title('Full HH vs Reduced HH');
legend('Full HH','Reduced HH');
grid on;
hold off;

% Reduced-model phase plane
figure('Name','Q16 Reduced HH Phase Plane');
plot(X16(:,1),X16(:,2),'LineWidth',1.5);
xlabel('V (mV)');
ylabel('n');
title('Reduced HH Phase Plane');
grid on;

%% ============================================================
%  Q17 - ANODE-BREAK EXCITATION
% =============================================================

hh17 = hh;
hh17.Iext = 0;

X0 = Xrest;

[t17,X17] = ode15s( ...
    @(t,X)hodgkinHuxleyPulse(t,X,hh17), ...
    [0 50],X0,odeOptsHH);

figure('Name','Q17 HH Anode-Break');
hV = plot(t17,X17(:,1),'LineWidth',1.5);
hold on;
hPulse = xline(20,'--');
xlabel('Time (ms)');
ylabel('V (mV)');
title('HH Anode-Break Excitation');
legend([hV,hPulse],'V(t)','Pulse end');
grid on;
hold off;

figure('Name','Q17 HH Gating Variables');
hN = plot(t17,X17(:,2),'LineWidth',1.3);
hold on;
hM = plot(t17,X17(:,3),'LineWidth',1.3);
hH = plot(t17,X17(:,4),'LineWidth',1.3);
hX = xline(20,'--');
xlabel('Time (ms)');
ylabel('Gating variables');
title('Gating Variables During Anode-Break Excitation');
legend([hN,hM,hH,hX],'n','m','h','Pulse end');
grid on;
hold off;

%% ============================================================
%  Q18 - REDUCED HH FAST SUBSYSTEM (V,m)
%  n and h fixed
% =============================================================

V_end = interp1(t17,X17(:,1),20);
n_end = interp1(t17,X17(:,2),20);
m_end = interp1(t17,X17(:,3),20);
h_end = interp1(t17,X17(:,4),20);

fprintf('\n================ Q18: Fast HH Subsystem ================\n');
fprintf('V_end = %.8f mV\n',V_end);
fprintf('n_end = %.8f\n',n_end);
fprintf('m_end = %.8f\n',m_end);
fprintf('h_end = %.8f\n',h_end);

n1 = nRest;
h1 = hRest;
n2 = n_end;
h2 = h_end;

Vnull18 = linspace(-80,60,3000);
mNull18 = hh_inf(Vnull18,2);

[mVnull1,mVnull2] = hhVmNullclines(Vnull18,hh,n1,h1,n2,h2);

% All equilibria
Fcase1 = @(V) ...
    -hh.GK*n1^4*(V-hh.EK) ...
    -hh.GNa*hh_inf(V,2)^3*h1*(V-hh.ENa) ...
    -hh.GL*(V-hh.EL);

Fcase2 = @(V) ...
    -hh.GK*n2^4*(V-hh.EK) ...
    -hh.GNa*hh_inf(V,2)^3*h2*(V-hh.ENa) ...
    -hh.GL*(V-hh.EL);

Vscan18 = linspace(-80,60,10000);
rootsCase1 = findAllRoots(Fcase1,Vscan18,arrayfun(Fcase1,Vscan18));
rootsCase2 = findAllRoots(Fcase2,Vscan18,arrayfun(Fcase2,Vscan18));

eqCase1 = [rootsCase1(:),hh_inf(rootsCase1(:),2)];
eqCase2 = [rootsCase2(:),hh_inf(rootsCase2(:),2)];

fprintf('\nCase 1: n,h fixed at resting values\n');
for k = 1:size(eqCase1,1)
    fprintf('Equilibrium %d: V = %.8f mV, m = %.8f\n', ...
        k,eqCase1(k,1),eqCase1(k,2));
end

fprintf('\nCase 2: n,h fixed at end-of-pulse values\n');
for k = 1:size(eqCase2,1)
    fprintf('Equilibrium %d: V = %.8f mV, m = %.8f\n', ...
        k,eqCase2(k,1),eqCase2(k,2));
end

% Simulations
X0case1 = [Vrest;mRest];
X0case2 = [V_end;m_end];

[tCase1,Xcase1] = ode15s( ...
    @(t,X)HH_Vm_reduced(t,X,hh,n1,h1), ...
    [0 50],X0case1,odeOptsHH);

[tCase2,Xcase2] = ode15s( ...
    @(t,X)HH_Vm_reduced(t,X,hh,n2,h2), ...
    [0 50],X0case2,odeOptsHH);

% Case 1 phase plane
figure('Name','Q18 Case 1');
hold on;
h1p = plot(Vnull18,mNull18,'LineWidth',2);
h2p = plot(Vnull18,mVnull1,'LineWidth',2);
h3p = plot(eqCase1(:,1),eqCase1(:,2),'ko','MarkerFaceColor','k','MarkerSize',8);
h4p = plot(Xcase1(:,1),Xcase1(:,2),'LineWidth',1.5);
xlabel('V (mV)');
ylabel('m');
title('Q18: V-m Phase Plane, n and h at Rest');
legend([h1p,h2p,h3p,h4p], ...
    'm-nullcline','V-nullcline','Equilibria','Trajectory');
xlim([-80 60]);
ylim([0 1]);
grid on;
hold off;

% Case 2 phase plane
figure('Name','Q18 Case 2');
hold on;
h1p = plot(Vnull18,mNull18,'LineWidth',2);
h2p = plot(Vnull18,mVnull2,'LineWidth',2);

if ~isempty(eqCase2)
    h3p = plot(eqCase2(:,1),eqCase2(:,2),'ko','MarkerFaceColor','k','MarkerSize',8);
end

h4p = plot(Xcase2(:,1),Xcase2(:,2),'LineWidth',1.5);
h5p = plot(X0case2(1),X0case2(2),'ks','MarkerFaceColor','k','MarkerSize',8);

xlabel('V (mV)');
ylabel('m');
title('Q18: V-m Phase Plane, n and h at End of Pulse');

if isempty(eqCase2)
    legend([h1p,h2p,h4p,h5p], ...
        'm-nullcline','V-nullcline','Trajectory','Initial condition');
else
    legend([h1p,h2p,h3p,h4p,h5p], ...
        'm-nullcline','V-nullcline','Equilibria', ...
        'Trajectory','Initial condition');
end

xlim([-80 60]);
ylim([0 1]);
grid on;
hold off;

fprintf('\nQ18 interpretation:\n');
fprintf(['Case 1 retains the low-voltage resting equilibrium.\n' ...
         'Case 2 loses the resting equilibrium near -60 mV after hyperpolarization.\n']);

%% ============================================================
%  HELPER FUNCTIONS
% =============================================================

function dxdt = morrisLecar(t,x,p)
    %#ok<INUSD>
    V = x(1);
    w = x(2);

    mInf = 0.5*(1+tanh((V-p.V1)/p.V2));
    wInf = 0.5*(1+tanh((V-p.V3)/p.V4));
    tauW = 1/cosh((V-p.V5)/p.V6);

    dV = ( ...
        p.Iext ...
        - p.gCa*mInf*(V-p.VCa) ...
        - p.gK*w*(V-p.VK) ...
        - p.gL*(V-p.VL)) / p.C;

    dw = p.phi*(wInf-w)/tauW;

    dxdt = [dV;dw];
end

function dx = ml_rhs(x,p,Ioverride)
    if nargin < 3
        I = p.Iext;
    else
        I = Ioverride;
    end

    V = x(1);
    w = x(2);

    mInf = 0.5*(1+tanh((V-p.V1)/p.V2));
    wInf = 0.5*(1+tanh((V-p.V3)/p.V4));
    tauW = 1/cosh((V-p.V5)/p.V6);

    dV = ( ...
        I ...
        - p.gCa*mInf*(V-p.VCa) ...
        - p.gK*w*(V-p.VK) ...
        - p.gL*(V-p.VL)) / p.C;

    dw = p.phi*(wInf-w)/tauW;

    dx = [dV;dw];
end

function [dV,dW] = ml_vector_field(V,W,p)
    mInf = 0.5*(1+tanh((V-p.V1)/p.V2));
    wInf = 0.5*(1+tanh((V-p.V3)/p.V4));
    tauW = 1./cosh((V-p.V5)/p.V6);

    dV = ( ...
        p.Iext ...
        - p.gCa*mInf.*(V-p.VCa) ...
        - p.gK*W.*(V-p.VK) ...
        - p.gL*(V-p.VL))/p.C;

    dW = p.phi*(wInf-W)./tauW;
end

function type = classify2DEquilibrium(lambda)
    r = real(lambda);

    if all(r < 0)
        if any(abs(imag(lambda)) > 1e-8)
            type = 'Stable focus';
        else
            type = 'Stable node';
        end
    elseif all(r > 0)
        if any(abs(imag(lambda)) > 1e-8)
            type = 'Unstable focus';
        else
            type = 'Unstable node';
        end
    elseif any(r > 0) && any(r < 0)
        type = 'Saddle';
    else
        type = 'Critical';
    end
end

function roots = findAllRoots(F,Vscan,Fscan)
    roots = [];

    for j = 1:numel(Vscan)-1
        fa = Fscan(j);
        fb = Fscan(j+1);

        if isfinite(fa) && isfinite(fb) && fa*fb < 0
            r = fzero(F,[Vscan(j),Vscan(j+1)]);

            if isempty(roots) || abs(r-roots(end)) > 1e-5
                roots(end+1) = r; %#ok<AGROW>
            end
        end
    end

    % Also catch roots that land numerically almost exactly on a scan point.
    closeIdx = find(abs(Fscan) < 1e-8);

    for j = closeIdx
        r = Vscan(j);

        if isempty(roots) || all(abs(r-roots) > 1e-5)
            roots(end+1) = r; %#ok<AGROW>
        end
    end

    roots = sort(roots);
end

function J = numericalJacobian(fun,x,h)
    x = x(:);
    n = numel(x);
    J = zeros(n);

    for k = 1:n
        xp = x;
        xm = x;
        xp(k) = xp(k)+h(k);
        xm(k) = xm(k)-h(k);

        J(:,k) = (fun(xp)-fun(xm))/(2*h(k));
    end
end

function [roots,weights] = mlEquilibria(I,p,Vscan)
    mInf = @(V) 0.5*(1+tanh((V-p.V1)/p.V2));
    wInf = @(V) 0.5*(1+tanh((V-p.V3)/p.V4));

    F = @(V) ...
        p.gCa*mInf(V).*(V-p.VCa) + ...
        p.gK*wInf(V).*(V-p.VK) + ...
        p.gL*(V-p.VL) - I;

    Fscan = arrayfun(F,Vscan);
    roots = findAllRoots(F,Vscan,Fscan);
    weights = wInf(roots);
end

function [value,isterminal,direction] = manifoldEvent(~,x,Vrange,wrange)
    V = x(1);
    w = x(2);

    value = min([ ...
        V-Vrange(1), ...
        Vrange(2)-V, ...
        w-wrange(1), ...
        wrange(2)-w]);

    isterminal = 1;
    direction = -1;
end

function [an,bn,am,bm,ah,bh] = hh_rates(V)

    % alpha_n
    an = -0.01*(V+50) ./ expm1(-(V+50)/10);

    % alpha_m
    am = -0.1*(V+35) ./ expm1(-(V+35)/10);

    % Handle the removable 0/0 cases
    an(abs(V+50) < 1e-8) = 0.1;
    am(abs(V+35) < 1e-8) = 1.0;

    % Remaining rate constants
    bn = 0.125*exp(-(V+60)/80);
    bm = 4*exp(-(V+60)/18);

    ah = 0.07*exp(-(V+60)/20);
    bh = 1 ./ (exp(-(V+30)/10) + 1);

end

function value = hh_inf(V,which)

    [an,bn,am,bm,ah,bh] = hh_rates(V);

    nInf = an./(an+bn);
    mInf = am./(am+bm);
    hInf = ah./(ah+bh);

    switch which
        case 1
            value = nInf;
        case 2
            value = mInf;
        case 3
            value = hInf;
        otherwise
            error('hh_inf: which must be 1, 2, or 3.');
    end
end

function dxdt = hodgkinHuxley(t,X,p)
    %#ok<INUSD>
    V = X(1);
    n = X(2);
    m = X(3);
    h = X(4);

    [an,bn,am,bm,ah,bh] = hh_rates(V);

    IK  = p.GK*n^4*(V-p.EK);
    INa = p.GNa*m^3*h*(V-p.ENa);
    IL  = p.GL*(V-p.EL);

    dV = (p.Iext-IK-INa-IL)/p.C;

    dn = p.phi*(an*(1-n)-bn*n);
    dm = p.phi*(am*(1-m)-bm*m);
    dh = p.phi*(ah*(1-h)-bh*h);

    dxdt = [dV;dn;dm;dh];
end

function dx = hh_rhs(X,p)
    dx = hodgkinHuxley(0,X,p);
end

function dxdt = hodgkinHuxleyPulse(t,X,p)

    V = X(1);
    n = X(2);
    m = X(3);
    h = X(4);

    [an,bn,am,bm,ah,bh] = hh_rates(V);

    if t < 20
        Iext = -3;
    else
        Iext = 0;
    end

    IK  = p.GK*n^4*(V-p.EK);
    INa = p.GNa*m^3*h*(V-p.ENa);
    IL  = p.GL*(V-p.EL);

    dV = (Iext-IK-INa-IL)/p.C;

    dn = p.phi*(an*(1-n)-bn*n);
    dm = p.phi*(am*(1-m)-bm*m);
    dh = p.phi*(ah*(1-h)-bh*h);

    dxdt = [dV;dn;dm;dh];
end

function dxdt = HH_reduced(~,X,p)
    V = X(1);
    n = X(2);

    [~,~,am,bm,~,~] = hh_rates(V);

    mInf = am/(am+bm);
    h = 0.8-n;

    IK  = p.GK*n^4*(V-p.EK);
    INa = p.GNa*mInf^3*h*(V-p.ENa);
    IL  = p.GL*(V-p.EL);

    dV = (p.Iext-IK-INa-IL)/p.C;

    [an,bn,~,~,~,~] = hh_rates(V);
    dn = an*(1-n)-bn*n;

    dxdt = [dV;dn];
end

function dxdt = HH_Vm_reduced(~,X,p,nFixed,hFixed)
    V = X(1);
    m = X(2);

    [~,~,am,bm,~,~] = hh_rates(V);

    IK  = p.GK*nFixed^4*(V-p.EK);
    INa = p.GNa*m^3*hFixed*(V-p.ENa);
    IL  = p.GL*(V-p.EL);

    dV = (-IK-INa-IL)/p.C;
    dm = am*(1-m)-bm*m;

    dxdt = [dV;dm];
end

function [mV1,mV2] = hhVmNullclines(V,p,n1,h1,n2,h2)

    mV1 = nan(size(V));
    mV2 = nan(size(V));

    for k = 1:numel(V)
        VV = V(k);

        numerator1 = ...
            -p.GK*n1^4*(VV-p.EK) ...
            -p.GL*(VV-p.EL);
        denominator1 = p.GNa*h1*(VV-p.ENa);

        if abs(denominator1) > 1e-12
            ratio = numerator1/denominator1;
            if ratio >= 0
                mv = ratio^(1/3);
                if mv >= 0 && mv <= 1
                    mV1(k) = mv;
                end
            end
        end

        numerator2 = ...
            -p.GK*n2^4*(VV-p.EK) ...
            -p.GL*(VV-p.EL);
        denominator2 = p.GNa*h2*(VV-p.ENa);

        if abs(denominator2) > 1e-12
            ratio = numerator2/denominator2;
            if ratio >= 0
                mv = ratio^(1/3);
                if mv >= 0 && mv <= 1
                    mV2(k) = mv;
                end
            end
        end
    end
end


