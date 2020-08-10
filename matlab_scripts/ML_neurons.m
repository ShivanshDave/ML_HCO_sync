function [Vm,nK,nCa,Iint] = ML_neurons(det,dt,Istim,Iint,omega,nK,nCa,Vm,PreSynVm,synapse,burst) 
%% Each passing var (except dt) can be an array

% --- Setting up neuron cellular parameters ---
% Model parameters for each neurons (from (Anderson et. al., 2015))
vCa = 120;                % Rev.Pot for Calcium channels
gCa = 4.4;                % Calcium conductance
vK  = -84;                % Rev.Pot for Potassium channels
gK  =   8;                % Potassium conductance
vL  = -60;                % Rev.Pot for leak channels
gL  =   2;                % Leak channels conductance
Cm  =  20;                % Membrane Conductance

% Defining synapse  (from Z.Yu, PJT, 2020)
vSyn = 100*synapse;         % Make Rev.pot. neg for inhibitory synapse
gSyn = 1*abs(synapse);    % Make conductance zero for no-synapse
synThr = 5; synSlope = 0.5;         % sigmoid activation fxn parameters 
Syn = @(v) 0.5*(1+tanh((v-synThr)/synSlope));
% Syn = @(v) 1./(1 + exp(-(synSlope/2).*(v-synThr)));

% Bursting 
% Slow current feedback for bursting
epsi = .01; v0 = -26;
dIint = @(v) epsi*(v0-v).*burst;

% K+ ion channel parameters
% vc = 12; vd = 17.4; phi_n = 0.23; % Adapted for bursting (from mlsqr)
% vc = 2; vd = 30; phi_n = 0.04; % For simple (non-bursting) neurons
vc = 2 + burst*10; vd = 30 - burst*12.6; phi_n = 0.04 + burst*0.19;
xi_n    = @(v) (v-vc)./vd;                   % scaled argument for n-gate input
ninf    = @(v) 0.5*(1+tanh(xi_n(v)));       % n-gate activation function
tau_n   = @(v) 1./(phi_n.*cosh(xi_n(v)/2));  % n-gate activation t-const
alpha_n = @(v) ninf(v)./tau_n(v);           % per capita opening rate
beta_n  = @(v) (1-ninf(v))./tau_n(v);       % per capita closing rate

% Ca2+ ion channel parameters
va = -1.2; vb = 18; phi_m = 2;
xi_m    = @(v) (v-va)/vb;                   % scaled argument for m-gate input
minf    = @(v) 0.5*(1+tanh(xi_m(v)));       % m-gate activation function
tau_m   = @(v) 1./(phi_m*cosh(xi_m(v)/2));  % m-gate time constant
alpha_m = @(v) minf(v)./tau_m(v);           % per capita opening rate
beta_m  = @(v) (1-minf(v))./tau_m(v);       % per capita closing rate

% --- Setting up Chemical Langevin Equation ---
tau = dt;
% stochasting K+ channels
omega_K = omega; 
No = nK.*omega_K;
dWk = randn(size(No)); % N(0,1)
h1k = @(v,X) alpha_n(v).*(omega_K-X);
h2k = @(v,X) beta_n(v).*(X);
N_next = @(v,X) X + tau.*(h1k(v,X) - h2k(v,X)) ...
    + sqrt(tau.*(h1k(v,X) + h2k(v,X))).*dWk;

% stochasting Ca2+ channels
omega_Ca = omega;  
Mo = nCa.*omega_Ca;
dWca = randn(size(Mo)); % N(0,1)
h1ca = @(v,X) alpha_m(v).*(omega_Ca-X);
h2ca = @(v,X) beta_m(v).*(X);
M_next = @(v, X) X + tau.*(h1ca(v,X) - h2ca(v,X)) ...
    + sqrt(tau.*(h1ca(v,X) + h2ca(v,X))).*dWca;

% Enforce valid limits
lim = @(n,Low,Max)(n>Max).*Max + (n<Low).*Low + (n<=Max & n>=Low).*n;

% Simulation using CLE or Deterministic eqns

% Update ion channels counts and fraction for each neuron
if det==1                   % ( FOR DETERMINISTIC SIM )
    nK = nK + dt*(ninf(Vm)-nK)./tau_n(Vm);
    nCa = minf(Vm);
else                        % ( FOR STOCHASTIC SIM )
    No = N_next(Vm,No);     
    Mo = M_next(Vm,Mo);  
    No = lim(No,zeros(size(No)),omega_K);
    Mo = lim(Mo,zeros(size(Mo)),omega_Ca);
    nK  = No./omega_K;
    nCa = Mo./omega_Ca;
end
Iint = Iint + dt*dIint(Vm);
    
% Update membrane voltage for each neuron
dV = (dt/Cm).*( Istim + Iint ...
    - gL*(Vm-vL) ...
    - gK*nK.*(Vm-vK) ...
    - gCa*nCa.*(Vm-vCa) ...
    - gSyn.*Syn(PreSynVm).*(Vm-vSyn) );  
Vm = Vm + dV;

end
