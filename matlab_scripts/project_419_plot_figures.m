tic

%% (Fig-0) Show neuron and ion channels
script_f0_stoch_channels_state;

%% (Fig-1) Simple ML Neuron Firing vs {SysSize, Current}
% F1-A : ML Plain {Omega/Det, Istm}
% F1-B : ML Bursting {Omega/Det}

Istim=[55,80,95,100];
Omega=[5000,500,50]; % Plus det
script_f1_compare_current_omega_single_neuron(Istim, Omega);
disp("f1"+toc)

%% (Fig-2) Synapses {ex/in} in ML Neuron & "excitable" state 
% F2-A : Stoch
% F2-B : Det
Omega=500;
script_f2_synapse_in_neurons(Omega);
disp("f2"+toc)

%% (Fig-3) Compare HCO types for system size
% F1-A : ML Plain {Omega/Det}
% F1-B : ML Bursting {Omega/Det}
Omega=[50,500,5000];  % (Length should be 3)
script_f3_compare_omega_hco_types(Omega);
disp("f2"+toc)

%% (Fig-4,5) Demo - Compare BinWidths
W=[50, 250, 500];
script_f4_f5_compare_binwidths(W);
disp("f4-5"+toc)

%% (Fig-Fx1,6) OLD - Single Run - Compare across diff omega
W=[50, 250, 500];
Omega = [5e4, 5e3, 5e2, 5e1, 5e0];  % Varying noise if stochastic sim
Plot = [1, 1,   1,   1,   1,   1]; % plot the det. and each omega (Length=omega+1)
script_f6_compare_omega_single_run(Omega,Plot,W);
disp("f6"+toc)

%% (Fig-7) Collect All Data (n runs)
[C,n] = script_f7_find_corr_all_data(25,100);
disp("f7"+toc)

%% ----------- TEMP -----------------
% script_fx1_single_HCO;