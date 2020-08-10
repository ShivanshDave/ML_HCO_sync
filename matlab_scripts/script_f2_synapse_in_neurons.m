function script_f2_synapse_in_neurons(Omega)
%% 
if ~exist('seed','var'), seed=100; end 

% Network
% neuron    :<1,2,3, 4,5,6,7, 8>
ML.net     = [0,0,1, 2,0,0,5, 6];  % Neuron connection
ML.synapse = [0,0,1,-1,0,0,1,-1];  % Synapse : inhi/exci/off
ML.burst = [zeros(1,4), ones(1,4)];   % Burst mode : on/off
ML.tmax = 2e3; ML.dt = 0.1; ML.demo = 2; 

% Stoch
det_flag = 0;
ML.Istim = [80*ones(1,2), zeros(1,6)]; % Current-clamp (nA) 
ML.system_size = Omega*ones(1,8);
[V.st,t.st,spikes.st] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>
% ttl_st = {
%     {"Neuron-1","( Omega = "+ML.system_size(1)+", Istm = "+ML.Istim(1)+" )"},... 
%     {"Neuron-1","( Omega = "+ML.system_size(2)+", Istm = "+ML.Istim(2)+" )"},...
%     {"Neuron-2 w/ Exc-Synapse","( Omega = "+ML.system_size(3)+", Istm = "+ML.Istim(3)+" )"},...
%     {"Neuron-2 w/ Inh-Synapse","( Omega = "+ML.system_size(4)+", Istm = "+ML.Istim(4)+" )"},...
%     {"Bursting Neuron-1","( Omega = "+ML.system_size(5)+", Istm = "+ML.Istim(5)+" )"},...
%     {"Bursting Neuron-1","( Omega = "+ML.system_size(6)+", Istm = "+ML.Istim(6)+" )"},...
%     {"Neuron-2 w/ Exc-Synapse","( Omega = "+ML.system_size(7)+", Istm = "+ML.Istim(7)+" )"},...
%     {"Neuron-2 w/ Inh-Synapse","( Omega = "+ML.system_size(8)+", Istm = "+ML.Istim(8)+" )" } };
ttl_st = {
    {"Neuron-1 ( Istm = "+ML.Istim(1)+" )"},... 
    {"Neuron-1 ( Istm = "+ML.Istim(2)+" )"},...
    "Neuron-2 w/ Exc-Synapse",...
    "Neuron-2 w/ Inh-Synapse",...
    "Bursting Neuron-1",...
    "Bursting Neuron-1",...
    "Neuron-2 w/ Exc-Synapse",...
    "Neuron-2 w/ Inh-Synapse"  };



% Det
det_flag = 1;
ML.Istim = [95*ones(1,2), zeros(1,6)]; % Current-clamp (nA) 
ML.system_size = ones(1,8);
[V.dt,t.dt,spikes.dt] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>
% ttl_dt = {
%     {"Neuron-1","( Istm = "+ML.Istim(1)+" )"},... 
%     {"Neuron-1","( Istm = "+ML.Istim(2)+" )"},...
%     {"Neuron-2 w/ Exc-Synapse","( Istm = "+ML.Istim(3)+" )"},...
%     {"Neuron-2 w/ Inh-Synapse","( Istm = "+ML.Istim(4)+" )"},...
%     {"Bursting Neuron-1","( Istm = "+ML.Istim(5)+" )"},...
%     {"Bursting Neuron-1","( Istm = "+ML.Istim(6)+" )"},...
%     {"Neuron-2 w/ Exc-Synapse","( Istm = "+ML.Istim(7)+" )"},...
%     {"Neuron-2 w/ Inh-Synapse","( Istm = "+ML.Istim(8)+" )"} };
ttl_dt = {
    {"Neuron-1 ( Istm = "+ML.Istim(1)+" )"},... 
    {"Neuron-1 ( Istm = "+ML.Istim(2)+" )"},...
    "Neuron-2 w/ Exc-Synapse",...
    "Neuron-2 w/ Inh-Synapse",...
    "Bursting Neuron-1",...
    "Bursting Neuron-1",...
    "Neuron-2 w/ Exc-Synapse",...
    "Neuron-2 w/ Inh-Synapse" };


%% Plots
cpsz = 14;

%% Activation function plot
synThr = 5; synSlope = 0.5; v=0:.01:10; Syn = 0.5*(1+tanh((v-synThr)/synSlope));
figure('Renderer', 'painters', 'Position', [0 0 400 320]); 
tiledlayout(1,1,'TileSpacing','none','Padding','none');
nexttile; plot(v,Syn,'r','Linewidth',2.5); yticks([0 .5 1]); xticks([0:2:10]);
ylabel('Post-synaptic Activation','fontsize', cpsz, 'FontWeight','Normal'); 
xlabel('Pre-synaptic Membrane Voltage (mV)','fontsize', cpsz, 'FontWeight','Normal');
% title('Synaptic Activation function for ML Neurons')

%% Plot det
figure('Renderer', 'painters', 'Position', [0 0 800 400])
tiledlayout(4,2,'TileSpacing','none','Padding','none')

% sgtitle('Deterministic : Neuron(1)---o Neuron(2)','Color',[.4 .4 1],...
%     'fontsize', 1.5*cpsz, 'FontWeight','Normal')
for i=1:8
    nexttile
    plot(t.dt, V.dt(:,i)); hold on
    if i>4
        scatter(spikes.dt{i,3}, 1.6*spikes.dt{i,2}, 25, 'v', 'filled')
        axis([-Inf Inf -80 45])
        yticks(-80:40:40)
    else
        scatter(spikes.dt{i,3}, 1.4*spikes.dt{i,2}, 25, 'v', 'filled')
        axis([-Inf Inf -80 60])
    end
    title(ttl_dt{i},'fontsize', cpsz, 'FontWeight','Normal');
    if i==7, ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal'); end
    if i>6, xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal'); end
end


%% Plot stoch
figure('Renderer', 'painters', 'Position', [0 0 800 400])
tiledlayout(4,2,'TileSpacing','none','Padding','none')

% sgtitle('Stochastic : Neuron(1)---o Neuron(2)','Color',[.4 .4 1],...
%     'fontsize', 1.5*cpsz, 'FontWeight','Normal')
for i=1:8
    nexttile
    plot(t.st, V.st(:,i)); hold on;
    if i>4
        scatter(spikes.st{i,3}, 1.6*spikes.st{i,2}, 25, 'v', 'filled')
        axis([-Inf Inf -80 45])
        yticks([-80:40:40])
    else
        scatter(spikes.st{i,3}, 1.4*spikes.st{i,2}, 25, 'v', 'filled')
        axis([-Inf Inf -80 60])
    end
    title(ttl_st{i},'fontsize', cpsz, 'FontWeight','Normal');
    if i>6, xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal'); end
    if i==7, ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal'); end
end

end
