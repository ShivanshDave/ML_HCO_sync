function script_f3_compare_omega_hco_types(Omega)
%%
% Make all neurons use same seed
if ~exist('seed','var'), seed=100; end 

% Simulation
ML.tmax = 2e3;  ML.dt = 0.1; ML.demo = 2; 
nO = length(Omega);

% --- Plain neuron ---
% Stoch
det_flag = 0;
ML.net     = [2,1, 4, 3, 6,5, 8, 7, 10,9,12,11];  % Neuron connection
ML.synapse = [1,1,-1,-1, 1,1,-1,-1, 1,1,-1,-1]; % Synapse : inhi/exci/off
ML.burst = 0; % Burst mode : on/off
ML.Istim = 95; % Current-clamp (nA)
ML.system_size = repmat(Omega,4,1); 
ML.system_size = ML.system_size(:)';
[V.st,t.st,spikes.st] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>
omga_st = ML.system_size;

% Det
det_flag = 1;
ML.net = [2,1,4,3];  % Neuron connection
ML.synapse = [1,1,-1,-1]; % Synapse : inhi/exci/off
ML.burst = 0; % Burst mode : on/off
ML.Istim = 95; % Current-clamp (nA)
ML.system_size = 1; 
[V.dt,t.dt,spikes.dt] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>

%
% --- Bursting neuron --- 
% Stoch
det_flag = 0;
ML.net     = [2,1, 4, 3, 6,5, 8, 7, 10,9,12,11];  % Neuron connection
ML.synapse = [1,1,-1,-1, 1,1,-1,-1, 1,1,-1,-1]; % Synapse : inhi/exci/off
ML.burst = 1; % Burst mode : on/off
ML.Istim = 95; % Current-clamp (nA)
ML.system_size = omga_st; 
[V.stB,t.stB,spikes.stB] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>

% Det
det_flag = 1;
ML.net = [2,1,4,3];  % Neuron connection
ML.synapse = [1,1,-1,-1]; % Synapse : inhi/exci/off
ML.burst = 1; % Burst mode : on/off
ML.Istim = 95; % Current-clamp (nA)
ML.system_size = 1; 
[V.dtB,t.dtB,spikes.dtB] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>

%% --- Plain neuron ----
cpsz = 14;
figure('Renderer', 'painters', 'Position', [0 0 800 400])
tiledlayout(4,2,'TileSpacing','none','Padding','none')

% sgtitle({'Comparing Omega and non-bursting HCO types',...
%     "(Stim.Current="+ML.Istim+")"},'Color',[.4 .4 1])
% Plot det
for i=1:2
    nexttile
    n1=2*i-1; n2=2*i;
    plot(t.dt, V.dt(:,n1)); hold on; plot(t.dt, V.dt(:,n2));
    scatter(spikes.dt{n1,3}, 1.4*spikes.dt{n1,2}, 25, 'bv', 'filled')
    scatter(spikes.dt{n2,3}, 1.4*spikes.dt{n2,2}, 25, 'rv', 'filled')
    if i==1
        title(["Excitatory-HCO","( Deterministic )"],'fontsize', cpsz, 'FontWeight','Normal'); 
    else
        title(["Inhibitory-HCO","( Deterministic )"],'fontsize', cpsz, 'FontWeight','Normal'); 
    end
%     ylabel('mV','fontsize', cpsz, 'FontWeight','Normal')
    axis([0 1e3 -60 60])
end

% Plot stoch
for i=1:2*nO
    nexttile
    n1=2*i-1; n2=2*i;
    plot(t.st, V.st(:,n1)); hold on; plot(t.st, V.st(:,n2));
    scatter(spikes.st{n1,3}, 1.4*spikes.st{n1,2}, 25, 'bv', 'filled')
    scatter(spikes.st{n2,3}, 1.4*spikes.st{n2,2}, 25, 'rv', 'filled')
    title("( Stoch: Omga = "+omga_st(2*i)+" )",'fontsize', cpsz, 'FontWeight','Normal'); 
    if i>4, xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal'); end
    if i==5, ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal'); end
    axis([0 1e3 -60 60])
end
% legend({'Neuron-1','Neuron-2'})


%% --- Bursting neuron ----
figure('Renderer', 'painters', 'Position', [0 0 800 400])
tiledlayout(4,2,'TileSpacing','none','Padding','none')

% sgtitle('Comparing Omega and bursting HCO types', 'Color',[.4 .4 1])
% Plot det
for i=1:2
    nexttile
    n1=2*i-1; n2=2*i;
    plot(t.dtB, V.dtB(:,n1)); hold on; plot(t.dtB, V.dtB(:,n2));
    scatter(spikes.dtB{n1,3}, 1.4*spikes.dtB{n1,2}, 25, 'bv', 'filled')
    scatter(spikes.dtB{n2,3}, 1.4*spikes.dtB{n2,2}, 25, 'rv', 'filled')
     if i==1
        title(["Excitatory-HCO","( Deterministic )"],'fontsize', cpsz, 'FontWeight','Normal'); 
    else
        title(["Inhibitory-HCO","( Deterministic )"],'fontsize', cpsz, 'FontWeight','Normal'); 
    end
    axis([-Inf Inf -80 40])
    yticks(-80:40:40)
end

% Plot stoch
for i=1:2*nO
    nexttile
    n1=2*i-1; n2=2*i;
    plot(t.stB, V.stB(:,n1)); hold on; plot(t.stB, V.stB(:,n2));
    scatter(spikes.stB{n1,3}, 1.4*spikes.stB{n1,2}, 25, 'bv', 'filled')
    scatter(spikes.stB{n2,3}, 1.4*spikes.stB{n2,2}, 25, 'rv', 'filled')
    title("( Stoch: Omga = "+omga_st(2*i)+" )",'fontsize', cpsz, 'FontWeight','Normal'); 
    if i>4, xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal'); end
    if i==5, ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal'); end
    axis([-Inf Inf -80 40])
    yticks(-80:40:40)
end


end
