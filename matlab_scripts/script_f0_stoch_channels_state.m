%% Deterministic
det=1;      
n.system_size=1;  
n.burst=[0,1];  
n.Istim=[95,0];
n.dt=0.1;   
n.demo=2;
n.tmax=2e3;  
n.net     = [0,0]; 
n.synapse = [0,0];
[V,t,spikes,X]=ML_network(2,det,n,100);

%% Plot
cpsz = 14;
figure('Renderer', 'painters', 'Position', [0 0 1000 500])
t1=tiledlayout(3,2,'TileSpacing','none','Padding','none');
% title(t1, 'Deterministic Simulation')

% Plot plain
nexttile(1); plot(t,V(:,1)); hold on; 
scatter(spikes{1,3}, 1.4*spikes{1,2}, 25, 'rv', 'filled');
title("ML neuron ( I = "+n.Istim(1)+" )",...
    'fontsize', cpsz, 'FontWeight','Normal')
ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal');

nexttile(3); plot(t,X.M(:,1)); axis([-Inf Inf 0 1]);
ylabel('M_{Open}/Omega','fontsize', cpsz, 'FontWeight','Normal');

nexttile(5); plot(t,X.N(:,1)); axis([-Inf Inf 0 1]);
ylabel('N_{Open}/Omega','fontsize', cpsz, 'FontWeight','Normal');
xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal'); 

% Plot Bursting
nexttile(2); plot(t,V(:,2)); hold on; 
scatter(spikes{2,3}, 1.6*spikes{2,2}, 25, 'rv', 'filled');
title("Bursting ML neuron",...
    'fontsize', cpsz, 'FontWeight','Normal')
% ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal');

nexttile(4); plot(t,X.M(:,2)); axis([-Inf Inf 0 1]);
% ylabel('M_{Open}/Omega','fontsize', cpsz, 'FontWeight','Normal');

nexttile(6); plot(t,X.N(:,2)); axis([-Inf Inf 0 1]);
% ylabel('N_{Open}/Omega','fontsize', cpsz, 'FontWeight','Normal');
xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal');

%% Stochastic
det=0;      
n.system_size=1e3;  
n.burst=[0,1];  
n.Istim=[95,0];
n.dt=0.1;   
n.demo=2;
n.tmax=2e3;  
n.net     = [0,0]; 
n.synapse = [0,0];
[V,t,spikes,X]=ML_network(2,det,n,100);

%%
cpsz = 14;
figure('Renderer', 'painters', 'Position', [0 0 1000 500])
t1=tiledlayout(3,2,'TileSpacing','none','Padding','none');

% Plot plain
nexttile(1); plot(t,V(:,1)); hold on; 
scatter(spikes{1,3}, 1.4*spikes{1,2}, 25, 'rv', 'filled');
title("ML neuron ( Omega = "+n.system_size+", I = "+n.Istim(1)+" )",...
    'fontsize', cpsz, 'FontWeight','Normal')
ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal');

nexttile(3); plot(t,X.M(:,1)); axis([-Inf Inf 0 1]);
ylabel('M_{Open}/Omega','fontsize', cpsz, 'FontWeight','Normal');

nexttile(5); plot(t,X.N(:,1)); axis([-Inf Inf 0 1]);
ylabel('N_{Open}/Omega','fontsize', cpsz, 'FontWeight','Normal');
xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal'); 

% Plot Bursting
nexttile(2); plot(t,V(:,2)); hold on; 
scatter(spikes{2,3}, 1.6*spikes{2,2}, 25, 'rv', 'filled');
title("Bursting ML neuron ( Omega = "+n.system_size+" )",...
    'fontsize', cpsz, 'FontWeight','Normal')
% ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal');

nexttile(4); plot(t,X.M(:,2)); axis([-Inf Inf 0 1]);
% ylabel('M_{Open}/Omega','fontsize', cpsz, 'FontWeight','Normal');

nexttile; plot(t,X.N(:,2)); axis([-Inf Inf 0 1]);
% ylabel('N_{Open}/Omega','fontsize', cpsz, 'FontWeight','Normal');
xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal');