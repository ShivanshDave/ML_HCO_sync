function script_f1_compare_current_omega_single_neuron(Istim, Omega, seed)
%%

% Make all neurons use same seed
if ~exist('seed','var'), seed=100; end 

% Simulation
ML.tmax = 2e3;  ML.dt = 0.1; ML.demo = 2; 
nI = length(Istim); nO = length(Omega);

% --- Plain neuron ---
% Stoch
n = nI*nO; % Total neurons
det_flag = zeros(1,n);
ML.net = zeros(1,n);  % Neuron connection
ML.synapse = zeros(1,n); % Synapse : inhi/exci/off
ML.burst = zeros(1,n); % Burst mode : on/off
ML.Istim = repmat(Istim,1,nO); % Current-clamp (nA)
ML.system_size = repmat(Omega,nI,1); 
ML.system_size = ML.system_size(:)';
[V.st,t.st,spikes.st] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>

% Det
n = nI; % Total neurons
det_flag = ones(1,n);
ML.net = zeros(1,n);  % Neuron connection
ML.synapse = zeros(1,n); % Synapse : inhi/exci/off
ML.burst = zeros(1,n); % Burst mode : on/off
ML.Istim = Istim; % Current-clamp (nA)
ML.system_size = 1; 
[V.dt,t.dt,spikes.dt] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>

%
% --- Bursting neuron --- 
% Stoch
n = nO; % Total neurons
det_flag = zeros(1,n);
ML.net = zeros(1,n);  % Neuron connection
ML.synapse = zeros(1,n); % Synapse : inhi/exci/off
ML.burst = ones(1,n); % Burst mode : on/off
ML.Istim = zeros(1,n); % Current-clamp (nA)
ML.system_size = Omega; 
[V.stB,t.stB,spikes.stB] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>

% Det
n = 1; % Total neurons
det_flag = ones(1,n);
ML.net = zeros(1,n);  % Neuron connection
ML.synapse = zeros(1,n); % Synapse : inhi/exci/off
ML.burst = ones(1,n); % Burst mode : on/off
ML.Istim = 0; % Current-clamp (nA)
ML.system_size = 1; 
[V.dtB,t.dtB,spikes.dtB] = ML_network(2, det_flag, ML, seed); % <<NOTE : FIXED SEED IS SET HERE>>


%% Plot X-> Det/Omega ; Y -> Current
cpsz = 14;
nX = nO + 1; nY = nI;

% --- Plain neuron ----
figure('Renderer', 'painters', 'Position', [0 0 1100 600])
tiledlayout(nX,nY,'TileSpacing','none','Padding','none')

% sgtitle('ML neurons with various Omega and Current', 'Color',[.4 .4 1],...
%     'fontsize', 1.5*cpsz, 'FontWeight','Normal')

%%%%
for i=1:nY
    nexttile
    plot(t.dt, V.dt(:,i))
    hold on
    scatter(spikes.dt{i,3}, 1.4*spikes.dt{i,2}, 55, 'v', 'filled')
    if i==1
        title({"Deterministic", "( Stim.Current = "+Istim(i)+" )"}, ...
            'fontsize', cpsz, 'FontWeight','Normal');
    else
        title({"( Stim.Current = "+Istim(i)+" )"},...
            'fontsize', cpsz, 'FontWeight','Normal'); 
    end
    ylabel('Vm (mV)', 'fontsize', cpsz)
    axis([-Inf Inf -60 60])
    if i==nY
        xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal')
    end
    
    for i2=1:nO
        j=i+(nY*(i2-1));
        nexttile
%         subplot(nX,nY,1+j+(nX*(i2-1)))
        plot(t.st, V.st(:,j))
        hold on
        scatter(spikes.st{j,3}, 1.4*spikes.st{j,2}, 55, 'v', 'filled')
        if i==1 
            title({"Stochastic : Omega="+Omega(i2), "( Stim.Current = "+Istim(i)+" )"},...
                'fontsize', cpsz, 'FontWeight','Normal');
        else
            title({"( Stim.Current = "+Istim(i)+" )"}, ...
                'fontsize', cpsz, 'FontWeight','Normal'); 
        end
        axis([-Inf Inf -60 60])
        if i==nY
            xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal')
        end
    end
    
       
end

 %% --- Bursting neuron ----
figure('Renderer', 'painters', 'Position', [0 0 500 500])
tiledlayout(nY,1,'TileSpacing','none','Padding','none')

% sgtitle('ML neurons with various Omega and Current', 'Color',[.4 .4 1],...
%     'fontsize', 1.5*cpsz, 'FontWeight','Normal')

%%%%
for i=1:1
    nexttile
    plot(t.dtB, V.dtB(:,i))
    hold on
    scatter(spikes.dtB{i,3}, 1.4*spikes.dtB{i,2}, 55, 'v', 'filled')
    title("Deterministic",'fontsize', cpsz, 'FontWeight','Normal');
    ylabel('Vm (mV)', 'fontsize', cpsz); yticks(-40:40:40)
    axis([-Inf Inf -60 40])
    for i2=1:nO
        j=i2;
        nexttile
        plot(t.stB, V.stB(:,j))
        hold on
        scatter(spikes.stB{j,3}, 1.4*spikes.stB{j,2}, 55, 'v', 'filled')
        title("Stochastic : Omega="+Omega(i2),...
                'fontsize', cpsz, 'FontWeight','Normal');
        axis([-Inf Inf -65 40]);yticks(-40:40:40)
        ylabel('Vm (mV)', 'fontsize', cpsz); yticks(-40:40:40)
    end
    xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal')
end

end
