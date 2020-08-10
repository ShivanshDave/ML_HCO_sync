function [PLCorr, BLCorr] = script_f6_compare_omega_single_run(Omega, Plot, W)
%%
if exist('Plot','var') 
    Plot(Plot==0) = nan;
else
    Plot = nan(1,length(Omega)+1); 
end

% network design
% (1,2)(7,8): e-HCO, (3,4)(9,10): i-HOC, -- (7:10 - bursting neurons)
ML.net     = [ 2, 1, 4, 3, 0, 0, 8, 7,10, 9]; % Neuron connection
ML.synapse = [ 1, 1,-1,-1, 0, 0, 1, 1,-1,-1]; % Synapse : inhi/exci/off
ML.burst   = [ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]; % Burst mode : on/off
ML.Istim   = [90,90,90,90,85, 0, 0, 0, 0, 0]; % Current-clamp (nA)
ML.tmax = 5e3;  ML.dt = 0.1;
% Use demo-style plots in ML_network 
    
% Simulating for each Omega
det_flag=[1, zeros(1,length(Omega))];
Omega=[0, Omega]; % first sim is deterministic
BLCorr = nan(length(Omega),4, length(W));
PLCorr = nan(length(Omega),4);

for i=1:length(Omega)
    % Get raw data
    ML.system_size = Omega(i);
    ML.demo = 2;    
    ML.demo = Plot(i);
    
%      [V,t,spikes] = ML_network(2, det_flag(i), ML, 100); % <<NOTE : FIXED SEED IS SET HERE>>
    [V,t,spikes] = ML_network(2, det_flag(i), ML, 'shuffle'); % (demo, n, seed)

%     Vplot(i) = V; spikesplot(i) = spikes; tplot(i) = t;

    % find spike-time binless correlation
    dt = ML.dt;
    for j=1:length(W)
        w=W(j); plt=nan;
        % Non-bursting
        ttl="Plain-eHCO--BW="+w+"ms--Om="+Omega(i);
        BLCorr(i,1,j) = GetBinlessCorr(spikes{1,1}, spikes{2,1}, w,t,dt,plt,ttl);
        ttl="Plain-iHCO--BW="+w+"ms--Om="+Omega(i);
        BLCorr(i,2,j) = GetBinlessCorr(spikes{3,1}, spikes{4,1}, w,t,dt,plt,ttl);
        % bursting
        ttl="Bursting-eHCO--BW="+w+"ms--Om="+Omega(i);
        BLCorr(i,3,j) = GetBinlessCorr(spikes{7,1}, spikes{8,1}, w,t,dt,plt,ttl);
        ttl="Bursting-iHCO--BW="+w+"ms--Om="+Omega(i);
        BLCorr(i,4,j) = GetBinlessCorr(spikes{9,1}, spikes{10,1}, w,t,dt,plt,ttl);
    end
   
    % find Pearson's linear correlation coefficient (raw signal)
    % Non-bursting
    PLCorr(i,1) = corr(V(:,1),V(:,2),'type','Pearson');
    PLCorr(i,2) = corr(V(:,3),V(:,4),'type','Pearson');
    % bursting
    PLCorr(i,3) = corr(V(:,7),V(:,8),'type','Pearson');
    PLCorr(i,4) = corr(V(:,9),V(:,10),'type','Pearson');
end

%% Plot
txts=["Excitatory Synapses\newline Current Clamp", ...
    "Inhibitory Synapses\newline Current Clamp", ...
    "Excitatory Synapses\newline Internal Bursting", ...
    "Inhibitory Synapses\newline Internal Bursting"];
ticks={'E-HCO(1)','I-HCO(1)','E-HCO(2)','I-HCO(2)'};
Lgnd='Deterministic';
for i=2:length(Omega)
    Lgnd=[Lgnd,"Omega:"+Omega(i)];
end

cpsz = 14;
figure('Renderer', 'painters', 'Position', [0 0 700 700])
tiledlayout(2,2,'TileSpacing','none','Padding','compact');

% Plot Pearson's linear correlation coefficient
% subplot(2,2,1)
nexttile
bar(1:4,PLCorr','facealpha',0.2)
% yloc = min(1.1*min(PLCorr(:)),-1.2);
% for i=1:4
%     text(i,yloc,txts(i),'FontSize',9,....
%         'HorizontalAlignment','center','VerticalAlignment','bottom')
% end
xticklabels(ticks)
title({'Raw Traces', '(no gaussian smoothing)'},'fontsize', cpsz, 'FontWeight','Normal')
% xlabel('Type of network','fontsize', cpsz, 'FontWeight','Normal')
ylabel('Correlation Coefficient', 'fontsize', cpsz, 'FontWeight','Normal')
axis([-Inf Inf -1 1])

% Plot spike-time binless correlation
for j=1:length(W)
    BLC = BLCorr(:,:,j);
%     subplot(2,2,1+j)
    nexttile
    bar(1:4,BLC','facealpha',0.2)
%     yloc = 1.1*min(BLC(:));
%     for i=1:4
%         text(i,yloc,txts(i),'FontSize',9,....
%             'HorizontalAlignment','center',...
%             'VerticalAlignment','bottom')
%     end
    xticklabels(ticks)
    if j>1
        title("( W = "+W(j)+" ms )",'fontsize', cpsz, 'FontWeight','Normal')
        xlabel('Type of network','fontsize', cpsz, 'FontWeight','Normal')
        if j==2
            ylabel('Correlation Coefficient','fontsize', cpsz, 'FontWeight','Normal')
        end
    else
        title({"Binless Corr. of Spikes", "( W = "+W(j)+" ms )"},...
            'fontsize', cpsz, 'FontWeight','Normal')
    end
    axis([-Inf Inf -1 1])
end
lgd=legend(Lgnd,'Location','southwest');
lgd.NumColumns = 2;

end
