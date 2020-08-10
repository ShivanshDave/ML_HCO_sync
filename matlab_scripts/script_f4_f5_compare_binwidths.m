function script_f4_f5_compare_binwidths(W)
%% Figure-4-5
det=0;      
n.system_size=1e3;  
n.burst=[0,0];  
n.Istim=[85,85];
n.dt=0.1;   
n.demo=2;
n.tmax=5e3;  
n.net     = [2,1]; 
n.synapse = [0,0];

[V,t,spikes]=ML_network(2,det,n,100);

% HCO Synchronization measure
dt=0.1;
ind_spk_1 = spikes{1,1};
ind_spk_2 = spikes{2,1};
C=nan(size(W));
Wstr = "w"+W';

for i=1:length(W)
    [C(i),G.(Wstr(i))] = GetBinlessCorr(ind_spk_1,ind_spk_2,W(i),t,dt);
end
disp("F5_GSCorr :::: for W="+W+"ms; C="+C);

%% Plot F4
cpsz = 14;
figure('Renderer', 'painters', 'Position', [0 0 800 200])
tiledlayout(1,1,'TileSpacing','none','Padding','none')
nexttile
plot(t,V); hold on; 
scatter(spikes{1,3}, 1.4*spikes{1,2}, 55, 'bv', 'filled');
scatter(spikes{2,3}, 1.4*spikes{2,2}, 55, 'rv', 'filled');
ylabel('Vm (mV)','fontsize', cpsz, 'FontWeight','Normal'); 
xlabel('time (ms)','fontsize', cpsz, 'FontWeight','Normal'); 
title({"Sample ML neurons (un-connected)","( Omega = "+n.system_size+", Istim = "+n.Istim(1)+")"},...
    'fontsize', 1*cpsz, 'FontWeight','Normal')

%% Plot F5
cpsz = 14;
figure('Renderer', 'painters', 'Position', [0 0 1000 650])
tiledlayout(3,3,'TileSpacing','none','Padding','none')

for i=1:3
    gs = G.(Wstr(i));
    nexttile([1 2])
%     subplot(3,3,(i*3-2):(i*3-1))
    plot(gs.t,gs.V1); hold on; plot(gs.t,gs.V2);
    if i==1
        title({"Gaussian Smoothed Neural spikes","( W = "+gs.W+" ms, S.D. = "+gs.sd+" ms)"},...
            'fontsize', cpsz, 'FontWeight','Normal')
    else
        title("( W = "+gs.W+" ms, S.D. = "+gs.sd+" ms)",'fontsize', cpsz, 'FontWeight','Normal')
    end
    if i==3
        xlabel('time(ms)','fontsize', cpsz, 'FontWeight','Normal'); 
        ylabel('Amp','fontsize', cpsz, 'FontWeight','Normal');
    end
    
%     subplot(3,3,i*3)
    nexttile
    plot(gs.tw,gs.win);
    if i==1
        title(["Gaussian Window (+/- 6*SD)"," "],...
            'fontsize', cpsz, 'FontWeight','Normal')
    else
        title(" ",'fontsize', cpsz, 'FontWeight','Normal')
    end
    if i==3
        xlabel('time(ms)','fontsize', cpsz, 'FontWeight','Normal'); 
    end

end
end
