function [C,n] = script_f7_find_corr_all_data(n,seed)
%%
% BinWidth = 500ms, 250ms, 50ms, Raw
% SystemSize = 5000, 500, 50, Deterministic

if exist('seed','var'), rng(seed); end
if ~exist('n','var'), n=25; end

tic  
ML.tmax = 5e3;  ML.dt = 0.1; dt=0.1;
ML.demo = 2; % no-plots
ML.net = [n+1:2*n, 1:n]; % Neuron connections

% --------- Deterministic ---------------
det_flag=1; ML.system_size=0;

% ex-NB : Det
disp("DET : Start 1/4 "+toc);
ML.synapse = 1*ones(size(ML.net)); % Synapse : inhi/exci/off
ML.burst = 0*ones(size(ML.net));
ML.Istim = 90*ones(size(ML.net));  
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.eNB.detw50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.eNB.detw250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.eNB.detw500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.eNB.detpc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% in-NB : Det
disp("DET : Start 2/4 "+toc);
ML.synapse = -1*ones(size(ML.net)); % Synapse : inhi/exci/off
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.iNB.detw50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.iNB.detw250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.iNB.detw500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.iNB.detpc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% ex-B : Det 
disp("DET : Start 3/4 "+toc);
ML.synapse = 1*ones(size(ML.net)); % Synapse : inhi/exci/off
ML.burst = 1*ones(size(ML.net));
ML.Istim = 0*ones(size(ML.net));
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.eB.detw50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.eB.detw250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.eB.detw500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.eB.detpc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% in-B : Det 
disp("DET : Start 4/4 "+toc);
ML.synapse = -1*ones(size(ML.net)); % Synapse : inhi/exci/off
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.iB.detw50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.iB.detw250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.iB.detw500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.iB.detpc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% --------- Omega=5000 ---------------
det_flag = 0; ML.system_size=5000;

% ex-NB : Omega=5000
disp("5000 : Start 1/4 "+toc);
ML.synapse = 1*ones(size(ML.net)); % Synapse : inhi/exci/off
ML.burst = 0*ones(size(ML.net));
ML.Istim = 90*ones(size(ML.net));  
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.eNB.s5000w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.eNB.s5000w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.eNB.s5000w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.eNB.s5000pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% in-NB : Omega=5000
disp("5000 : Start 2/4 "+toc);
ML.synapse = -1*ones(size(ML.net)); % Synapse : inhi/exci/off
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.iNB.s5000w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.iNB.s5000w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.iNB.s5000w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.iNB.s5000pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% ex-B : Omega=5000 
disp("5000 : Start 3/4 "+toc);
ML.synapse = 1*ones(size(ML.net)); % Synapse : inhi/exci/off
ML.burst = 1*ones(size(ML.net));
ML.Istim = 0*ones(size(ML.net));
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.eB.s5000w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.eB.s5000w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.eB.s5000w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.eB.s5000pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% in-B : Omega=5000 
disp("5000 : Start 4/4 "+toc);
ML.synapse = -1*ones(size(ML.net)); % Synapse : inhi/exci/off
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.iB.s5000w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.iB.s5000w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.iB.s5000w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.iB.s5000pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% --------- Omega=500 ---------------
ML.system_size=500;

% ex-NB : Omega=500
disp("500 : Start 1/4 "+toc);
ML.synapse = 1*ones(size(ML.net)); % Synapse : inhi/exci/off
ML.burst = 0*ones(size(ML.net));
ML.Istim = 90*ones(size(ML.net));  
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.eNB.s500w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.eNB.s500w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.eNB.s500w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.eNB.s500pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% in-NB : Omega=500
disp("500 : Start 2/4 "+toc);
ML.synapse = -1*ones(size(ML.net)); % Synapse : inhi/exci/off
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.iNB.s500w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.iNB.s500w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.iNB.s500w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.iNB.s500pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% ex-B : Omega=500 
disp("500 : Start 3/4 "+toc);
ML.synapse = 1*ones(size(ML.net)); % Synapse : inhi/exci/off
ML.burst = 1*ones(size(ML.net));
ML.Istim = 0*ones(size(ML.net));
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.eB.s500w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.eB.s500w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.eB.s500w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.eB.s500pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% in-B : Omega=500 
disp("500 : Start 4/4 "+toc);
ML.synapse = -1*ones(size(ML.net)); % Synapse : inhi/exci/off
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.iB.s500w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.iB.s500w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.iB.s500w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.iB.s500pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% --------- Omega=50 ---------------
ML.system_size=50;

% ex-NB : Omega=50
disp("50 : Start 1/4 "+toc);
ML.synapse = 1*ones(size(ML.net)); % Synapse : inhi/exci/off
ML.burst = 0*ones(size(ML.net));
ML.Istim = 90*ones(size(ML.net));  
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.eNB.s50w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.eNB.s50w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.eNB.s50w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.eNB.s50pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% in-NB : Omega=50
disp("50 : Start 2/4 "+toc);
ML.synapse = -1*ones(size(ML.net)); % Synapse : inhi/exci/off
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.iNB.s50w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.iNB.s50w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.iNB.s50w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.iNB.s50pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% ex-B : Omega=50 
disp("50 : Start 3/4 "+toc);
ML.synapse = 1*ones(size(ML.net)); % Synapse : inhi/exci/off
ML.burst = 1*ones(size(ML.net));
ML.Istim = 0*ones(size(ML.net));
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.eB.s50w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.eB.s50w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.eB.s50w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.eB.s50pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

% in-B : Omega=50 
disp("50 : Start 4/4 "+toc);
ML.synapse = -1*ones(size(ML.net)); % Synapse : inhi/exci/off
[V,t,spikes] = ML_network(2, det_flag, ML, 'shuffle'); % (demo, n, seed)
for i=1:n
    C.iB.s50w50(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 50, t, dt);
    C.iB.s50w250(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 250, t, dt);
    C.iB.s50w500(i) = GetBinlessCorr(spikes{i,1}, spikes{i+n,1}, 500, t, dt);
    C.iB.s50pc(i) = corr(V(:,i),V(:,i+n),'type','Pearson');
end

disp(toc);
%--------------------

%%
% BinWidth = 500ms, 250ms, 50ms, Raw
% SystemSize = 5000, 500, 50, Deterministic
cpsz = 14;

if ~exist('n','var'), n=25; end
sz=55; alp=.5; jtr=.2;
clr={[0 0.4470 0.7410];
    [0.8500 0.3250 0.0980]; 
    [0.9290 0.6940 0.1250]; 
    [0.4940 0.1840 0.5560]};
type={'eNB';'iNB';'eB';'iB'};

figure('Renderer', 'painters', 'Position', [0 0 800 1000])
% sgtitle({'Correlation of HCO-neurons spike-train with different system size',...
%     "(For each cluster n="+n+" repeats ; Line connects the mean values)"}, 'Color',[.4 .4 1])

% PC traces
tiledlayout(2,2,'TileSpacing','none','Padding','compact')

nexttile
for i=1:4
    xval=[log10(5000)*ones(1,n),  log10(500)*ones(1,n),  log10(50)*ones(1,n), log10(1)*ones(1,n) ];
    yval=[C.(type{i}).s5000pc, C.(type{i}).s500pc, C.(type{i}).s50pc, C.(type{i}).detpc];
    ymean=[mean(C.(type{i}).s5000pc), mean(C.(type{i}).s500pc), mean(C.(type{i}).s50pc)];
    p(i)=scatter(xval,yval,sz,'MarkerEdgeColor', clr{i},'jitter', 'on', 'jitterAmount', jtr); 
    alpha(alp); hold on;
    plot([log10(5000),log10(500),log10(50)], ymean, '--', 'color', clr{i}, 'LineWidth', 2);
end
% xlabel('log_1_0(Omega)','fontsize', cpsz, 'FontWeight','Normal'); 
ylabel('Correlation Coeff','fontsize', cpsz, 'FontWeight','Normal');
xticks([0, log10(50), log10(500), log10(5000)])
xticklabels({'Det','log(50)','log(500)','log(5000)'})
yticks(-1:.2:1); axis([-1 4.5 -1.2 1.2]);
title({'Raw Traces','(no gaussian smoothing)'},'fontsize', 1.25*cpsz, 'FontWeight','Normal')

% Binless correlation
W=["w50","w250","w500"];
W2=[50,250,500];
for j=1:length(W)
    s={"s5000"+W(j), "s500"+W(j), "s50"+W(j), "det"+W(j)};
    nexttile
    for i=1:4
        xval=[log10(5000)*ones(1,n),  log10(500)*ones(1,n),  log10(50)*ones(1,n), log10(1)*ones(1,n) ];
        yval=[C.(type{i}).(s{1}), C.(type{i}).(s{2}), C.(type{i}).(s{3}), C.(type{i}).(s{4})];
        ymean=[mean(C.(type{i}).(s{1})), mean(C.(type{i}).(s{2})), mean(C.(type{i}).(s{3}))];
        p(i)=scatter(xval,yval,sz,'MarkerEdgeColor', clr{i},'jitter', 'on', 'jitterAmount', jtr); 
        alpha(alp); hold on;
        plot([log10(5000),log10(500),log10(50)], ymean, '--', 'color', clr{i}, 'LineWidth', 2);
    end
    if j>1 
        xlabel('log_1_0(Omega)','fontsize', cpsz, 'FontWeight','Normal'); 
        title({"( W = "+W2(j)+" ms )"}, 'fontsize', 1.25*cpsz, 'FontWeight','Normal')
    else
        title({"Binless Corr. of Spikes", "( W = "+W2(j)+" ms )"}, 'fontsize', 1.25*cpsz, 'FontWeight','Normal')
    end
    
    if j==2, ylabel('Correlation Coeff','fontsize', cpsz, 'FontWeight','Normal'); end
    xticks([0, log10(50), log10(500), log10(5000)])
    xticklabels({'Det','log(50)','log(500)','log(5000)'})
    yticks(-1:.2:1); axis([-1 4.5 -1.2 1.2]);
    
end

set(gcf,'unit','inches');
figure_size = get(gcf,'position');
leg_nm = {
    sprintf('Excitatory-HCO'),...
    sprintf('Inhibitory-HCO'),...
    sprintf('Excitatory-HCO\n(Bursting Mode)'),...
    sprintf('Inhibitory-HCO\n(Bursting Mode)')};
[leg,icons]=legend(p(1:4), leg_nm{1:4},'Location','NorthEastOutside','fontsize', .75*cpsz); 
title(leg,'HCO Type','fontsize', cpsz, 'FontWeight','Normal');
set(leg, 'unit', 'inches')
legend_size = get(leg, 'position');
figure_size(3) = figure_size(3) + legend_size(3);
set(gcf, 'position', figure_size)


end