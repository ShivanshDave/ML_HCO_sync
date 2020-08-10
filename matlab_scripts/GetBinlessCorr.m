function [C,gs] = GetBinlessCorr(spk_ind_1,spk_ind_2,W,t,dt,plt,ttl)
%% 

gs.t=t;
gs.W=W;

% Get spike train in bool values
sz = size(gs.t);
if isempty(spk_ind_1) || isempty(spk_ind_2)
    C=0; return;
end
spike_train_1 = false(sz);
spike_train_1(spk_ind_1) = true;
spike_train_2 = false(sz);
spike_train_2(spk_ind_2) = true;

% Defining the gaussian kernel
gs.sd = gs.W/sqrt(12)/dt;   % SAMPLES; SD of Gaussian Curve
t_window = -gs.sd*6:1:gs.sd*6; % 6*sd window
gs.win = 1/(sqrt(2*pi)*gs.sd) * exp(-(t_window.^2)/(2*gs.sd^2));
gs.tw = t_window*dt;
gs.sd = gs.sd*dt;

% Convolute it with a gaussian kernel
gs.V1 = conv(spike_train_1,gs.win,'same');
gs.V2 = conv(spike_train_2,gs.win,'same');

% Find correlation coeff for gaussian smoothed signal
meanV1 = mean(gs.V1); meanV2 = mean(gs.V2);
Cov12 = (gs.V1 - meanV1).*(gs.V2 - meanV2); 
Var1 = (gs.V1 - meanV1).^2; 
Var2 = (gs.V2 - meanV2).^2;
C = sum(Cov12) / ( sqrt(sum(Var1))*sqrt(sum(Var2)) );

if exist('plt','var')
    if plt==0 || plt==2
        figure;
        plot(gs.tw,gs.win);
        title("Gaussian Window (+/- 6*SD) W="+gs.W+"ms");
        xlabel('time(ms)'); ylabel('Amp');
    end
    if plt==1 || plt==2
        figure;
        plot(gs.t,gs.V1); hold on; plot(gs.t,gs.V2);
        xlabel('time(ms)'); ylabel('Amp');
        title("Gaussian Smoothed Neural spikes  W="+gs.W+"ms")
        legend({'Neuron-1','Neuron-2'})
    end
    if plt==5
        figure;
        plot(gs.t,gs.V1); hold on; plot(gs.t,gs.V2);
        xlabel('time(ms)'); ylabel('Amp');
        title(ttl)
        legend({'Neuron-1','Neuron-2'})
    end
    
end

end
