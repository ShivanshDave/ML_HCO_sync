%% HCO 
det=0;      
n.system_size=50;  
n.burst=[0,0];  
n.Istim=[100,95];
n.dt=0.1;   
n.demo=2;
n.tmax=5e3;  
n.net     = [0,0]; 
n.synapse = [0,0];

[V,t,spikes]=ML_network(2,det,n,100);

figure
plot(t,V); hold on; 
scatter(spikes{1,3}, 1.4*spikes{1,2}, 55, 'bv', 'filled');
scatter(spikes{2,3}, 1.4*spikes{2,2}, 55, 'rv', 'filled');
ylabel('Vm (mV)'); xlabel('time (ms)'); 
title("ML neuron: Det="+det+", Omg="+n.system_size+", Burst="+n.burst+...
        ", I="+n.Istim+", Syn="+ n.synapse(1))
legend({'Neuron-1','Neuron-2'})

%% HCO Synchronization measure
W = 50; plt=2; dt=0.1;
ind_spk_1 = spikes{1,1};
ind_spk_2 = spikes{2,1};
[C] = GetBinlessCorr(ind_spk_1,ind_spk_2,W,t,dt,plt);
disp("GSCorr :::: C : " + C );