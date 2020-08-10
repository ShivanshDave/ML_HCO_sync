function [V,t,spikes,X] = ML_network(demo, det, n, seed)
%%    (for Demo) : [V,t,spikes] = ML_network(1,0);  % Stochastic demo
%     (for Demo) : [V,t,spikes] = ML_network(1,1);  % Deterministic demo
%     (from fxn) : [V,t,spikes] = ML_network(2,det,n);
%     (other..)  : [V,t,spikes] = ML_network(0);  % (Uses network @ line-27)
% 
% Actual ML neuron is the function "ML_neurons" (at the very end; line# ~150)

% Read Network Design
if demo == 1
    % Demo: 1,2 : HCO-excitatory w/ stim current
    %       3,4 : HCO-inhibitory w/ stim current
    %       5   : Un-connected w/ low stim current
    %       6   : Busrting neuron (w/ NO stim-current) w/ NO synapse
    %       8,7 : HCO-excitatory w/ Busrting neuron
    %      10,9 : HCO-inhibitory w/ Busrting neuron
    % Neuron #  1  2  3  4  5  6  7  8  9  10
    net     = [ 2, 1, 4, 3, 0, 0, 8, 7,10, 9]; % Neuron connection
    synapse = [ 1, 1,-1,-1, 0, 0, 1, 1,-1,-1]; % Synapse : inhi/exci/off
    burst   = [ 0, 0, 0, 0, 0, 1, 1, 1, 1, 1]; % Burst mode : on/off
    if det==1
        Istim=[90,90,90,90,88, 0, 0, 0, 0, 0]; % Current-clamp (nA)
    else
        Istim=[85,85,85,85,80, 0, 0, 0, 0, 0]; % Current-clamp (nA)
    end
    dt = 0.1;  tmax = 3e3; system_size = 1e3; seed = 10;
    
elseif demo == 0
    % Network Design
    net     = []; % Neuron connection
    synapse = []; % Synapse : inhi/exci/off
    Istim   = []; % Current-clamp (nA)
    burst   = []; % Burst mode : on/off
    system_size = []; % noise (array or single value)
    tmax = 3e3;
    dt = 0.1;
else 
    % read data from 'n' variable
    net=n.net; synapse=n.synapse; Istim=n.Istim; burst=n.burst;  
    system_size=n.system_size; tmax=n.tmax; dt=n.dt; demo=n.demo;
end

% Simulation
% Initialization
if exist('seed','var'), rng(seed); end
sz = size(net);
omega = system_size.*ones(sz);
Vm   = -100*rand(sz); % Start with random voltage
nK   = rand(sz); % Start with random fraction
nCa  = rand(sz); % Start with random fraction
Iint = 0*ones(sz);
t = (0:dt:tmax)';
synapse(net==0) = 0; % turn-off synapse for unconnected neurons
net(net==0) = find(net==0); % (note-1)
V = nan(length(t),sz(2));
M = nan(length(t),sz(2));
N = nan(length(t),sz(2));

% Find membrane voltage at each time-step
for i=1:length(t)
    [Vm,nK,nCa,Iint] = ML_neurons(det, dt, Istim, Iint, omega, ...
        nK, nCa, Vm, Vm(net), synapse, burst);
    V(i,:) = Vm; M(i,:) = nCa; N(i,:) = nK; 
end
% Save Ion-channels state
X.M = M; X.N = N; X.Omega = system_size; X.t = t;

% Identify Spikes
spikeThr = 15;  % mV; Neural spike detection threshold
refTh = 15/dt; % ms->Steps; Refractory period for spikes
spikes=cell(sz(2),3);

for i=1:sz(2)
    % Find Peaks
    xval=V(:,i);
    len_x = length(xval); val=[]; index=[]; xi=2;  % start at second data point
    while xi < len_x-1
        if xval(xi) > xval(xi-1)
            if xval(xi) > xval(xi+1)                    % definite max
                val =[val xval(xi)];
                index = [ index xi];
            elseif xval(xi)==xval(xi+1) && xval(xi)==xval(xi+2)	% 'long' flat spot
                xi = xi + 2;                      % skip 2 points
            elseif xval(xi)==xval(xi+1)                 % 'short' flat spot
                xi = xi + 1;                      % skip one point
            end
        end
        xi = xi + 1;
    end
    % Filter Peaks
    index(val<=spikeThr) = [];          % Apply Vm threshold
    val(val<=spikeThr) = [];
    ind=2; 
    while ind <= length(index)          % Apply refractory period      
       if abs(index(ind)-index(ind-1)) < refTh      
           index(ind)=[];
           val(ind)=[];
       else
           ind=ind+1;
       end
    end
    % Save Spikes
    spikes{i,1} = index;                % Save Spike index
    spikes{i,2} = val;                  % Save Spike amplitude (mV)
    spikes{i,3} = t(index);             % Save Spike time (ms)
end

% Plot traces
if demo == 1
        % plot Demo network
    figure
    suptitle("Plain HCO Demo (Omg="+system_size+", I="+Istim(1)+")")
    for i = 1:5
        subplot(5,1,i)
        plot(t,V(:,i))
        hold on
        scatter(spikes{i,3}, 1.4*spikes{i,2}, 55, 'v', 'filled')
        if i==1, title('HCO-excitatory w/ stim current');
        elseif i==3, title('HCO-inhibitory w/ stim current');
        elseif i==5, title('Un-connected w/ low stim current');
        end
    end
    ylabel('Vm (mV)')
    xlabel('time (ms)')

    figure
    suptitle("Bursting HCO Demo (Omg="+system_size+", I="+Istim(1)+")")
    for i = 1:5
        subplot(5,1,i)
        plot(t,V(:,i+5))
        hold on
        scatter(spikes{i+5,3}, 1.4*spikes{i+5,2}, 55, 'v', 'filled')
        if i==2, title('HCO-excitatory w/ internal current');
        elseif i==4, title('HCO-inhibitory w/ internal current');
        elseif i==1, title('Un-connected neuron');
        end
    end
    ylabel('Vm (mV)')
    xlabel('time (ms)')
    
elseif demo == 0 
    % all neurons in one plot
    figure
    suptitle('Vm (mV) vs  time for each neuron')
    for i = 1:sz(2)
        subplot(sz(2),1,i)
        plot(t,V(:,i))
        hold on
        scatter(spikes{i,3}, 1.4*spikes{i,2}, 55, 'v', 'filled')
    end
    xlabel('time (ms)')
end

% NOTES:
% 
% (1) un-connected neurons get a dummy connection to Neuron1, for valid indexing.
% I have made sure that this dummy connection does not affect in anyway,
% because the synapse is forcefully "turned-off" in the previous line.
end
