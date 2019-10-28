%% Part 1A
clear all
close all
T = .5; % time interval in seconds of the experiment
dt = 1*10^-5;% time step in .1 ms
NT = ceil(T/dt); % Number of steps
I0 = 1; %I0template set to 1 at 100ms for 300ms and then go back to 0
Iext = zeros(NT,1); %Injected current starts at all zeros

t1=0.1;t2 = 0.4; % Create timepoints at 100 and 400 ms

nt1 = ceil(t1/dt); nt2 = ceil(t2/dt); % create array with step location of t1 and t2
Iext(nt1:nt2) = I0; %Create time point in which we input the current from 100 to 400ms at 1 amp
figure
plot(Iext)
ylim([-.1,1.1])
xlabel("time steps (.01ms increments)")
ylabel("Amps")

Iext = 7*10^-9*Iext; %Set new current injection to 7nA

figure
plot(Iext)
ylim([-1*10^-8, 1*10^-8])
xlabel("time steps (.01ms increments)")
ylabel("Amps")
%% Part 1B
Vt = zeros(NT,1);
tau =.01; % Set Tau to 10 ms
R = 2.5*10^6; % Set the resistance to 2.5*10^-6 Ohms
Vthresh = 15*10^-3; % Set threshold to spike
Vspike = 60*10^-3; % Set value if spikes
Vreset = -5*10^-3; % Set reset value if spikes
for t=1:NT-1
 if Vt(t)==Vspike % If current voltage = spike voltage than make the next value a spike reset
     Vt(t+1)=Vreset;
 elseif Vt(t)>Vthresh % If current voltage is greater than Vthresh than set next value to spike value
     Vt(t+1)=Vspike;
 else % If neither of the above occur than calculate a new dvdt based on the voltage and add to get next voltage.
  dVdt = (-Vt(t)+R*Iext(t))/tau;   
  Vt(t+1) = Vt(t)+dVdt*dt;  
     
 end
% if Vt(t+1) > Vthresh
%     Vt(t+1) = Vspike
%     Vt(t+2) = Vreset
end

plot(Vt)
ylim([-1*10^-2, 6*10^-2])
xlabel("time steps .01ms increments")
ylabel("Voltage (Volts)")


%% Part 1C
Vi=zeros(40,1); % create a new matrix with 40 zeros
for I = 1:1:40
    Iext(nt1:nt2)=(I-1)*10^-9*.5; % Set the current to nanometers and steps of .5 starting at 0
    for t=1:NT-1;
        if Vt(t)==Vspike; 
           Vt(t+1)=Vreset;
        elseif Vt(t)>Vthresh;
           Vt(t+1)=Vspike;
        else
           dVdt = (-Vt(t)+R*Iext(t))/tau;   
           Vt(t+1) = Vt(t)+dVdt*dt;  
        end
    end
    current = (I-1)*10^-9*.5; % Create a value for applied current during the iteration
    Vi(I,1)= current; % Store in column 1
    cnt = simple_spike_count(Vt, 15*10^-3); % Create a value for number of spikes based on that current injection
    firingrate = cnt/NT; % Create firing rate value
    Vi(I,2)= firingrate; % Add firing rates to column 2 

end
plot(Vi(:,1),Vi(:,2))
xlabel("Current Injected (Amps)")
ylabel("Firing rate")
%% Part 2A
ts=(1:NT)*dt;
Iext(nt1:nt2)=7;
[V,n,m,h] = HHmodel(ts,Iext);
figure
plot(V)
xlabel("time")
ylabel("voltage")
%% Part 2B
figure
subplot(4,1,1)
plot(V)
xlim([10000 13500])
subplot(4,1,2)
plot(m)
xlim([10000 13500])
subplot(4,1,3)
plot(n)
xlim([10000 13500])
subplot(4,1,4)
plot(h)
xlim([10000 13500])

%m increases first and runs pretty in line with the increase in voltage of
%the first spike. This makes sense because as the number of sodium channels
%opening increases we will see a quick depolarization in the cell. n
%increases soon after m as the K+ channels start to open closely after the spike
%starts. This would help repolarize the cell and thus end the the
%current spike. However, n occurs slightly later and movement of K+ is
%slightly slower than Na+ and thus the repolarization is timed so that it starts after the spike hase peaked. Finally Na+ inactivation occurs in a
%similar time frame as K+ opening. Both of these factors help repolarize
%similar time frame as the K+. These two factors contribute to the quick
%repolarization of a neuron.
%% Part 2C
FI=zeros(20,1);
for I = 1:20
    Iext(nt1:nt2)=I; % Create currents from 1-20nA
    [V,n,m,h] = HHmodel(ts,Iext);% Run HHmodel for different currents
    cnt = simple_spike_count(V, 0); % Create a value for number of spikes based on that current injection
    firingrate = cnt/NT; % Create firing rate value
    FI(I,1)=I; % Set first column of FI to current
    FI(I,2)= firingrate; % Add firing rates to column 2 
end
figure
plot(FI(:,1),FI(:,2))
xlabel("Current in nA")
ylabel("Firing Rate")

%This curve is less linear in its increase for the IF than the IF curve in
%question 1. This difference would be a result of the different dynamics
%and parameters created in the two models we used.
%% Part 3A
Inoise = lowpass_gwn( 1000, 0.5, 1e-5 );
Inoise = Inoise*3;

ts=(1:NT)*dt;
Iext(nt1:nt2)=7;
[V1,n,m,h] = HHmodel(ts,Iext);
figure
plot(V1)
xlabel("time")
ylabel("voltage")
hold on
Iext(nt1:nt2)=Iext(nt1:nt2)+Inoise(nt1:nt2);
[V2,n,m,h] = HHmodel(ts,Iext);
plot(V2)
%% Part 3B

FI=zeros(32,1);
for t=1:32
  Inoise = lowpass_gwn(1000, 0.5, 1*10^-5);
  Inoise = Inoise*3;
  Iext(nt1:nt2)=7;
  Iext(nt1:nt2)=Iext(nt1:nt2)+Inoise(nt1:nt2);
  [V3,n,m,h] = HHmodel(ts,Iext); 
  cnt = simple_spike_count(V3, 0); % Create a value for number of spikes based on that current injection
  firingrate = cnt/NT; % Create firing rate value
  FI(t,1)=t; % Set first column of FI to current
  FI(t,2)= firingrate; % Add firing rates to column 2   
end
figure
plot(FI(:,1),FI(:,2))
xlabel("Trial number")
ylabel("Firing rate")
%% Part 3C
mean(FI(:,2))
var(FI(:,2))
FanoFactor=var(FI(:,2))/mean(FI(:,2))

%No, the fano factor is very far from 1 and therefore it is unlikely that
%this is a poisson neuron. The HHmodel should cause a spike with an uptick
%in in current and the current fluctuation should be relatively consistent.
%Meaning the spike train should be relatively deterministic and have low
%variancce.
