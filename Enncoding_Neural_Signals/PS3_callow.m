%% PS3
clear all
close all
load("LGN_FFdata.mat")

%% Part 1A

theta=[1:1:180];
ro=4
rmax=46;
sigma=30;
target_ori = 90;

for theta=1:180
    orituning(theta)=ro+rmax*exp(1).^(-(theta-target_ori).^2/(2*sigma^2));
end
figure
plot(orituning)
ylabel('Firing Rate (Hz)', 'Fontsize', 14)
xlabel('Orientation (Degrees)', 'Fontsize', 14)
title('Tuning Curve')

%% Part 1B


mean_spk_cnt = orituning*.2; % orituning is in Hz (spikes per second) and we are sampling from .2 seconds
stdev_spk_cnt = sqrt(mean_spk_cnt); %know variance equals mean and variance = standard deviation ^2
figure
plot(stdev_spk_cnt)
ylabel('Spike count standard deviation', 'Fontsize', 14)
xlabel('Orientation (Degrees)', 'Fontsize', 14)
title('Standard deviation of spike count')


%% Part 1C
figure
plot(mean_spk_cnt)
hold on
plot(mean_spk_cnt+stdev_spk_cnt)
hold on 
plot(mean_spk_cnt+-1*stdev_spk_cnt)
hold on
yline(50*.2)
ylabel('Average spike count', 'Fontsize', 14)
xlabel('Orientation (Degrees)', 'Fontsize', 14)
title('Spike count with standard deviation included')
legend('mean spike count', 'positive standard deviation', 'negative standard deviation')

% Angles between 65 and 115 degrees correspond to a firing rate of 50Hz.
%% Part 1D

% The mean spike count at 60 degrees is 6.38 vs 10 at 90 degrees
% orientation. It would seem that 90 degrees is better 


%% Part 2A
z = numel(FFstim);
stimtime = [];
for t=1:z
    r=t*DTstim;
    stimtime(t,1)=FFstim(t);
    stimtime(t,2)= r;
end
firingrate = histc(FFspks, -1e-8:DTstim:120.0081);

stim_spk=[]
for t=1:numel(firingrate)
    if firingrate(t)>0
        stim_spk=[stim_spk, FFstim(t,1)];
    end
    t;
end

mean(stim_spk)

%This spike triggered average is not surprising because the stimulus seems pretty close to 0, which
% indicates that the stimulus leading to the spikes are both negatiive
%and positive. This makes sense because the stimulus with no lag is
%probably not actually related to the firing of the spike and thus the
%stimulus would be somewhat random. This would result in an average close
%to 0 as the average of a large number of random stimuli would likely hover
%near 0.

%% Part 2B
TS=.100/.0083; %Number of timesteps for stimulus (About 12 timepoints)
mean_stim_spk_lag=[]
for lag=1:TS
stim_spk=[];
mean_stim_spk=[];
for t=10:numel(firingrate) % start after first 2 spikes
    if firingrate(t)>0
        stim_spk=[stim_spk, FFstim(t-lag,1)]; % add timing of spike stimulus   
    else
    t;
    end
end
mean_stim_spk=mean(stim_spk); % average mean spike stim for this version of lag
mean_stim_spk_lag =[mean_stim_spk_lag, mean(mean_stim_spk)];
end
 figure
 plot(flip(mean_stim_spk_lag))
 xlabel("lag (.0083ms steps)")
 ylabel("Stimulus")
 
 %% Part 3A

load("FFrf.mat");
Kdan1=(Kdan)/norm(Kdan);
g = g_convolve(FFstim,Kdan1,1);
figure
plot(g)
xlim([800/8.3 1000/8.3])

%% Part 3B
x=histc(g,-6:0.2:6);
for t = 1:61
    v=-6+(t*.2);
    x(t,2)=v;
end
figure
plot(x(:,2),x(:,1))
%% Part 3C


f=[]
for t=1:14391 % Create an array for filtered stimulus just when the spikes occur
    if firingrate(t)>0;
        f=[f, g(t)];
    else
    t;
    end
end

f=transpose(f);



y=histc(f,-6:0.2:6);
for t = 1:61
    v=-6+(t*.2);
    y(t,2)=v;
end

figure
plot(x(:,2),x(:,1))
hold on
plot(y(:,2),y(:,1))

%% Part 3D
% Set g and f to just when above zero, if not, set to 0
for t=1:numel(g)
    if g(t)<0
        g(t)=0;
    else
        t;
    end
end

for t=1:numel(f)
    if f(t)<0
        f(t)=0;
    else
        t;
    end
end

h=histc(f,0:.2:4);
j=histc(g,0:.2:4);

spk_nonlinearity = h./j;
for t = 1:21
    v=0+(t*.2);
    spk_nonlinearity(t,2)=v;
end
figure
plot(spk_nonlinearity(:,2), spk_nonlinearity(:,1))
xlabel("stimulus")
ylabel("Probability")
xlim([0 5])
ylim([0 1.1])



%% Part 4A

%Was not able to complete part 4

firingrate = histc(FFspksR, -1e-4:DTstim:10);
firingrate = transpose(firingrate)

figure
plot(firingrate)
xlim([800 1000])
xlabel("time (ms)")
ylabel("Firing rate (Hz)")

%% Part 4B

g = g_convolve(FFstimR,Kdan1,1);
var(g)
plot(g)
xlim([800 1000])

%% Part 4C
filt_stim = histc(g,0:.2:4);
pred_firingrate = spk_nonlinearity(:,1).*filt_stim;


figure
plot(firingrate)
xlim([800 1000])
figure
plot(spk_nonlinearity(:,2),pred_firingrate,':')
xlim([.8 1])
xlabel("time (seconds)")
ylabel("Firing rate (Hz)")
hold off
