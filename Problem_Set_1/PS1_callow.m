%% Cell 1 clear all and close out
clear all
close all
rng(5)
%% Main variabless defined
bin = .001; %bin size for spike generation
f = 20; %firing rate in Hz
trial_length = 1; % 1 second trial length
num_bins = ceil(trial_length/bin); %num bins
p = f*bin;     % assumed probability of a spike
%% Part 1a Create 100 trials of spikes
trials = []; %Create empty array for the 100 trials
spks0 = []; %Create empty array for each individual trial
spks_per_trial = []; %create empty array to put # spikes per trial
spks100 = []; %Create empty array to store spikes occuring in firt 100 ms ***for part 1d***
for t = 1:100  
    for x = 1:num_bins; %for loop to append spikes into the 1000 bins
     spks_rand(x) = rand(1); %create random spikes for each bin
     if spks_rand(x) < p;  %if spikes less than .02 then append current index to spike
       spks0 = [spks0 x];
     end
    end
    trials = spks0; % append the leftover array with the bins that each spike occurs in to trials array
    spks2(t) = [numel(trials)]; % adds number of spikes in each trial **for part 1d**
    spks100 = spks0;
    spks100(spks100>100)=[]; %create a new array that removes all spikes that don't occur in first 100ms
    spks3(t) = [numel(spks100)];
    spks0 = [spks0 -1]; % add a -1 to indicate end of an individual trial for raster
end 

%% Part 1b Display first 20 repeats with raster.m

raster(spks0(1,1:435), [0 1000]) % Display first 20 spike repeats with raster ****MUST RUN raster.m SCRIPT FIRST***

%% Part 1c Determine firing rate and create histogram

spks_ms = spks0/1000; %Convert spike bins into ms units
firingrate = histc(spks_ms, -1e-8:0.01:1)/100; % create firing rate (with probability values) using histc bins of .01 (so probability of firing every 10 ms)
spks_sec =firingrate*100; % Multiply the probability times 100 to get firing rate in spikes/sec (Hz)
figure % graph firing rate
bar(0:10:1000,spks_sec,'histc')
xlabel("Trial Period")
ylabel("Spikes/sec (Hz)")

%% Part 1d Find average number of spikes in each trial and put in a histogram

spks_per_trial = [];
spks_per_trial = spks2(1); %add first number of spikes to the array spks per trial
for t = 1:99
    x = spks2(t+1) - spks2(t)-1; % set x equal to the difference between the total spikes and previous number of spikes must subtract 1 because after first array there is a negative one for each trial
    spks_per_trial = [spks_per_trial x]; %append the difference between total spikes and the previous total to get number of spikes for that trial
end

figure
bar(1:100, spks_per_trial)
histogram(spks_per_trial, 0:40)
xlabel("Number of Spikes")
ylabel("Counts")
%% Part 1e mean number of spikes accross all trials
average_number_spikes = mean(spks_per_trial); %calculate mean number of spikes across all trials in spks_per_trial
spike_variance = var(spks_per_trial); %calculate variance across all trials in spks_per_trial

fano_factor = spike_variance/average_number_spikes; %fano factor =variance/mean and a fano factor of one indicates neuron may be using rate coding

%% Part 1f recteate histogram and recalculate mean, variance and fanofactor, based on only the firing rate in first 100 ms.
spks_per_trial_100 = [];
spks_per_trial_100 = spks3(1); %add first number of spikes to the array spks per trial
for t = 1:99
    x = spks3(t+1) - spks3(t)-1; % set x equal to the difference between the total spikes and previous number of spikes must subtract 1 because after first array there is a negative one for each trial
    spks_per_trial_100 = [spks_per_trial_100 x]; %append the difference between total spikes and the previous total to get number of spikes for that trial
end

figure
bar(1:100, spks_per_trial_100)
histogram(spks_per_trial_100, 0:10)
xlabel("Number of spikes")
ylabel("Counts")

mean_100 = mean(spks_per_trial_100)
variance_100 = var(spks_per_trial_100)
fano_factor_100 = variance_100/mean_100

%% Part 1G Do answers match what you would expect.
%Yes! Considering this modeling of spike firing rates is based off of a
%poison neuron we would expect to see a fano factor of approximately 1. The
%fano_factor over the entire second for the 100 trials is approximately .97
%while the fano factor for the first 100ms of all 100 trials is
%approximately .91. These are relatively close to 1, with the number
%approaching 1 as we sample from a longer time period.

%% Part 2a How many repeat in the data, mean firit rate and show a raster plot for first 200 ms of data for all repeats
load('LGN_FFdata.mat')

num_repeats = sum(FFspksR == -1) %Sum number of -1 which will give us the number of repeats.
mean_firing_rate = mean(FFspksR>0) %Mean of FFspksR excluding -1 values.
raster(FFspksR(FFspksR<.2), [0 .2])

%% Part 2b
FFspksR_no_neg =FFspksR(FFspksR>0);% Create FFspksR with no negative 1
for i = 1:numel(FFspksR_no_neg)
    FF_ISI(i) = FFspksR(i+1)-FFspksR(i); % look at the difference in time between each spike
end
FF_ISI_50 =FF_ISI(FF_ISI<.05); % create an array of difference in time between spikes that only includes ISI of less than .05
figure
nbins=100;
edges = [0:.001:.05];
histogram(FF_ISI_50, edges) % plot a histogram of interspike intervals that occur within 50 ms.
xlabel("Inter Spike Interval")
ylabel("Occurances")

% The results surprised me. I would have expected either a more normal
% distribution of the interspike interval that would be consistent with a
% poisson neuron.I was also surprised by the number of ISI's that were so
% small <.005. I would have thought this would be too short of an ISI based
% on what we discussed in lecture.

%% Part 2c
FFspksR_160_180 = FFspksR(FFspksR>.160 & FFspksR<.180 | FFspksR==-1); %Find all spikes that occur from 160 ms to 180 ms. Leave -1's


loc_neg = find(FFspksR_160_180==-1); %find location of all -1's

for w = 2:numel(loc_neg)
    spks_per_trial_160_180(1) = [loc_neg(1)-1]; %Add the first trials spike number which is the first -1 (which occurs at 1) to get a total of 0 spikes in first trial (can't use the equation below beccause can't have a 0 index
    spks_per_trial_160_180(w) = [loc_neg(w)-loc_neg(w-1)]; %Subtract the next -1 location from the previous -1 location to get the number of spikes for each trial.
end

var(spks_per_trial_160_180)
mean(spks_per_trial_160_180)
Fano_Factor_LGN = var(spks_per_trial_160_180)/mean(spks_per_trial_160_180)

%No I would not have expected to see such a high number of spikes nor suuch
%a low variance in the number of spikes in this section. The mean number of
%spikes would indicate 2.6 spikes per 20ms period which seems very high
%(more than 100Hz) comparedd to the neurons we have previously discussed
%that are poisson neurons. The variance is also very low. This must mean
%that this neuron is similar to the ones we have discussed in class that
%respond very consistently to naturally stimuli. We would expect a
%Fano Factor of about 1 if it was a poisson neuron however the Fano Factor
%was very low for this neuron.