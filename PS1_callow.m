%% Cell 1 clear all and close out
clear all
close all
rng(5)
%% Main variabless defined
bin = .001; %bin size
f = 20; %firing rate in Hz
trial_length = 1 % 1 second trial length
num_bins = ceil(trial_length/bin) %num bins
p = f*bin     % assumed probability of a spike
%% Part 1a Create 100 trials of spikes
trials = [] %Create empty array for the 100 trials
spks0 = [] %Create empty array for each individual trial
spks_per_trial = []
for t = 1:100  
    for x = 1:num_bins %for loop to append spikes into the 1000 bins
     spks_rand(x) = rand(1); %create random spikes for each bin
     if spks_rand(x) < p  %if spikes less than .02 then append current index to spike
       spks0 = [spks0 x];
     end
    end
    trials = spks0; % append the leftover array with the bins that each spike occurs in to trials array
    x = numel(spks0);
    %spks_per_trial = [spks0 numel(spks0(t))]
    spks0 = [spks0 -1]; % add a -1 to indicate end of an individual trial for raster
end 

%% Part 1b Display first 20 repeats with raster.m

raster(spks0(1,1:435), [0 1000]) % Display first 20 spike repeats with raster ****MUST RUN raster.m SCRIPT FIRST***

%% Part 1c Determine firing rate and create histogram

spks_seconds = spks0/1000; %Convert spike bins into ms units
firingrate = histc(spks_seconds, -1e-8:0.01:1)/100; % create firing rate (with probability values) using histc bins of .01
firingrate =firingrate*100; % Multiply the probability times 100 trials to get firing rate in Hz
figure % graph firing rate
bar(0:10:1000,firingrate,'histc')
xlabel("Trial Period")
ylabel("Spikes/sec (Hz)")

%% Part 1d Find average number of spikes in each trial and put in a histogram
spikes = [] %Create empty array for the 100 trials
spks1 = [] %Create empty array for each individual trial
spks_per_trial = []
for t = 1:100  
    for x = 1:num_bins %for loop to append spikes into the 1000 bins
     spks_rand(x) = rand(1); %create random spikes for each bin
     if spks_rand(x) < p  %if spikes less than .02 then append current index to spike
       spks1 = [spks1 x];
     end
    end
    spks2(t)=[numel(spks1)]; %create array spks 2 that has the number of firing with no -1 
end
spks_per_trial = spks2(1,1); %add first number of spikes to the array spks per trial
for t = 1:99
    x = spks2(t+1) - spks2(t); %set x equal to the difference between the total spikes and previous number of spikes
spks_per_trial = [spks_per_trial x]; %append the difference between total spikes and the previous total to get number of spikes for that trial
end

figure
bar(1:100, spks_per_trial)
xlabel("trials")
ylabel("Number of Spikes Per Trial")
%% Part 1e mean number of spikes accross all trials
average_number_spikes = mean(spks_per_trial) %
trial_spike_variance = 
fano_factor = 

%% Part 2a How many repeat in the data, mean firit rate and show a raster plot for first 200 ms of data for all repeats


