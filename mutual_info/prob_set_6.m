%% PROBLEM SET #6 SOLUTIONS
clear all
close all
%% PROBLEM 1A. Single-bin information capacity
dt = 0.001; % sec
r = 40; % Hz

% Calculate information capacity (single-bin response entropy)
HR = myentropy([r*dt 1-r*dt])
% = 0.2423 bits

% Convert to rate
HR/dt
% = 242.3 bits/s

%% 1B. Optimal entropy rate
hR = zeros(1,1001);
for r = 0:1000
  hR(r+1) = myentropy([r*dt 1-r*dt])/dt;
end

figure
plot((0:1000),hR,'LineWidth',1)
box off
axis([0 1000 0 1005])
xlabel( 'Firing rate (Hz)')
ylabel( 'Information capacity rate (bits/s)')

% The optimal firing rate is at 500 Hz, because this is were the neuron has
% an equal probability of firing a spike or silence, giving the maximum
% response entropy (1 bit) possible with only two responses.

%% 1C. Time resolution
dt = 0.0005; % sec
r = 40; % Hz

hR = myentropy([r*dt 1-r*dt])/dt
% = 0.282.9 bits/s

% Yes: this is a big increase, and does concern me, because it does not
% depend on the neuron, but rather the time resolution of my analysis.
% In other words, the overall information capacity does depend on how
% finely I analyze the data, but if I'm doing things right, these factors
% will cancel out when analyzing the information (and not just entropy).

%% 1D. Information rate

% The information rate is zero, because the firing rate is constant across
% the experiment, and as a result H[R|s] = H[R]. Thus, the noise entropy is
% equal to the total response entropy, and no information is conveyed about
% the stimulus.

%% PROBLEM 2
% Code from last problem set (setup)
rmax = 40;
r0 = 10;
sigma = 30;
theta0 = 180;
thetas = 2:2:360;
rs = 0:1:80;

% Tuning curve
f = (r0 + rmax*exp(-(thetas-theta0).^2/(2*sigma^2)));
jpd = zeros(length(thetas),length(rs));

% Calulate prior
ps = 1/length(thetas);   % const prob of stimulus

for a = 1:length(thetas)
  prs = exp(-(rs-f(a)).^2/(2*f(a)));
  prs = prs/sum(prs);  % normalize
  jpd(a,:) = ps*prs; 
end

% Check normalization
sum(jpd(:))

%% 2A. Entropy of prior
HS = myentropy(sum(jpd'))

% Also can calulate using ps = 1/180 is uniform in 180 bins
HS = 180 * -(1/180) * log2(1/180)
% = 7.4919 bits

%% 2B. Specific information of responses
pr = sum(jpd);
Isp = zeros(1,length(rs));

for n = 1:length(rs)
%  r = rs(n);   % dont need to know this
  psr = jpd(:,n) / pr(n);
  Isp(n) = HS - myentropy( psr );
end

figure
plot(rs,Isp)
xlabel('Firing Rate (Hz)');  ylabel( 'I_{sp} (bits)')

%% 2C. Mutual information
% Calculate mutual information using Isp
I = sum( pr .* Isp )
% I = 1.1307

%% 2D. Stimulus-specific information
Iss = zeros(1,180);

for n = 1:180
  prs = jpd(n,:) / sum(jpd(n,:));
  Iss(n) = sum( prs .* Isp );
end

figure
plot(thetas,Iss)
xlabel('Angle');  ylabel( 'I_{SSI} (bits)')
axis([0 360 0.5 3])

%% 2E. Maximum of the SSI for [high] noise.

% The SSI is maximum at the peak of the tuning curve.  This is because
% the respones at the peak are most informative (from Isp).

%% 2F. Recalculate the specific information for rmax = 80 Hz higher firing rates

% Recalculate new tuning curve
rmax = 80;
rs = 0:150;

f = (r0 + rmax*exp(-(thetas-theta0).^2/(2*sigma^2)));
jpd = zeros(length(thetas),length(rs));

% Recalculate new JPD
for a = 1:length(thetas)
  prs = exp(-(rs-f(a)).^2/(2*f(a)));
  prs = prs/sum(prs);  % normalize
  jpd(a,:) = ps*prs; 
end

% Check normalization
sum(jpd(:))

% Recalculate new marginal distribution and specific informations
pr = sum(jpd);
Isp = zeros(1,length(rs));

for n = 1:length(rs)
  r = rs(n);   % dont need to know this
  psr = jpd(:,n) / pr(n);
  Isp(n) = HS - myentropy( psr );
end

figure
plot(rs,Isp)
xlabel('Firing Rate (Hz)');  ylabel( 'I_{sp} (bits)')

%% 2G. Recalculate the stimulus specific information 
Iss = zeros(1,180);

for n = 1:180
  prs = jpd(n,:) / sum(jpd(n,:));
  Iss(n) = sum( prs .* Isp );
end

figure
plot(thetas,Iss)
xlabel('Angle');  ylabel( 'I_{SSI} (bits)')
axis([0 360 0.5 3.5])

%% 2H. Why has the peak changed?

% Raising the overall firing rate means that there is effectively less 
% noise in the response, because the variability in spike count goes as the
% square root of the standard deviation (under the Poisson assumption). 
% This means that there is more specific information in responses on the
% flank of the tuning curve, and correspondingly more SSI there.
