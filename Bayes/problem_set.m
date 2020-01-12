%% Problem Set 5
clear all
close all
%% Par1 A Prior
x=-1:0.01:4;
mu = 1
sigma = .5
f=exp(-(x-mu).^2/(2*sigma^2));
f=f./sum(f,'all');
figure
subplot(3,1,1)
plot(f)
xlabel("feedback")
ylabel("probability")

%% Part1 B Visual Uncertainty

x=-1:0.01:4;
mu = 2
sigma1 = .1
sigma2 = .2
sigma3 = .5
f1=exp(-(x-mu).^2/(2*sigma1^2));
f1=f1./sum(f1,'all');
f2=exp(-(x-mu).^2/(2*sigma2^2));
f2=f2./sum(f2,'all');
f3=exp(-(x-mu).^2/(2*sigma3^2));
f3=f3./sum(f3,'all');
subplot(3,1,2)
plot(f1)
hold on
plot(f2)
hold on 
plot(f3)
xlabel("feedback")
ylabel("probability")
legend('.1 cm', '.2 cm', '.5 cm')

%% Part1 C Posterior

f_1=f.*f1;
f_1=f_1./sum(f_1,'all');
f_2=f.*f2;
f_2=f_2./sum(f_2,'all');
f_3=f.*f3;
f_3=f_3./sum(f_3,'all');
subplot(3,1,3)
plot(f)
hold on
plot(f_1)
hold on
plot(f_2)
hold on
plot(f_3)
xlabel('feedback')
ylabel('probability')
legend('no feedback','.1 cm', '.2 cm', '.5 cm')

%% Error Calculation

% Mean Squared Error using visual feedback and prior experience
x=-1:0.01:4;
MSE1=[]
for x=1:numel(x)
    MSE1(x) = f_1(x)*(x-2)^2;
end
MSE1_error = sqrt(sum(MSE1,'all'))

x=-1:0.01:4;
MSE2=[]
for x=1:numel(x)
    MSE2(x) = f_2(x)*(x-2)^2;
end
MSE2_error = sqrt(sum(MSE2,'all'))

x=-1:0.01:4;
MSE3=[]
for x=1:numel(x)
    MSE3(x) = f_3(x)*(x-2)^2;
end
MSE3_error = sqrt(sum(MSE3,'all'))


% Mean Squared Error just using visual feedback
x=-1:0.01:4;
MSE1_v=[] 
for x=1:numel(x)
    MSE1_v(x) = f1(x)*(x-2)^2;
end
MSE1_v_error = sqrt(sum(MSE1_v,'all'))

x=-1:0.01:4;
MSE2_v=[]
for x=1:numel(x)
    MSE2_v(x) = f2(x)*(x-2)^2;
end
MSE2_v_error = sqrt(sum(MSE2_v,'all'))

x=-1:0.01:4;
MSE3_v=[]
for x=1:numel(x)
    MSE3_v(x) = f3(x)*(x-2)^2;
end
MSE3_v_error = sqrt(sum(MSE3_v,'all'))
% The mean squared errors are larger when just using the visual feedback
% without priors

%% Advantage of Bayesian Estimation
x=-1:0.01:4;
mu = 1
sigma1 = .1
sigma2 = .2
sigma3 = .5
f1=exp(-(x-mu).^2/(2*sigma1^2));
f1=f1./sum(f1,'all');
f2=exp(-(x-mu).^2/(2*sigma2^2));
f2=f2./sum(f2,'all');
f3=exp(-(x-mu).^2/(2*sigma3^2));
f3=f3./sum(f3,'all');
figure 
subplot(3,1,1)
plot(f)
subplot(3,1,2)
plot(f1)
hold on
plot(f2)
hold on 
plot(f3)
xlabel("feedback")
ylabel("probability")
legend('.1 cm', '.2 cm', '.5 cm')

% Posterior

f_1=f.*f1;
f_1=f_1./sum(f_1,'all');
f_2=f.*f2;
f_2=f_2./sum(f_2,'all');
f_3=f.*f3;
f_3=f_3./sum(f_3,'all');
subplot(3,1,3)
plot(f)
hold on
plot(f_1)
hold on
plot(f_2)
hold on
plot(f_3)
xlabel('feedback')
ylabel('probability')
legend('no feedback','.1 cm', '.2 cm', '.5 cm')

% Error Calculation

% Mean Squared Error using visual feedback and prior experience
x=-1:0.01:4;
MSE1=[]
for x=1:numel(x)
    MSE1(x) = f_1(x)*(x-2)^2;
end
MSE1_error = sqrt(sum(MSE1,'all'))

x=-1:0.01:4;
MSE2=[]
for x=1:numel(x)
    MSE2(x) = f_2(x)*(x-2)^2;
end
MSE2_error = sqrt(sum(MSE2,'all'))

x=-1:0.01:4;
MSE3=[]
for x=1:numel(x)
    MSE3(x) = f_3(x)*(x-2)^2;
end
MSE3_error = sqrt(sum(MSE3,'all'))

% Mean Squared Error just using visual feedback
x=-1:0.01:4;
MSE1_v=[] 
for x=1:numel(x)
    MSE1_v(x) = f1(x)*(x-2)^2;
end
MSE1_v_error = sqrt(sum(MSE1_v,'all'))

x=-1:0.01:4;
MSE2_v=[]
for x=1:numel(x)
    MSE2_v(x) = f2(x)*(x-2)^2;
end
MSE2_v_error = sqrt(sum(MSE2_v,'all'))

x=-1:0.01:4;
MSE3_v=[]
for x=1:numel(x)
    MSE3_v(x) = f3(x)*(x-2)^2;
end
MSE3_v_error = sqrt(sum(MSE3_v,'all'))
% The results show that with the new position of 1 cm for the visual
% feedback the mean squared error is significantly lower than when the mean
% visual feedback was 2 cm.

%% Part F

%% Part 2 Entropy Functions

%% Part 2A
%function out = mlog2(x)
%out = log2(x);
%end
mlog2(.5)

%% Part 2B

%function final=myentropy(x)
%out=[];
%for i=1:numel(x)
%    out(i)=x(i)*mlog2(x(i))
%end
%final=-sum(out);
%end
p = [0.25 0.25 0.25 0.25];
myentropy(p)

%% Part 3A

jpd=zeros(2,5);
jpd(2,:)=1
jpd(:,1)=.5
jpd=jpd./sum(jpd,'all');

ps = [.2 .2 .2 .2 .2];
% Specific information given a response of 1
prob1 = jpd(1,:)./sum(jpd(1,:)); % Probability of stimulu
H_SR_1 = myentropy(prob1);
H_S_1 = myentropy(ps);
specific_info_1 = H_S_1-H_SR_1

% Specific information given a response of 0
prob_0 = jpd(2,:)./sum(jpd(2,:));
H_SR_0 = myentropy(prob_0);
H_S_0 = myentropy(ps);
specific_info_0 = H_S_0-H_SR_0

%% Part B
% Mutual information from specific information
MI1 = sum(.1*specific_info_1+.9*specific_info_0);


% Mutual information
I=[]
pr =[.1 .9];
ps = [.2 .2 .2 .2 .2];
I=zeros(2,5);
for r=1:2
    for s = 1:5
        I(r,s)=jpd(r,s)*mlog2((jpd(r,s))./(pr(r)*ps(s)));
    end
end
MI2 = sum(I,'all');

%% Part 3C
% Calculate stimulus specific information
SSI = []

for s=1:5
SSI(s)=sum(jpd(1,s))*specific_info_1+sum(jpd(2,s))*specific_info_0;
end

SSI
sum(SSI)
