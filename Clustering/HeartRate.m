clear all
close all
clc

% loading the HeartData.mat
rng(1)
load('HeartData.mat')
subject_10 = HeartData(:, 10);

% plotting the histogram in normalized PDF of subject number 10 
figure(1)
histogram(subject_10, 30, 'Normalization', 'pdf')

gm_10 = fitgmdist(subject_10,2,'Options',statset('MaxIter',1000));
x = linspace(min(subject_10),max(subject_10),500);
x = x';
y = pdf(gm_10, x); 

% fitting the PDF curve of the subject number 10 over the hitogram plot
hold on
plot(x, y, 'r');
title('Probability Density of Heart Rate (bpm) of Subject No.10')
ylabel('Probability Density')
xlabel('Heart rate in bpm')
legend('Normalized histogram of heart rate', 'Probability density function')
hold off
% The pdf fits the data well

% categorizing the first half of the data as the athletic group, and the latter half as the non-athletic group
rng(1)
athletes = HeartData(:, 1:50);
non_athletes = HeartData(:, 51:100); 

% initializing variables
athlete_lambda = zeros(1, 50);
non_athlete_lambda = zeros(1, 50);
athlete_mu_min = zeros(1,50);
non_athlete_mu_min = zeros(1,50);
athlete_mu_max = zeros(1,50);
non_athlete_mu_max = zeros(1,50);

% getting lambda of athletes
for i = 1:50
    gm_athlete = fitgmdist(athletes(:,i),2,'Options',statset('MaxIter',1000));
    [mu, c] = min(gm_athlete.mu);
    athlete_lambda(:,i) = gm_athlete.ComponentProportion(1, c);
    athlete_mu_min(:,i) = min(gm_athlete.mu);
    athlete_mu_max(:,i) = max(gm_athlete.mu);
end

% getting the lambda of non-athletes
for i = 1:50
    gm_non_athlete = fitgmdist(non_athletes(:,i),2,'Options',statset('MaxIter',1000));
    [mu, c] = min(gm_non_athlete.mu);
    non_athlete_lambda(:,i) = gm_non_athlete.ComponentProportion(1, c);
    non_athlete_mu_min(:,i) = min(gm_non_athlete.mu);
    non_athlete_mu_max(:,i) = max(gm_non_athlete.mu);
end

% plotting the histogram mixture proportion of the athlete group
figure(2)
subplot(2,1,1);
histogram(athlete_lambda, 15)
title('Histogram of Mixture Proportion of the Athlete Group')
ylabel('Frequency of mixture proportion')
xlabel('Mixture proportion')

  
% plotting the histogram mixture proportion of the non-athlete group
subplot(2,1,2);
histogram(non_athlete_lambda, 15)
title('Histogram of Mixture Proportion of the Non-Athlete Group')
ylabel('Frequency of mixture proportion')
xlabel('Mixture proportion')

% using VARTEST to check if variance between two groups is equal
[h_var,p_var,ci_var,stats_var] = vartest2(athlete_lambda, non_athlete_lambda,'alpha',0.05);
%Since p-value is less than the alpha, their variances are UNEQUAL

% using TTEST2 with UNEQUAL variance to check if there is statistically significant difference
% between two groups
[h_lambda, p_lambda] = ttest2(athlete_lambda, non_athlete_lambda,'alpha',0.05,'vartype','unequal');
% since p-value is higher than alpha (0.05), accept the null hypothesis
% no statistically significant difference

% plotting the histogram of the mean heart rate of the athlete and non-athlete group
figure(3)
subplot(2,1,1);
histogram([athlete_mu_min; athlete_mu_max], 30)
title('Histogram of Mean Heart Rate of the Athlete Group')
ylabel('Frequency of Mean Heart Rate')
xlabel('Mean Heart Rate in bpm')

subplot(2,1,2);
histogram([non_athlete_mu_min;non_athlete_mu_max], 30)
title('Histogram of Mean Heart Rate of the Non-Athlete Group')
ylabel('Frequency of Mean Heart Rate')
xlabel('Mean Heart Rate in bpm')

% finding which cluster bpm 89 belongs to
all_mu = [athlete_mu_min'; athlete_mu_max'; non_athlete_mu_min';non_athlete_mu_max'];
gm_all = fitgmdist(all_mu,4,'Options',statset('MaxIter',1000));
x_all = linspace(min(all_mu), max(all_mu),200);
x_all = x_all';
y_all = pdf(gm_all, x_all);

% getting the Guassian Mixture Distribution between the max and min in each group
gm_athlete = fitgmdist([athlete_mu_min'; athlete_mu_max'],2,'Options',statset('MaxIter',1000));
gm_non_athlete = fitgmdist([non_athlete_mu_min';non_athlete_mu_max'],2,'Options',statset('MaxIter',1000));

% plotting the histogram of the mean of the two groups combined and fitting a PDF curve over it
figure(4)
histogram(all_mu, 30, 'Normalization', 'pdf');
hold on
plot(x_all, y_all)
title('Probability Density Function of the Mean Heart Rate of Athlete and Non-Athlete Groups')
xlabel('Mean Heart Rate in bpm')
ylabel('Probability Density')
hold off

% getting the posterior probability of the Gaussian Mixture Proportion to find which cluster it belongs to
probabilities = posterior(gm_all, 89);
disp(max(probabilities))
disp('The probability 0.97 corresponds to the mean of the non-athlete & exercise cluster.')
