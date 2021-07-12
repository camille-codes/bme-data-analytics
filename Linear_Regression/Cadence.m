clear all
close all
clc

load('GaitData.mat')

%Categorizing participants according to their gender, presence of osteoarthritis, and lifestyle
%female with OA and desk job
count_FAD = 0;
for i = 1:(length(A))
    if (G(i) == 1) && (A(i) == 1) && (D(i) == 1)
        count_FAD = count_FAD + 1;
    end
end

%female with no OA and desk job
count_FOD = 0;
for i = 1:(length(A))
    if (G(i) == 1) && (A(i) == 0) && (D(i) == 1)
        count_FOD = count_FOD + 1;
    end
end

%female with OA and active job
count_FAA = 0;
for i = 1:(length(A))
    if (G(i) == 1) && (A(i) == 1) && (D(i) == 0)
        count_FAA = count_FAA + 1;
    end
end

%female with no OA and active job
count_FOA = 0;
for i = 1:(length(A))
    if (G(i) == 1) && (A(i) == 0) && (D(i) == 0)
        count_FOA = count_FOA + 1;
    end
end

%male with OA and desk job
count_MAD = 0;
for i = 1:(length(A))
    if (G(i) == 0) && (A(i) == 1) && (D(i) == 1)
        count_MAD = count_MAD + 1;
    end
end

%male with no OA and deskjob
count_MOD = 0;
for i = 1:(length(A))
    if (G(i) == 0) && (A(i) == 0) && (D(i) == 1)
        count_MOD = count_MOD + 1;
    end
end

%male with OA and active job
count_MAA = 0;
for i = 1:(length(A))
    if (G(i) == 0) && (A(i) == 1) && (D(i) == 0)
        count_MAA = count_MAA + 1;
    end
end

%male with no OA and active job
count_MOA = 0;
for i = 1:(length(A))
    if (G(i) == 0) && (A(i) == 0) && (D(i) == 0)
        count_MOA = count_MOA + 1;
    end
end

% Using ANOVA to analyze data and plotting using multicompare
[p,tbl,stats] = anovan(Cadence,{G,A,D},'model','interaction','varnames',{'Gender','Osteoarthritis','Job'});
figure(1)
c = multcompare(stats,'Dimension',[1 2 3]);

% Using Linear Regression to analyze data
column_1 = ones(132,1); %design matrix
column_2 = G;
column_3 = A;
column_4 = D;
column_5 = G.*A;
column_6 = G.*D;
column_7 = A.*D;
X = [column_1,column_2,column_3,column_4,column_5, column_6, column_7];

[b_x,bint_x,r_x,rint_x,stats_x] = regress(Cadence,X);

% Demonstrating that the R2 values from ANOVA and Linear Regression are the same,
anova_R = 1 - (3076.1/17864.4);
%Question 5 - Create a new vector that is the Cadence predicted by the regression model in Q2.
%Create a figure that plots the measured Cadence overlaid by predicted Cadence.
theta = inv(X'*X)*(X'*Cadence);
predict_cadence = X*theta;

% Plotting predicted cadence over actual cadence
figure(2)
plot(1:132, predict_cadence)
xlabel('Subject')
ylabel('Cadence [steps/min]')
title('Comparison of Predicted and Actual Cadence')
hold on
plot(1:132, Cadence)
legend('Predicted cadence','Actual cadence')
hold off
