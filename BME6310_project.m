clc
clear
close all

%% read xls data
num = xlsread('FirstPassData');

for i = 1:6
    t_all(i,:) = num(:,3*(i-1)+1);
    AIF_all(i,:) = num(:,3*(i-1)+2);
    TF_all(i,:) = num(:,3*(i-1)+3);
%     figure, plot(t_all(1,:),AIF_all(1,:)); xlabel('t/s');ylabel('AIF Gd');
end

group = 3;
t = t_all(group,:);
TF = TF_all(group,:);
AIF = AIF_all(group,:);

%% STEP 1: use gamma variate to fit C_p(t)


points = 55;
t1 = t(1:points);
AIF1 = AIF(1:points);
TF1 = TF(1:points);

%AIF1 = filloutliers(AIF1,'pchip','percentiles',[10 90])

arrival = 3.3;
alpha = 2.5;
beta = 0.2;
A = 7;

z=zeros(size(AIF1));
z(t1>arrival) = 1;

func_t_AIF = @(paras,t1) z*paras(1).*((t1-paras(4)).^paras(2)).*exp(-(t1-paras(4))./paras(3));

initials = [A alpha beta arrival];

gammafitAIF = lsqcurvefit(func_t_AIF,initials,t1,AIF1);
g = real(gammafitAIF);

AIF_gammafit = real(func_t_AIF(g,t1));

figure, 
plot(t1,AIF1,'.');
hold on, plot(t1,AIF_gammafit)
title('Cp vs time fit');
xlabel('t (s)');
ylabel('Gd (mg/ul)');

%% STEP 2: fit TF curve

t1 = t1(20:55);
TF1 = TF1(20:55);

tissuefun = @(x,t1) x(1)*conv(exp(-x(2)*t1),AIF_gammafit,'same');
% tissuefun = @(x,t2) x(1)*AIF_gammafit_long + (1-x(1))*x(2)*conv( AIF_gammafit_long, exp(-x(2)*t2),'valid');

convfitTF = lsqcurvefit(tissuefun,[0.001, 0.01],t1,TF1);
TF_fit = tissuefun(convfitTF,t1);

figure,
plot(t1,TF1, 'r.'), title('Ce vs time fit');
hold on, plot(t1, TF_fit);
xlabel('t (s)');
ylabel('Gd (mg/ul)');
% %----------------------------------------%
%plot Cp and Ce together
figure,
plot(t1,TF1, 'r.'), title('Cp and Ce curves');
hold on, plot(t1, TF_fit,'r-');
plot(t1,AIF1(20:55),'g.');
plot(t1,AIF_gammafit(20:55),'g-')
xlabel('t (s)');
ylabel('Gd (mg/ul)');
legend('Cp raw data', 'Cp fitted curve', 'Ce raw data', 'Ce fitted curve')
%%------------------------------------------%
%% STEP 3: outliers replacement

%use outlier filter to filter out data more than 90 percentile and
%replace with spline interpolation
TF1_out = filloutliers(TF1,'pchip','percentiles',[0 90])

tissuefun = @(x,t1) x(1)*conv(exp(-x(2)*t1),AIF_gammafit,'same');
% tissuefun = @(x,t2) x(1)*AIF_gammafit_long + (1-x(1))*x(2)*conv( AIF_gammafit_long, exp(-x(2)*t2),'valid');

convfitTF_out = lsqcurvefit(tissuefun,[0.001, 0.01],t1,TF1_out);
TF_fit_out = tissuefun(convfitTF_out,t1);
%plot new Ce curve after outlier replacement
figure,
plot(t1,TF1_out, 'r.'), title('Ce vs time fit');
hold on, plot(t1, TF_fit_out);
xlabel('t (s)');
ylabel('Gd (mg/ul)');
%plot Ce curve compare between before and after outlier replacement
figure,
plot(t1, TF_fit_out,'r-'), title('Ce curves of outlier replaced vs original data')
hold on, plot(t1, TF_fit,'b-');
xlabel('t (s)');
ylabel('Gd (mg/ul)');
legend('outlier replaced', 'original data');
legend('Location','southeast')
