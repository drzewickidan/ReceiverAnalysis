clc
clear all
close all
data = xlsread('shankar_project_spring#2.xls').';
noTarg = sort(data(1:70))';
targ = sort(data(71:100))';



ksdensity(noTarg)
hold on
ksdensity(targ)
[fnoTarg,xnoTarg] = ksdensity(noTarg);
[fTarg,xTarg] = ksdensity(targ);


pdNakagamiNT = fitdist(noTarg,'Nakagami');
pdfNakaNT=pdf(pdNakagamiNT,xnoTarg);

 
pdGammaNT = fitdist(noTarg,'Gamma');
pdfGammaNT=pdf(pdGammaNT,xnoTarg);


pdRaylieghNT = fitdist(noTarg,'Rayleigh');
pdfRaylieghNT = pdf(pdRaylieghNT,xnoTarg);

pdRicianNT = fitdist(noTarg,'Rician');
pdfRicianNT = pdf(pdRicianNT,xnoTarg);

[h_Weibull,p_Weibull,stats_Weibull]=chi2gof(noTarg,'CDF',pdRaylieghNT,'NBins',8,'EMin',1)
[h_Gamma,p_Gamma,stats_Gamma]=chi2gof(noTarg,'CDF',pdGammaNT,'NBins',8,'EMin',1)
[h_Nakagami,p_Nakagami,stats_Nakagami]=chi2gof(noTarg,'CDF',pdNakagamiNT, 'NBins',8,'EMin',1)
[h_Rician,p_Rician,stats_Rician]=chi2gof(noTarg,'CDF',pdRicianNT, 'NBins',8,'EMin',1)

plot(xnoTarg,pdfGammaNT,'--')


pdNakagamiT = fitdist(targ,'Nakagami');
pdfNakaT=pdf(pdNakagamiT,xTarg);



 
pdGammaT = fitdist(targ,'Gamma');
pdfGammaT=pdf(pdGammaT,xTarg);


pdWeibullT = fitdist(targ,'Rayleigh');
pdfWeibullT = pdf(pdWeibullT,xTarg);

pdRicianT = fitdist(targ,'Rician');
pdfRicianT = pdf(pdRicianT,xTarg);

[h_Weibull,p_Weibull,stats_Weibull]=chi2gof(targ,'CDF',pdWeibullT,'NBins',8,'EMin',1)
[h_Gamma,p_Gamma,stats_Gamma]=chi2gof(targ,'CDF',pdGammaT,'NBins',8,'EMin',1)
[h_Nakagami,p_Nakagami,stats_Nakagami]=chi2gof(targ,'CDF',pdNakagamiT, 'NBins',8,'EMin',1)
[h_Rician,p_Rician,stats_Rician]=chi2gof(targ,'CDF',pdRicianT, 'NBins',8,'EMin',1)

plot(xTarg,pdfRicianT,'--')

legend('Data (Target Absent)','Data (Target Present)','Theoretical Fit (Target Absent)','Theoretical Fit (Target Present)');
xlabel('Values')
ylabel('PDF')

figure

z = [zeros(70,1).'  ones(30,1).']; %array of 70's for noTarg and 30 1's for Targ
dataMat=[z;data];
B=(sortrows(dataMat',-2))'; % 2x100 matrix in descending order with each value assigned 1 for Targ and 0 for noTarg

num0=[];  % a cumulative array for number of 0's 
num1=[];  % a cumulative array for number of 1's 
zeroMat=B(1,:);
for i=1:numel(data)+1    % This loop is used to find the values for ROC
    temp =zeroMat(1:i-1);
    num0=[num0 sum(temp<1)];
    num1=[num1 sum(temp>0)];
end
num0 = num0./70;    
num1 = num1./30;
plot(num0,num1)   % This is the plot for the ROC
hold on
%figure

pd = 0.8333;
pfa = 0.1143;

snr_min = albersheim(pd,pfa);
[theor_pd,theor_pfa] = rocsnr(snr_min,'SignalType',...
    'Real');
semilogx(theor_pfa,theor_pd)

data2 = [pdfNakaNT pdfRicianT];
z2 = [zeros(100,1).'  ones(100,1).']; %array of 70's for noTarg and 30 1's for Targ
dataMat2=[z2;data2];
B2=(sortrows(dataMat2',-2))'; % 2x100 matrix in descending order with each value assigned 1 for Targ and 0 for noTarg
 
num0=[];  % a cumulative array for number of 0's 
num1=[];  % a cumulative array for number of 1's 
zeroMat2=B2(1,:);
for i=1:numel(data2)+1    % This loop is used to find the values for ROC
     temp =zeroMat2(1:i-1);
     num0=[num0 sum(temp<1)];
     num1=[num1 sum(temp>0)];
end

% pd=[];
% pf=[];
% for i=i:numel(xTarg))
%     temp1 = pdfRicianT(pdfRicianT>)
%    pd = [pd trapz(fliplr(pdfRicianT(1:i)))];
%    pf = [pf trapz(fliplr(pdfNakaNT(1:i)))];
% end
% pd = pd./max(pd);
% pf = pf./max(pf);
% 
% num0 = num0./100;    
% num1 = num1./100;
% plot(pd,pf)   % This is the plot for the ROC
% 
% 
% 





