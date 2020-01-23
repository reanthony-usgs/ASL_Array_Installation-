% This code plots the median PSDs for an identical month for the array
% sites

% Load all the data into a giant matrix 

PSD_M = zeros(10,32768);

clear all

% frequencies 
load('HHZ_freqs_1Hr.mat')

% Periods = 1./freq;
freqs = flipud(freq);

%freqs = freq;

% Load the reference spectra

load Median_PSDs/ASA3_HHZ_Full.mat
Reference = MedianPSD;


load SAGE_Month_Medians/ASA1_HHZ.mat
PSD_M(1,:) = MedianPSD;

load SAGE_Month_Medians/ASA2_HHZ.mat
PSD_M(2,:) = MedianPSD;

load SAGE_Month_Medians/ASA3_HHZ.mat
PSD_M(3,:) = MedianPSD;

load SAGE_Month_Medians/ASA4_HHZ.mat
PSD_M(4,:) = MedianPSD;

load SAGE_Month_Medians/ASA5_HHZ.mat
PSD_M(5,:) = MedianPSD;

load SAGE_Month_Medians/ASA6_HHZ.mat
PSD_M(6,:) = MedianPSD;

load SAGE_Month_Medians/ASL9_HHZ.mat
PSD_M(7,:) = MedianPSD;

load SAGE_Month_Medians/VEA1_HHZ.mat
PSD_M(8,:) = MedianPSD;

load SAGE_Month_Medians/ALQ1_HHZ.mat
PSD_M(9,:) = MedianPSD;

load SAGE_Month_Medians/ANMO_HHZ.mat
PSD_M(10,:) = MedianPSD;


Diff_M = PSD_M - Reference;



% ANMO -Blue
% ALQ - Magenta 
% VEA1 - Red
% ASL9 - Blue
% Postholes - grey to black (0.2, 0.4, 0.6, 0.8, 1)
% ASA 4 - CYAN 




figure(5); clf

semilogx(freqs,PSD_M(10,:),'b-','LineWidth',2)
hold on
semilogx(freqs,PSD_M(9,:),'m-','LineWidth',2)
semilogx(freqs,PSD_M(8,:),'r-','LineWidth',2)

xlabel('Frequency (Hz)')
ylabel('dB')
%legend('Mean: 1 Hour Windows','Median: 1 Hour Windows','Mean: 3 Hour Windows','Median: 3 Hour Windows')


set(gca,'FontSize',20) 
xlim([0.002 40])



%%
figure(9); clf
semilogx(freqs,Diff_M(1,:),'-','LineWidth',3, 'color', [0.1,0.1,0.1])
hold on
semilogx(freqs,Diff_M(2,:),'-','LineWidth',2, 'color', [0.25,0.25,0.25])
semilogx(freqs,Diff_M(3,:),'-','LineWidth',4, 'color', [0.40,0.40,0.40])
%semilogx(freqs,Diff_M(3,:),'w-','LineWidth',1)
semilogx(freqs,Diff_M(4,:),'c-','LineWidth',2)
semilogx(freqs,Diff_M(5,:),'r-','LineWidth',2)
semilogx(freqs,Diff_M(6,:),'-','LineWidth',2, 'color', [0.8,0.8,0.8])
semilogx(freqs,Diff_M(10,:),'b-','LineWidth',2)
semilogx(freqs,Diff_M(9,:),'m-','LineWidth',2)
semilogx(freqs,Diff_M(8,:),'r-','LineWidth',2)
semilogx(freqs,Diff_M(7,:),'g-','LineWidth',2)
semilogx(freqs,zeros(length(freqs),1),'k--','LineWidth',5)


xlabel('Frequency (Hz)')
ylabel('dB Difference')
%legend('ASA1', 'ASA2', 'ASA3', 'ASA4', 'ASA5', 'ASA6', 'ANMO 10', 'ALQ1', 'VEA1', 'ASL9')


set(gca,'FontSize',20) 
xlim([0.002 40])
ylim([-20 40])