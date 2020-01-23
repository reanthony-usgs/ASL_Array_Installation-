% This code takes two individual horizontal component PSDs recorded at 
% identical time steps and makes them into single "horizontal power" PSD

clear all

net = 'IU';
station = 'ANMO';

% Load in the data 
path1 = ['./PSD_Database/',net,station,'10HH1_Month.mat'];
load(path1); 
PSD_M1 = PSD_M;

% for ANNMO only 

PSD_M1(719:720,:) = [];

path2 = ['./PSD_Database/',net,station,'10HH2_Month.mat'];
load(path2); 
PSD_M2 = PSD_M;

clear PSD_M


% Process the data to be a single horizontal PSD

Power_1 = 10.^(PSD_M1./10);
Power_2 = 10.^(PSD_M2./10);

H_Power = Power_1 + Power_2;

PSD_M_H = 10*log10(H_Power); 

% Save the metrics that we want to store 

save_path = ['./PSD_Database/Horizontal_PSDs/',net,station,'10HHH_Month.mat'];
save(save_path,'PSD_M_H');

MedianPSD = median(PSD_M_H); 

save_path2 = ['./SAGE_Month_Medians/',station,'_HHH.mat'];

save(save_path2,'MedianPSD');












