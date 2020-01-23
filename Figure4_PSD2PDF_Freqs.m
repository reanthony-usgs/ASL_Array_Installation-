% This code takes the PSDs output by HHZ_Calculations_1Hr.py and plots them
% as a PSD PDF. Note we must manually make PSD_M by catting together 2
% separate time periods as it is too large to save in memory 

PSD_M = fliplr(PSD_M);

% frequencies 
load('HHZ_freqs_1Hr.mat')

b = length(freq); 
Freqs = flipud(freq);


% Get the percentile statistics 

MedianPSD = prctile(PSD_M,50);
MeanPSD = mean(PSD_M);

p2 = prctile(PSD_M,2.5);
p97 = prctile(PSD_M,97.5);


runs = size(PSD_M,1);


histcent = [-200:.1:-80];
[counts] = hist(PSD_M(:,:), histcent);

%Periods = Periods(1:10:end);


% Make the Peterson curves
fs=250;
dlP=.05;
PSDTOL=15;
[LNMA,HNMA,lpd1,lpd2]=peterson_acc(dlP,fs);

%Smoothed Peterson curves for plotting
NMplotind=(0.001:dlP:10);
LNMAp=spline(lpd1,LNMA,NMplotind);
HNMAp=spline(lpd2,HNMA,NMplotind);

% Remove the Illegal Portions of the High-Frequency Noise Model (Above 10
% Hz)

lpd2 = lpd2(22:end);
HNMA = HNMA(22:end);

pd1 = 10.^(lpd1);
pd2 = 10.^(lpd2);
%% 
figure(21); clf

h = pcolor(Freqs,histcent,log10(counts))
cmap = viridis;
% Make values 0-5 black:
cmap(1,:) = 0.3*ones(1,3);
colormap(cmap);
c=colorbar
xlim([0.002 40])
ylim([-200 -90])
caxis([0 2.3])
set(gca,'FontSize',20)


hold on


H5 = plot(Freqs,MedianPSD,'k');
H2 = plot(1./10.^(lpd1),LNMA,'w:');
H3 = plot(1./10.^(lpd2),HNMA,'w:');
H6 = plot(Freqs,p2,'k:');
H7 = plot(Freqs,p97,'k:');


set(H2,'LineWidth',5.0);
set(H3,'LineWidth',5.0);
set(H5,'LineWidth',3.0);
set(H6,'LineWidth',3.0);
set(H7,'LineWidth',3.0);


%lgd = legend([H5 H6 H2], 'Median PSD', '2.5%/97.5% PSD', 'NHNM/NLNM');
%lgd.Color = [0.7 0.7 0.7];



set(h, 'EdgeColor', 'none');


grid off


axis off
axbot = gca;
set(axbot, 'XScale', 'log', 'YScale', 'linear');
axtop = axes('Position',get(axbot,'Position'),'Color','none',...
            'Xlim',get(axbot,'XLim'), 'Ylim',get(axbot,'YLim'),...
            'XScale', 'log', 'YScale', 'linear' , ...
            'YMinorTick','off' , 'YMinorGrid','off'....
            ) ;
set(gca,'FontSize',20) 

%ticks = get(axtop,'XTickLabel')
%Xticklabels = cellstr(ticks, '10^%d');

%set(axtop,'Xticklabel',Xticklabels)
%XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^%d'));

axtop.LineWidth = 3;  

xlabel('Frequency (Hz)')
ylabel('dB (rel. 1 (m/s^2)^2/Hz)')
ylabel(c,'Log_{10}(Counts)') 











