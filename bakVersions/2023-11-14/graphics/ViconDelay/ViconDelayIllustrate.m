clear;

Ts = 5e-3;
meanTD = 0.034;
cPC = 0.995;

%load ViconDelayIdentArtError.mat
load ViconDelayIdentPropellerOn.mat
AgentNr = 1;

t = meas{AgentNr}.vSerial(1,:);
kk = find( (t>=10.5) & (t<=500) );
t = t(kk);

pENC = meas{AgentNr}.vSerial([2;4;6],kk);
pIMU = meas{AgentNr}.vSerial([3;5;7],kk);
qVIC = meas{AgentNr}.vSerial(8:11,kk);

tVICp = meas{AgentNr}.vSerial(15,kk) - 25;
tREC = meas{AgentNr}.vSerial(17,kk);

% remove redundant vicon frames
[tREC, kk, ~] = unique(tREC);
tVICp = tVICp(kk);

TPC = zeros(size(tREC));
TPC(1) = tREC(1) - tVICp(1);
for kp = 2:length(TPC)
  TPC(kp) = cPC*TPC(kp-1) + (1-cPC)*(tREC(kp) - tVICp(kp));
end

meanTPC = mean(tREC - tVICp);
meanTsV = mean(diff(tVICp));

fig = figure(1); clf;
sp(1) = subplot(3,1,1); 
xlim([30 100]);
for i=1:6
  line(xlim, (meanTsV+i*Ts)*[1 1], 'Linestyle', ':', 'Color', 0.75*[1 1 1]);
  line(xlim, (meanTsV-i*Ts)*[1 1], 'Linestyle', ':', 'Color', 0.75*[1 1 1]);
end
line(xlim, meanTsV*[1 1], 'Linestyle', '-.', 'Color', 'k');

l(1) = line(tREC(2:end), diff(tREC), 'Color', 'r', 'Linestyle', 'none', 'Marker', '.');
l(2) = line(tREC(2:end), diff(tVICp), 'Color', 'g', 'Linestyle', 'none', 'Marker', '.');
ylim([0.07 0.13]);
legend(l(1:2), {'$\rm{diff}(t_{\rm{REC}})$', '$\rm{diff}(t_{\rm{VIC}})$'}, 'Interpreter', 'latex');
ylabel('s');

sp(2) = subplot(3,1,2); 
xlim([30 100]);
for i=1:6
  line(xlim, (meanTPC+i*Ts)*[1 1], 'Linestyle', ':', 'Color', 0.75*[1 1 1]);
  line(xlim, (meanTPC-i*Ts)*[1 1], 'Linestyle', ':', 'Color', 0.75*[1 1 1]);
end
l(1) = line(tREC, TPC, 'Color', 'b', 'Linewidth', 1);
l(2) = line(tREC, tREC - tVICp, 'Color', 'r', 'Linestyle', 'none', 'Marker', '.');
l(3) = line(xlim, meanTPC*[1 1], 'Linestyle', '-.', 'Color', 'k');
ylim([6.666 6.7]);
legend(l([2, 3, 1]), {'$t_{\rm{REC}} - t_{\rm{VIC}}^\prime$', '$\rm{mean}(t_{\rm{REC}} - t_{\rm{VIC}}^\prime)$', '$\hat{T}_{\rm{PC}} + \bar{T}_{\rm{D}}$'}, 'Interpreter', 'latex');
ylabel('s');

sp(3) = subplot(3,1,3);
xlim([30 100]);
for i=1:12
  line(xlim, i*Ts*[1 1], 'Linestyle', ':', 'Color', 0.75*[1 1 1]);
end
l(1) = line(tREC, tREC - (tVICp+TPC) + meanTD, 'Color', 'b', 'Linestyle', 'none', 'Marker', '.');
l(2) = line(xlim, meanTD*[1 1], 'Linestyle', '-.', 'Color', 'k');
ylim([0.025 0.06]);
legend(l(1:2), {'$\hat{T}_{\rm{D}}$', '$\bar{T}_{\rm{D}}$'}, 'Interpreter', 'latex');
ylabel('s');
xlabel('t in s');

%% export figure
figSize = [16 12];
fig.Units = 'centimeters';
fig.Position = [0 0 figSize(1) figSize(2)];
tightfig(fig);
fig.Position = [0 0 figSize(1) figSize(2)];

JokerPrintFig( fig, 'ViconDelayIllustrate', 'pdf', 0 );

