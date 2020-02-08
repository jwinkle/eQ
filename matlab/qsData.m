%%
%  QS DATA VS. FLOW RATE

load '~/Dropbox/xps/eQ/matlab/flowData.mat'

flows = {"5", "10", "50", "100"};
xScale = 5;%trap width is 500 um

dataRange = 2:4;
data_uM = flowData(:,dataRange)/1000;
xData = flowData(:,1)*xScale;
signalThresh = 3.5 * ones(length(xData), 1);


figure(1); clf;

p=plot(xData, data_uM, 'LineWidth', 4);  hold on;
p(1).Color = 'blue';
p(2).Color = 'green';
p(3).Color = 'red';
plot(xData, signalThresh, 'k--', 'LineWidth', 3);
a = area(xData, signalThresh);
a.FaceColor = 'r';
a.FaceAlpha = 0.1;
hold off;

lgdFS = 16;
lgd=legend(flows(1:length(dataRange)), 'Location', 'northwest', 'FontSize', lgdFS);
lgd.Title.String = {'Channel Flow', 'Rate ($\mu m/sec$)'};
lgd.Title.FontSize = lgdFS;
lgd.Interpreter = 'Latex';


tFS = 20;
ax=gca;
ax.FontSize = tFS;
title('HSL CONCENTRATION VS. FLOW RATE', 'FontSize', tFS);
ylabel('PEAK HSL H_P(x) [\muM]', 'FontSize', tFS);
xlabel('Trap x-position [\mum]', 'FontSize', tFS);



%%
%  SENDER-RECEIVER DATA VS. D, h

D = [3e4; 1.6e4;];   %C4, C14 diffusion coeffs.
h = [150; 100; 75; 50; 25;];  %trap y-heighhts

hScales = (D(1)./D);  %how peak HSL will scale with D

gamma_e = 8*D*(h.^-2)';  %effective degradation (peak): 8D/h^2
lambdas = sqrt(h.^2/8);  %effective signaling depth lambda^2 = D/gamma_e

x = linspace(0,250, 100);  %250 microns
% range = 1:3;
range = 2;

l = lambdas(range);

decays = exp(-x' * (1./l)');  %column vector
decay2 = hScales(2) * exp(-x' * (1./l)' * hScales(2));  %column vector
s = hScales;
data = decays * s';  %decay x # scales
% data = decays;

figure(2); clf;
% plot(x, decays);
plot(x, data, 'LineWidth', 4);  hold on;
plot(x, decay2, 'k--', 'LineWidth', 4, 'color', [0 0 0 0.5]);

%primitive search for 1/e points:
d1e1 = find(data(:,1) < data(1,1)/exp(1), 1);
d2e1 = find(data(:,2) < data(1,2)/exp(1), 1);
d3e1 = find(decay2 < decay2(1)/exp(1), 1);
mlw = 2;
plot(x(d1e1), data(d1e1,1), 'ks', 'MarkerSize', 15, 'LineWidth', mlw);
plot(x(d2e1), data(d2e1,2), 'ks', 'MarkerSize', 15, 'LineWidth', mlw);
plot(x(d3e1), decay2(d3e1,1), 'ks', 'MarkerSize', 15, 'LineWidth', mlw);
hold off;

lgdFS = 16;
legs = ["D = 3.0e4 [$\mu m^2/min$]",  "D = 1.8e4", {"D = 1.8e4 with fixed $\gamma_e$"}, '$e^{-1}$ location'];
lgd=legend(legs, 'Location', 'northeast', 'FontSize', lgdFS);
lgd.Title.String = {'DIFFUSION RATES'};
lgd.Title.FontSize = lgdFS;
lgd.Interpreter = 'Latex';


tFS = 20;
ax=gca;
ax.FontSize = tFS;
title('HSL CONCENTRATION VS. D', 'FontSize', tFS);
ylabel('PEAK HSL H_P(x) [\muM]', 'FontSize', tFS);
xlabel('Trap x-position [\mum]', 'FontSize', tFS);



% get(gca,'colororder')
% ans =
%          0    0.4470    0.7410
%     0.8500    0.3250    0.0980
%     0.9290    0.6940    0.1250
%     0.4940    0.1840    0.5560
%     0.4660    0.6740    0.1880
%     0.3010    0.7450    0.9330
%     0.6350    0.0780    0.1840



