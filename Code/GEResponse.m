clear; close all

filename = '/Users/yitaek/Dropbox/Spring 2016/Biology 311/Biology311FinalProject/GEdata.csv';
[geneid,x_40min_H2O2_ura3,x_20min_H2O2_ura3,x0min_H2O2_ura3,x10min_H2O2_ura3,x20min_H2O2_ura3,x40min_H2O2_ura3,x60min_H2O2_ura3,x80min_H2O2_ura3,x_40min_noH2O2_DrosR,x_20min_noH2O2_DrosR,x0min_H2O2_DrosR,x10min_H2O2_DrosR,x20min_H2O2_DrosR,x40min_H2O2_DrosR,x60min_H2O2_DrosR,x80min_H2O2_DrosR] = ...
    importGEdata(filename);

RosRidx = find(ismember(geneid, 'VNG0258H'));

parent = [x_40min_H2O2_ura3,x_20min_H2O2_ura3,x0min_H2O2_ura3,x10min_H2O2_ura3,x20min_H2O2_ura3,x40min_H2O2_ura3,x60min_H2O2_ura3,x80min_H2O2_ura3];
knockout = [x_40min_noH2O2_DrosR,x_20min_noH2O2_DrosR,x0min_H2O2_DrosR,x10min_H2O2_DrosR,x20min_H2O2_DrosR,x40min_H2O2_DrosR,x60min_H2O2_DrosR,x80min_H2O2_DrosR];

figure; clf
plot(parent(205, :), 'k.-')
hold on
plot(knockout(205, :), 'r.-')
set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});
xlabel('Time after H_2O_2 Addition (min)')
legend('Parent', 'RosR Deletion')
ylabel('Normalized Log_{10} Expression Ratio')

hold off

%% Replicating Figure 5
fig5genes = {'VNG0132C', 'VNG0233H', 'VNG0359C', 'VNG0420H', 'VNG0438G', ...
    'VNG0652H', 'VNG0720G', 'VNG0733H', 'VNG0812G', 'VNG0815G', 'VNG0926H',...
    'VNG1042H', 'VNG1182H', 'VNG1262G', 'VNG1615G', 'VNG1926G', 'VNG2108G',...
    'VNG2260H', 'VNG2341H', 'VNG2366C', 'VNG2535H'};

idx = zeros(1, length(fig5genes));

for i=1:length(fig5genes)
    idx(i) = find(ismember(geneid, fig5genes(i)));
end

% %% Drawing the Plot
% for x=1:length(fig5genes)
%     figure;
%     plot(parent(idx(x), :), 'k.-')
%     hold on
%     plot(knockout(idx(x), :), 'r.-')
%     h2o2label = get(gca, 'ylim');
%     plot([3 3], h2o2label, 'b--');
%     set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});
%     xlabel('Time after H_2O_2 Addition (min)')
%     %legend('Parent', '\DeltaVNG0258H', 'H_2O_2 Addition', 'Location', 'Best')
%     ylabel('Log_{2} Gene Expression Ratio')
%     
%     hold off
% end

%% Side Section
figure(50); clf

h1=subplot(3,2,1)
plot(parent(idx(6), :), 'k.-')
hold on
plot(knockout(idx(6), :), 'r.-')
title('A')
h2o2label = get(gca, 'ylim');
plot([3 3], h2o2label, 'b--');
set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});

h2=subplot(3,2,2)
plot(parent(idx(16), :), 'k.-')
hold on
plot(knockout(idx(16), :), 'r.-')
title('D')
h2o2label = get(gca, 'ylim');
plot([3 3], h2o2label, 'b--');
set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});

h3=subplot(3,2,3)
plot(parent(idx(10), :), 'k.-')
hold on
plot(knockout(idx(10), :), 'r.-')
title('F')
h2o2label = get(gca, 'ylim');
plot([3 3], h2o2label, 'b--');
set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});
ylabel('Log_{2} Gene Expression Ratio')

h4=subplot(3,2,4)
plot(parent(idx(19), :), 'k.-')
hold on
plot(knockout(idx(19), :), 'r.-')
title('G')
h2o2label = get(gca, 'ylim');
plot([3 3], h2o2label, 'b--');
set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});
legend('Parent', '\DeltaVNG0258H', 'H_2O_2 Addition', 'Location', 'EastOutside')


h5=subplot(3, 2,5)
plot(parent(idx(11), :), 'k.-')
hold on
plot(knockout(idx(11), :), 'r.-')
title('H')
h2o2label = get(gca, 'ylim');
plot([3 3], h2o2label, 'b--');
set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});

h6=subplot(3,2,6)
plot(parent(idx(8), :), 'k.-')
hold on
plot(knockout(idx(8), :), 'r.-')
title('J')
h2o2label = get(gca, 'ylim');
plot([3 3], h2o2label, 'b--');
set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});

p1=get(h5,'position');
p2=get(h6,'position');
h3=axes('position',[p1(1) p1(2) 0.75  p2(3) ],'visible','off');             
h_label=xlabel('Time after H_2O_2 Addition (min)','visible','on');


%% blurp


figure;
ax1 = subplot(3,1,1);
plot(rand(30,2))
ax2 = subplot(3,1,2);
plot(rand(30,2))
legend({'first','second'},'Location','EastOutside')
ax3 = subplot(3,1,3);
plot(rand(30,2))
p1 = get(ax1,'Position');
p2 = get(ax2,'Position');
p3 = get(ax3,'Position');
set(ax1,'Position',[p2(1) p1(2) p2(3) p1(4)])
set(ax3,'Position',[p2(1) p3(2) p2(3) p3(4)])

