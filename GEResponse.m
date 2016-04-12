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

for x=1:length(fig5genes)
    figure; 
    plot(parent(idx(x), :), 'g.-')
    hold on
    plot(knockout(idx(x), :), 'r.-')
    set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});
    xlabel('Time after H_2O_2 Addition (min)')
    legend('Parent', 'RosR Deletion', 'Location', 'Best')
    ylabel('Log_{2} Gene Expression Ratio')

    hold off
end
