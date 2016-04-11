clear; close all

filename = '/Users/yitaek/Dropbox/Spring 2016/Biology 311/Biology311FinalProject/GEdata.csv';
[geneid,x_40min_H2O2_ura3,x_20min_H2O2_ura3,x0min_H2O2_ura3,x10min_H2O2_ura3,x20min_H2O2_ura3,x40min_H2O2_ura3,x60min_H2O2_ura3,x80min_H2O2_ura3,x_40min_noH2O2_DrosR,x_20min_noH2O2_DrosR,x0min_H2O2_DrosR,x10min_H2O2_DrosR,x20min_H2O2_DrosR,x40min_H2O2_DrosR,x60min_H2O2_DrosR,x80min_H2O2_DrosR] = ...
    importGEdata(filename);

RosRidx = find(ismember(geneid, 'VNG0258H'));

parent = [x_40min_H2O2_ura3,x_20min_H2O2_ura3,x0min_H2O2_ura3,x10min_H2O2_ura3,x20min_H2O2_ura3,x40min_H2O2_ura3,x60min_H2O2_ura3,x80min_H2O2_ura3];
knockout = [x_40min_noH2O2_DrosR,x_20min_noH2O2_DrosR,x0min_H2O2_DrosR,x10min_H2O2_DrosR,x20min_H2O2_DrosR,x40min_H2O2_DrosR,x60min_H2O2_DrosR,x80min_H2O2_DrosR];

figure; clf
plot(parent(205, :), 'g.-')
hold on
plot(knockout(205, :), 'r.-')
set(gca,'XTickLabel', {'-40', '-20', '0', '10', '20', '40', '60', '80'});
xlabel('Time after H_2O_2 Addition (min)')
legend('Parent', 'RosR Deletion')
ylabel('Normalized Log 10 Expression Ratio')

hold off

