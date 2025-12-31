function plotwps(solved, solvernums, titlename, makecodetime, scdtimelimits, scdstepslimits)
% Plot work-precision-steps data in solved structure from wpsStructure run
% solvernums is vector of solver indices to plot (usually 1:5 or 1:6)
% Richard Neidinger 12/3/25

linespec = ["-","-o","-^","-*","-.+","-square"]; % up to 6 curves
close all % new figures in new window

% Figure 1: work (time) on y axis vs precision (scd) on x axis
hold on; 
box on;
for i = solvernums
    plot(solved(i).scd,solved(i).time,linespec(i));
end
plot([0,15],[makecodetime,makecodetime],"--");
legend(solved(solvernums).name,'psm coding','Location','northwest');
% plot degree used as point marker on odepsmJZ 
text(solved(1).scd,solved(1).time,arrayfun(@num2str,solved(1).deg,'UniformOutput',false));
title(titlename);
ylim("padded");
xlim([0,1.07*max(10,solved(1,end).scd(end))])
if nargin >= 5
    axis(scdtimelimits);
end
xticks(0:1:15);
xlabel('precision: scd (significant correct digits)');
ylabel('work: CPU time in seconds');

% Figure 2:  precision (scd) on x axis vs number of steps on y axis
figure
hold on; 
box on;
for i = solvernums
    plot(solved(i).scd,solved(i).steps,linespec(i));
end
legend(solved(solvernums).name,'Location','northwest');
% plot degree used as point marker on odepsmJZ 
text(solved(1).scd,solved(1).steps,arrayfun(@num2str,solved(1).deg,'UniformOutput',false));
title(titlename);
ylim("padded");
xlim([0,1.07*max(10,solved(1,end).scd(end))])
if nargin >= 6
    axis(scdstepslimits);
end
xticks(0:1:15);
xlabel('precision: scd (significant correct digits)');
ylabel('number of solver steps')

% Figure 3: requested scd (-log10(tolerances)) on xaxis vs scd - (requested scd) on yaxis
figure
hold on; 
box on;
scdrequested = -log10(solved(1).tols);
for i = solvernums
    plot(scdrequested,solved(i).scd - scdrequested,linespec(i));
end
yline(0,"--k") 
legend(solved(solvernums).name,'Location','northeast');
% plot degree used as point marker on odepsmJZ 
text(scdrequested,solved(1).scd - scdrequested,arrayfun(@num2str,solved(1).deg,'UniformOutput',false));
title(titlename);
ylim([-5,3]);
xlim([2,14])
if nargin >= 7
    axis(scdexcesslimits);
end
xticks(2:14);
xlabel('scd requested');
ylabel('scd - (scd requested)');

end