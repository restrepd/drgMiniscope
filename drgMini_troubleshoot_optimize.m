%Display the lane 1 trials
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])

hold on
seconds_to_plot=5;
ii_seconds_to_plot=ceil(seconds_to_plot/dt_miniscope);
%Lane 1
iil1=6;
iilw1=1;
plot(arena.xsync(trials.ii_laneodor1(iil1):trials.ii_lanewater1(iilw1)),arena.ysync(trials.ii_laneodor1(iil1):trials.ii_lanewater1(iilw1)),'-b')
plot(arena.xsync(trials.ii_laneodor1(iil1)),arena.ysync(trials.ii_laneodor1(iil1)),'ob')
plot(arena.xsync(trials.ii_lanewater1(iilw1)),arena.ysync(trials.ii_lanewater1(iilw1)),'xb')
pffft=1;

fprintf(1,['\nStart lane 1 x ' num2str(arena.xsync(trials.ii_laneodor1(iil1))) ' y ' num2str(arena.ysync(trials.ii_laneodor1(iil1)))  '\n\n'])

%Lane 4
iil4=1;
iilw4=1;
plot(arena.xsync(trials.ii_laneodor4(iil4):trials.ii_lanewater4(iilw4)),arena.ysync(trials.ii_laneodor4(iil4):trials.ii_lanewater4(iilw4)),'-r')
plot(arena.xsync(trials.ii_laneodor4(iil4)),arena.ysync(trials.ii_laneodor4(iil4)),'or')
plot(arena.xsync(trials.ii_lanewater4(iilw4)),arena.ysync(trials.ii_lanewater4(iilw4)),'xr')
pffft=1;
fprintf(1,['\nStart lane 4 x ' num2str(arena.xsync(trials.ii_laneodor1(iil4))) ' y ' num2str(arena.ysync(trials.ii_laneodor1(iil4)))  '\n\n'])

set(gca, 'YDir', 'reverse');
xlabel('x (pixels)')
ylabel('y (pixels)')
title('Lane 1 trajectory1 (blue) and lane 4 trjectory 1 (red)')
