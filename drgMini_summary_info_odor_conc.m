%The output of drgMini_summary_info_odor_conc 
save_PathOCMI='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/';
save_FileOCMI='outputPerROIInformationContentOC.mat';

load([save_PathOCMI save_FileOCMI])

figureNo=0;

%Plot histogram for SSI mutual info for both lanes for xy x bindFF
ssi_all_spatial_info=handles_outic.ssi_all_spatial_info;

figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

edges=[-4:2:34];
histogram(ssi_all_spatial_info)

plot([3 3],[0 450],'-k','LineWidth',2)

xlim([-5 35])
title('Histogram for SSI')
xlabel('Sigma')
ylabel('Count')


%Plot histogram for SSI mutual info for op x bindFF
ssi_all_op_info=handles_outic.ssi_all_op_info;

%Plot histogram for SSI mutual info for op x bindFF
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on

edges=[-4:1:20];
histogram(ssi_all_op_info,edges)

plot([3 3],[0 900],'-k','LineWidth',2)

xlim([-5 20])
title('Histogram for SOCI')
xlabel('Sigma')
ylabel('Count')


%Now plot SSI vs SOCI for SSI>+3 or SOCI>=3
include_these=(logical(ssi_all_op_info>=3))|(logical(ssi_all_spatial_info>=3));

%Plot SSI vs. SOCI
figureNo=figureNo+1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
set(hFig, 'units','normalized','position',[.1 .1 .3 .3])
hold on



% light gray rectangle from [0 0] to [0.5 handles_anal_XY.r2ERy_sig_thr]
patch([-5 3 3 -5], ...
      [3 3 25 25], ...
      [0 0.8 0.8], ...          % cyan
      'EdgeColor','none', ...
      'FaceAlpha',0.3);           % optional transparency

% light gray rectangle from [0 0] to [handles_anal_conc.r2ERconc_sig_thr 0.5]
patch([3 3 35 35], ...
      [-5 3 3 -5], ...
      [0 0.8 0.8], ...          %  cyan
      'EdgeColor','none', ...
      'FaceAlpha',0.3);           % optional transparency

plot([-5 35], [3 3],'-','Color',[0 0.7 0.7],'LineWidth',2)
plot([3 3],[-5 25],'-','Color',[0 0.7 0.7],'LineWidth',2)

plot(ssi_all_spatial_info(include_these),ssi_all_op_info(include_these),'.k' )



plot([-5 25],[-5 25],'-','LineWidth',2,'Color',[1 0 1])


ylabel('SOCI')
xlabel('SSI')
xlim([-5 35])
ylim([-5 25])

title('SOCI vs SSI')


pffft=1;
