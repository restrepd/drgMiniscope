%drgMini_dFFprediction_analysis
clear all
warning('off')

if exist('handles_choices')==0
    clear all
    close all


    %First troubleshooting files
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220804_FCM22/';
    pred_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm_deczdFFopt2_711103.mat';

    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220526_FCM6_withodor_lane4/';
    % pred_file='20220526_FCM6withodor_odorarena_L4_sync_mm_deczdFFopt2_711103.mat';

    % pred_file='20220804_FCM22withodor_odorarena_L1andL4_sync_mm_deczdFFopt_7012103.mat';


    % %Second troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220713_FCM6/';
    % pred_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm_deczdFF_1012103.mat';

    % %Third troubleshooting files
    % this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/PreProcessed/20220727_FCM19/';
    % % pred_file='20220713_FCM6withodor_odorarena_L1andL4_syn_mm_decdFF_1012103.mat';
    % pred_file='20220727_FCM19withodor_odorarena_L1andL4_sync_mm_deczdFF_1012103.mat';

    
 
else
end
figNo=0;
no_to_plot=20;
ii_repeat=1;

load([this_path pred_file])

Rxys_above_95=handles_out.Rxys_above_95;
iiROI_Rxys=handles_out.iiROI_Rxys;
Rops_above_95=handles_out.Rops_above_95;
iiROI_Rops=handles_out.iiROI_Rops;
Rxyops_above_95=handles_out.Rxyops_above_95;
iiROI_Rxyops=handles_out.iiROI_Rxyops;
per_ROI=handles_out.per_ROI;
per_ROI_sh=handles_out.per_ROI_sh;
Rops_for_Rxys_above_95=handles_out.Rops_for_Rxys_above_95;
Rxys_for_Rops_above_95=handles_out.Rxys_for_Rops_above_95;
Rxys_for_Rxyops_above_95=handles_out.Rxys_for_Rxyops_above_95;
Rxyops_for_Rxys_above_95=handles_out.Rxyops_for_Rxys_above_95;

information_content=handles_out.information_content;
sparsity=handles_out.sparsity;
information_contentl1=handles_out.information_contentl1;
sparsityl1=handles_out.sparsityl1;
information_contentl4=handles_out.information_contentl4;
sparsityl4=handles_out.sparsityl4;
spatial_rhol1l4=handles_out.spatial_rhol1l4;
delta_center_of_mass=handles_out.delta_center_of_mass;
sh_information_content=handles_out.sh_information_content;
sh_sparsity=handles_out.sh_sparsity;
sh_information_contentl1=handles_out.sh_information_contentl1;
sh_sparsityl1=handles_out.sh_sparsityl1;
sh_information_contentl4=handles_out.sh_information_contentl4;
sh_sparsityl4=handles_out.sh_sparsityl4;
sh_spatial_rhol1l4=handles_out.sh_spatial_rhol1l4;
sh_delta_center_of_mass=handles_out.sh_delta_center_of_mass;

no_neurons=length(information_content);
p_value_lane_trial=[];
p_value_lanexy_trial=[];
p_value_xy=[];
p_value_xyl1=[];
p_value_xyl4=[];
for ii_ROI=1:no_neurons
    p_value_lane_trial=[p_value_lane_trial handles_out.glm_l14_pvalues.ROI(ii_ROI).pValues(2)];
    p_value_lanexy_trial=[p_value_lanexy_trial handles_out.glm_pvalues.ROI(ii_ROI).pValues(3)];
    p_value_xy=[p_value_xy handles_out.glm_pvalues.ROI(ii_ROI).pValues(2)];
    p_value_xyl1=[p_value_xy handles_out.glm_l1_pvalues.ROI(ii_ROI).pValues(2)];
    p_value_xyl4=[p_value_xy handles_out.glm_l4_pvalues.ROI(ii_ROI).pValues(2)];
end

%Calculate the significance of spatial information (SSI).

if ~isfield(handles_choices,'process_these_ROIs')
    handles_choices.process_these_ROIs=[1:no_neurons];
end
SSI=zeros(1,no_neurons);
SSIl1=zeros(1,no_neurons);
SSIl4=zeros(1,no_neurons);
for iiROI=1:no_neurons
    SSI(iiROI)=(information_content(iiROI)-mean(sh_information_content(iiROI,:)))/std(sh_information_content(iiROI,:));
    SSIl1(iiROI)=(information_contentl1(iiROI)-mean(sh_information_contentl1(iiROI,:)))/std(sh_information_contentl1(iiROI,:));
    SSIl4(iiROI)=(information_contentl4(iiROI)-mean(sh_information_contentl4(iiROI,:)))/std(sh_information_contentl4(iiROI,:));
end

fprintf(1,['Number of ROIs with SSI>1 ' num2str(sum(SSI>=3)) ' out of ' num2str(no_neurons) '\n\n'])
fprintf(1,['Number of ROIs with SSIl1>1 ' num2str(sum(SSIl1>=3)) ' out of ' num2str(no_neurons) '\n\n'])
fprintf(1,['Number of ROIs with SSIl4>1 ' num2str(sum(SSIl4>=3)) ' out of ' num2str(no_neurons) '\n\n'])

%SSI1 vs SSI4
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);
hold on

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])

for ii_ROI=1:no_neurons
    if p_value_lane_trial(ii_ROI)<drsFDRpval(p_value_lane_trial)
        if (p_value_xy(ii_ROI)<drsFDRpval(p_value_xy))
            plot(log10(SSIl1(ii_ROI)),log10(SSIl4(ii_ROI)),'oy')
        else
            plot(log10(SSIl1(ii_ROI)),log10(SSIl4(ii_ROI)),'or')
        end
    else
        if (p_value_xy(ii_ROI)<drsFDRpval(p_value_xy))
            plot(log10(SSIl1(ii_ROI)),log10(SSIl4(ii_ROI)),'ob')
        else
            plot(log10(SSIl1(ii_ROI)),log10(SSIl4(ii_ROI)),'ok')
        end
    end
end

title('SSI lane 4 vs lane 1, red lane trial sig, blue space sig')
ylabel('SSI4')
xlabel('SSI1')


%Do histograms for spatial information content to define the space cells
%We use SSI as defined by Fusi https://doi.org/10.1016/j.neuron.2020.05.022
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on
edges=[-5:30];
% histogram(SSI,edges)
histogram(SSIl1,edges)
histogram(SSIl4,edges)

this_ylim=ylim;
plot([3 3],this_ylim,'-k')

title('Significance of spatial information')
ylabel('No ROIs')
xlabel('SSI')

%Plot a cumulative histogram
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on

% histogram(SSI,edges)
[f_ssil1,x_ssil1] = drg_ecdf(log10(SSIl1+5));
plot(x_ssil1,f_ssil1,'Color',[1 0 0],'LineWidth',3)

[f_ssil4,x_ssil4] = drg_ecdf(log10(SSIl4+5));
plot(x_ssil4,f_ssil4,'Color',[0 0 1],'LineWidth',3)


this_ylim=ylim;
plot([3 3],this_ylim,'-k')

title('Cumulative histogram for SSI (red SSIl1, blue SSIl4')
ylabel('Probability')
xlabel('log10(SSI+5)')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on
edges=[-5:30];
% % histogram(SSI,edges)
% histogram(SSIl1(p_value_lane_trial<drsFDRpval(p_value_lane_trial)),edges)
% histogram(SSIl4(p_value_lane_trial<drsFDRpval(p_value_lane_trial)),edges)

histogram(SSIl1(p_value_xy<drsFDRpval(p_value_xy)),edges)
histogram(SSIl4(p_value_xy<drsFDRpval(p_value_xy)),edges)

this_ylim=ylim;
plot([3 3],this_ylim,'-k')

% title('Significance of spatial information (p_lane<pFDR')
title('Significance of spatial information (p_y<pFDR')
ylabel('No ROIs')
xlabel('SSI')

%Now get the correlation between predicted and actual dFFs
Rop=[];
Rxy=[];
Rxyop=[];

for iiROI=handles_choices.process_these_ROIs
    Rop(iiROI)=per_ROI(iiROI).results.Rop;
    Rxy(iiROI)=per_ROI(iiROI).results.Rxy;
    Rxyop(iiROI)=per_ROI(iiROI).results.Rxyop;
end

Ropl1=[];
Rxyl1=[];
Rxyopl1=[];

for iiROI=handles_choices.process_these_ROIs
    Ropl1(iiROI)=per_ROI(iiROI).results.Ropl1;
    Rxyl1(iiROI)=per_ROI(iiROI).results.Rxyl1;
    Rxyopl1(iiROI)=per_ROI(iiROI).results.Rxyopl1;
end

Ropl4=[];
Rxyl4=[];
Rxyopl4=[];

for iiROI=handles_choices.process_these_ROIs
    Ropl4(iiROI)=per_ROI(iiROI).results.Ropl4;
    Rxyl4(iiROI)=per_ROI(iiROI).results.Rxyl4;
    Rxyopl4(iiROI)=per_ROI(iiROI).results.Rxyopl4;
end

%Do a histogram for spatial rho
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on
edges=[-.5:0.05:1];
histogram(spatial_rhol1l4(handles_choices.process_these_ROIs),edges)
title('Spatial rho l1 vs l4')
xlabel('rho')
ylabel('Number of ROIs')

rho95=prctile(sh_spatial_rhol1l4(handles_choices.process_these_ROIs),95);
rho5=prctile(sh_spatial_rhol1l4(handles_choices.process_these_ROIs),5);

% Define the region for the faded blue background
these_ylims=ylim;
xRegion = [rho95, 1, 1, rho95];
yRegion = [these_ylims(1), these_ylims(1), these_ylims(2), these_ylims(2)]; 

% Add a faded blue patch
patch(xRegion, yRegion, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Define the region for the faded blue background
these_ylims=ylim;
xRegion = [-0.6, rho5, rho5, -0.6];
yRegion = [these_ylims(1), these_ylims(1), these_ylims(2), these_ylims(2)]; 

% Add a faded blue patch
patch(xRegion, yRegion, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

%Plot histogram of center of mass differences 
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on
edges=[0:600/20:600];
histogram(delta_center_of_mass(handles_choices.process_these_ROIs),edges)
title('Difference in center of mass')
xlabel('delta CM')
ylabel('Number of ROIs')

rho95=prctile(sh_delta_center_of_mass(handles_choices.process_these_ROIs),95);
rho5=prctile(sh_delta_center_of_mass(handles_choices.process_these_ROIs),5);

% Define the region for the faded blue background
these_ylims=ylim;
xRegion = [rho95, 600, 600, rho95];
yRegion = [these_ylims(1), these_ylims(1), these_ylims(2), these_ylims(2)]; 

% Add a faded blue patch
patch(xRegion, yRegion, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% Define the region for the faded blue background
these_ylims=ylim;
xRegion = [0, rho5, rho5, 0];
yRegion = [these_ylims(1), these_ylims(1), these_ylims(2), these_ylims(2)]; 

% Add a faded blue patch
patch(xRegion, yRegion, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

%Do histogram of relative dFF change
rel_dFF_change=[];

for this_ROI=handles_choices.process_these_ROIs
    rel_dFF_change(this_ROI)=abs(mean(per_ROI(this_ROI).results.all_dFFl4)-mean(per_ROI(this_ROI).results.all_dFFl1))/...
        (mean(per_ROI(this_ROI).results.all_dFFl4)+mean(per_ROI(this_ROI).results.all_dFFl1));
end


%Do center of mass differences vs rho
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on

for ii_ROI=1:no_neurons
    if p_value_lane_trial(ii_ROI)<drsFDRpval(p_value_lane_trial)
        plot(delta_center_of_mass(ii_ROI),spatial_rhol1l4(ii_ROI),'or')
    else
        plot(delta_center_of_mass(ii_ROI),spatial_rhol1l4(ii_ROI),'ob')
    end

    % if p_value_lane_trial(ii_ROI)<drsFDRpval(p_value_lane_trial)
    %     if (p_value_x(ii_ROI)<drsFDRpval(p_value_x))||(p_value_y(ii_ROI)<drsFDRpval(p_value_y))
    %         plot(delta_center_of_mass(ii_ROI),spatial_rhol1l4(ii_ROI),'oy')
    %     else
    %         plot(delta_center_of_mass(ii_ROI),spatial_rhol1l4(ii_ROI),'or')
    %     end
    % else
    %     if (p_value_x(ii_ROI)<drsFDRpval(p_value_x))||(p_value_y(ii_ROI)<drsFDRpval(p_value_y))
    %         plot(delta_center_of_mass(ii_ROI),spatial_rhol1l4(ii_ROI),'ob')
    %     else
    %         plot(delta_center_of_mass(ii_ROI),spatial_rhol1l4(ii_ROI),'ok')
    %     end
    % end
end

% plot(delta_center_of_mass(handles_choices.process_these_ROIs),spatial_rhol1l4(handles_choices.process_these_ROIs),'ob')
title('Rho vs difference in center of mass')
xlabel('delta CM')
ylabel('Rho')


%Do Rop vs Rho
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on

plot(Rop(handles_choices.process_these_ROIs),spatial_rhol1l4(handles_choices.process_these_ROIs),'ok')
title('Spatial rho l1 vs l4 vs Rho op')
ylabel('Spatial rho')
xlabel('Rho op')

%SSI vs spatial rho
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on

plot(SSI(handles_choices.process_these_ROIs),spatial_rhol1l4(handles_choices.process_these_ROIs),'ok')
title('Spatial rho l1 vs SSI')
ylabel('Spatial rho')
xlabel('SSI')


figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on
plot(delta_center_of_mass(handles_choices.process_these_ROIs), Rop(handles_choices.process_these_ROIs),'ob')
hold on
% plot(sh_delta_center_of_mass, sh_spatial_rhol1l4,'or')
title('Rop vs delta center of mass')
ylabel('Rop')
xlabel('Delta center of mass (mm)')

figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on
plot(SSIl1(handles_choices.process_these_ROIs), SSIl4(handles_choices.process_these_ROIs),'ob')
hold on

title('SSI l1 vs l4 ')
ylabel('SSI l4')
xlabel('SSI l1')



%Display Rho op vs Rho xy figure for this ii_wf and fileNo
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

%Plot Rops above 95
for ii_rho1=1:length(Rops_above_95)
    match_found=0;
    for ii_rho2=1:length(Rxys_above_95)
        if (iiROI_Rxys(ii_rho2)==iiROI_Rops(ii_rho1))
            match_found=1;
        end
    end
    if match_found==1
        plot(Rxys_for_Rops_above_95(ii_rho1),Rops_above_95(ii_rho1),'ok')
    else
        plot(Rxys_for_Rops_above_95(ii_rho1),Rops_above_95(ii_rho1),'ob')
    end
end

%Plot Rxys above 95
for ii_rho1=1:length(Rxys_above_95)
    match_found=0;
    for ii_rho2=1:length(Rops_above_95)
        if (iiROI_Rops(ii_rho2)==iiROI_Rxys(ii_rho1))
            match_found=1;
        end
    end
    if match_found==1
        plot(Rxys_above_95(ii_rho1),Rops_for_Rxys_above_95(ii_rho1),'ok')
    else
        plot(Rxys_above_95(ii_rho1),Rops_for_Rxys_above_95(ii_rho1),'or')
    end
end

this_xlim=xlim;
this_ylim=ylim;
text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.15*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rop>95%','Color',[0 0 1])
text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.10*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy>95%','Color',[1 0 0])
text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.05*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy and Rop>95%','Color',[0 0 0])


plot([this_xlim(1) this_xlim(2)],[this_xlim(1) this_xlim(2)],'-k')

xlabel('Rho xy')
ylabel('Rho op')
title(['Rho op vs. Rho xy log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha) ])

%Display Rho xyop vs Rho xy figure for this ii_wf and fileNo
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.1 .1 .3 .3])


hold on

%Plot Rxyops above 95
for ii_rho1=1:length(Rxyops_above_95)
    match_found=0;
    for ii_rho2=1:length(Rxys_above_95)
        if (iiROI_Rxys(ii_rho2)==iiROI_Rxyops(ii_rho1))
            match_found=1;
        end
    end
    if match_found==1
        plot(Rxys_for_Rxyops_above_95(ii_rho1),Rxyops_above_95(ii_rho1),'ok')
    else
        plot(Rxys_for_Rxyops_above_95(ii_rho1),Rxyops_above_95(ii_rho1),'ob')
    end
end

%Plot Rxys above 95
for ii_rho1=1:length(Rxys_above_95)
    match_found=0;
    for ii_rho2=1:length(Rxyops_above_95)
        if (iiROI_Rxyops(ii_rho2)==iiROI_Rxys(ii_rho1))
            match_found=1;
        end
    end
    if match_found==1
        plot(Rxys_above_95(ii_rho1),Rxyops_for_Rxys_above_95(ii_rho1),'ok')
    else
        plot(Rxys_above_95(ii_rho1),Rxyops_for_Rxys_above_95(ii_rho1),'or')
    end
end

this_xlim=xlim;
this_ylim=ylim;
text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.15*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxyop>95%','Color',[0 0 1])
text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.10*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy>95%','Color',[1 0 0])
text(this_xlim(1)+0.7*(this_xlim(2)-this_xlim(1)), 0.05*(this_ylim(2)-this_ylim(1))+this_ylim(1),'Rxy and Rxyop>95%','Color',[0 0 0])

plot([this_xlim(1) this_xlim(2)],[this_xlim(1) this_xlim(2)],'-k')

xlabel('Rho xy')
ylabel('Rho xyop')
title(['Rho xyop vs. Rho xy log10= ' num2str(handles_choices.weber_fechner) ' alpha= ' num2str(handles_choices.alpha) ])



%Plot the traces and predicted traces
%Rxys above 95 percentile
iiROI_Rxys_to_plot=[];
if ~isempty(iiROI_Rxys)
    to_sort=[Rxys_above_95' iiROI_Rxys'];
    sorted_rows=sortrows(to_sort,1);
    iiROI_Rxys_to_plot=zeros(1,length(iiROI_Rxys));
    iiROI_Rxys_to_plot(1,:)=sorted_rows(:,2);

    these_alldFF=[];
    
    for ii=1:length(iiROI_Rxys)
        these_alldFF=[these_alldFF per_ROI(iiROI_Rxys(ii)).results.all_dFF];
    end

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
    hold on

    % Determine the y spacing of the traces
    y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));


    for iiROI=1:length(iiROI_Rxys)
        this_ROI=iiROI_Rxys_to_plot(iiROI);
        no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
        no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);

        plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
        plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)

        plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xy+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
        plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xy+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
    end

    %Show the last few
    ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])

    xlabel('time(sec)')
    title(['All dFF timecourses for ROIs with xy predictions above 95 percentile'])
end


%Show the fits for xyrops
iiROI_Rxyops_to_plot=[];
if ~isempty(iiROI_Rxyops)
    to_sort=[Rxyops_above_95' iiROI_Rxyops'];
    sorted_rows=sortrows(to_sort,1);
    iiROI_Rxyops_to_plot=zeros(1,length(iiROI_Rxyops));
    iiROI_Rxyops_to_plot(1,:)=sorted_rows(:,2);

    

    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
    hold on

    % Determine the y spacing of the traces
    y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));


    for iiROI=1:length(iiROI_Rxyops_to_plot)
        this_ROI=iiROI_Rxyops_to_plot(iiROI);
        no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
        no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);

        plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
        plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)

        plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xyop+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
        plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xyop+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
    end

    %Show the last few
    ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])

    xlabel('time(sec)')
    title(['All dFF timecourses for ROIs with xyop predictions above 95 percentile'])
end



%Show the fits for ops
iiROI_Rops_to_plot=[];
if ~isempty(iiROI_Rops)
    to_sort=[Rops_above_95' iiROI_Rops'];
    sorted_rows=sortrows(to_sort,1);
    iiROI_Rops_to_plot=zeros(1,length(iiROI_Rops));
    iiROI_Rops_to_plot(1,:)=sorted_rows(:,2);


    figNo=figNo+1;
    try
        close(figNo)
    catch
    end

    hFig = figure(figNo);

    set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
    hold on

    % Determine the y spacing of the traces
    y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));


    for iiROI=1:length(iiROI_Rops_to_plot)
        this_ROI=iiROI_Rops_to_plot(iiROI);
        no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
        no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);

        plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
        plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)

        plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_op+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
        plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_op+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
    end

    %Show the last few
    ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])

    xlabel('time(sec)')
    title(['All dFF timecourses for ROIs with op predictions above 95 percentile'])
end

% these_alldFF=[];
% ii_repeat=1;
% for ii=1:length(iiROI_Rxyops)
%     these_alldFF=[these_alldFF per_ROI(iiROI_Rxyops(ii)).results.all_dFF];
% end
% 
% %Show the fits for xyrops
% to_sort=[Rxyops_above_95' iiROI_Rxyops'];
% sorted_rows=sortrows(to_sort,1);
% iiROI_Rxyops_to_plot=zeros(1,length(iiROI_Rxyops));
% iiROI_Rxyops_to_plot(1,:)=sorted_rows(:,2);
% 
% to_sort=[Rxys_above_95' iiROI_Rxys'];
% sorted_rows=sortrows(to_sort,1);
% iiROI_Rxys_to_plot=zeros(1,length(iiROI_Rxys));
% iiROI_Rxys_to_plot(1,:)=sorted_rows(:,2);
% 
% to_sort=[Rops_above_95' iiROI_Rops'];
% sorted_rows=sortrows(to_sort,1);
% iiROI_Rops_to_plot=zeros(1,length(iiROI_Rops));
% iiROI_Rops_to_plot(1,:)=sorted_rows(:,2);
% 
% % ii_repeat=this_Rxyop_repeat_ii;
% 
% 
% figNo=figNo+1;
% try
%     close(figNo)
% catch
% end
% 
% hFig = figure(figNo);
% 
% set(hFig, 'units','normalized','position',[.05 .1 .85 .8])
% hold on
% 
% % Determine the y spacing of the traces
% y_shift=6*(prctile(these_alldFF(:),95)-prctile(these_alldFF(:),5));
% 
% 
% for iiROI=1:length(iiROI_Rxyops_to_plot)
%     this_ROI=iiROI_Rxyops_to_plot(iiROI);
%     no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
%     no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);
% 
%     plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*iiROI,'-k','LineWidth',1.5)
%     plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*iiROI,'-k','LineWidth',1.5)
% 
%     plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xyop+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
%     plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xyop+y_shift*iiROI-0.3*y_shift,'-r','LineWidth',1)
% end
% 
% %Now place lines between trials
% 
% %Lane 1
% for ii_tr=1:length(per_ROI(this_ROI).results.all_trl1)
%     this_x=per_ROI(this_ROI).results.all_trl1(ii_tr);
%     plot([this_x this_x], [0 y_shift*(iiROI+1)],'-','Color',[0.6 0.6 0.6])
% end
% 
% %Lane 4
% for ii_tr=1:length(per_ROI(this_ROI).results.all_trl4)
%     this_x=no_time_bins1+100+per_ROI(this_ROI).results.all_trl4(ii_tr);
%     plot([this_x this_x], [0 y_shift*(iiROI+1)],'-','Color',[0.6 0.6 0.6])
% end
% 
% %Show the last few
% ylim([y_shift*(iiROI-10) y_shift*(iiROI+2)])
% 
% xlabel('time(sec)')
% title(['All dFF timecourses for ROIs with xyop predictions above 95 percentile'])

%Now show the pseudocolor activity plots
x=25:50:475;
y=24:48:456; 

%Keep track of all iiROIs
% all_iiROIs_above_95=unique([iiROI_Rxyops_to_plot iiROI_Rxys_to_plot iiROI_Rops_to_plot]);
% all_Rs_above_95=zeros(1,length(all_iiROIs_above_95));
% ii_all_iiROIs=0;
% groups_all_iiROIs_above_95=zeros(length(all_iiROIs_above_95),3); %First column op, second xy, third xyop
% information_content=zeros(length(all_iiROIs_above_95),3); %First column xy, sencond lane1, third lane 2
% sparsity=zeros(length(all_iiROIs_above_95),3); %First column xy, sencond lane1, third lane 2



%Now plot all ROIs
figNo=figNo+1;
no_place_cells=0;
place_cells=[];
no_lane_trial_cells=0;
lane_trial_cells=[];

for this_ROI=handles_choices.process_these_ROIs

    %Initialize variables
    this_dFF_activity=zeros(10,10);
    this_dFF_activity_n=zeros(10,10);
    sum_dFF_activity=0;
    
    this_dFFl1_activity=zeros(10,10);
    this_dFFl1_activity_n=zeros(10,10);
    sum_dFFl1_activity=0;

    this_dFFl4_activity=zeros(10,10);
    this_dFFl4_activity_n=zeros(10,10);
    sum_dFFl4_activity=0;

    these_all_dFF=(per_ROI(this_ROI).results.all_dFF)-min(per_ROI(this_ROI).results.all_dFF);
    these_all_lanes=per_ROI(this_ROI).results.all_lanes;
    for ii_t=1:length(per_ROI(this_ROI).results.all_dFF)
        
        this_x_ii=ceil(per_ROI(this_ROI).results.all_XY(ii_t,1)/50);
        if this_x_ii==11
            this_x_ii=10;
        end

        this_y_ii=ceil(per_ROI(this_ROI).results.all_XY(ii_t,2)/48);
        if this_y_ii==11
            this_y_ii=10;
        end

        this_dFF_activity(this_x_ii,this_y_ii)=this_dFF_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
        this_dFF_activity_n(this_x_ii,this_y_ii)=this_dFF_activity_n(this_x_ii,this_y_ii)+1;
        sum_dFF_activity=sum_dFF_activity+these_all_dFF(ii_t);

        if these_all_lanes(ii_t)==1
            this_dFFl1_activity(this_x_ii,this_y_ii)=this_dFFl1_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
            this_dFFl1_activity_n(this_x_ii,this_y_ii)=this_dFFl1_activity_n(this_x_ii,this_y_ii)+1;
            sum_dFFl1_activity=sum_dFFl1_activity+these_all_dFF(ii_t);
        else
            this_dFFl4_activity(this_x_ii,this_y_ii)=this_dFFl4_activity(this_x_ii,this_y_ii)+these_all_dFF(ii_t);
            this_dFFl4_activity_n(this_x_ii,this_y_ii)=this_dFFl4_activity_n(this_x_ii,this_y_ii)+1;
            sum_dFFl4_activity=sum_dFFl4_activity+these_all_dFF(ii_t);
        end
    end

    for ii_x=1:10
        for ii_y=1:10
            if this_dFF_activity_n(ii_x,ii_y)~=0
                this_dFF_activity(ii_x,ii_y)=this_dFF_activity(ii_x,ii_y)/this_dFF_activity_n(ii_x,ii_y);
            end
        end
    end

    for ii_x=1:10
        for ii_y=1:10
            if this_dFFl1_activity_n(ii_x,ii_y)~=0
                this_dFFl1_activity(ii_x,ii_y)=this_dFFl1_activity(ii_x,ii_y)/this_dFFl1_activity_n(ii_x,ii_y);
            end
        end
    end

    for ii_x=1:10
        for ii_y=1:10
            if this_dFFl4_activity_n(ii_x,ii_y)~=0
                this_dFFl4_activity(ii_x,ii_y)=this_dFFl4_activity(ii_x,ii_y)/this_dFFl4_activity_n(ii_x,ii_y);
            end
        end
    end
    %
    % this_dFF_activity=this_dFF_activity/sum_dFF_activity;
    % this_dFFl1_activity=((no_trials/2)/trials.lane1)*this_dFFl1_activity/(sum_dFFl1_activity+sum_dFFl4_activity);
    % this_dFFl4_activity=((no_trials/2)/trials.lane4)*this_dFFl4_activity/(sum_dFFl1_activity+sum_dFFl4_activity);

    
    try
        close(figNo)
    catch
    end

    %Plot the shifted odor plume
    hFig = figure(figNo);
    set(hFig, 'units','normalized','position',[.1 .1 .8 .6])

    %Plot the fits
    no_time_bins1=length(per_ROI(this_ROI).results.all_dFFl1);
    no_time_bins4=length(per_ROI(this_ROI).results.all_dFFl4);

    y_shift=max(per_ROI(this_ROI).results.all_dFF);

    if y_shift==0
        y_shift=1;
    end

    %Plot dFF timecourses
    subplot(2,3,1:3)
    hold on
    ii_plot=0;
   

    

    plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*ii_plot,'-k','LineWidth',1.5)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*ii_plot,'-k','LineWidth',1.5)

    plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_op+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_op+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)
    

    ii_plot=1;
    plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*ii_plot,'-k','LineWidth',1.5)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*ii_plot,'-k','LineWidth',1.5)

    plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xy+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xy+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)

    ii_plot=2;
    plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFl1+y_shift*ii_plot,'-k','LineWidth',1.5)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFl4+y_shift*ii_plot,'-k','LineWidth',1.5)

    plot([1:no_time_bins1], per_ROI(this_ROI).results.all_dFFpredl1_xyop+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)
    plot([no_time_bins1+100:no_time_bins1+99+no_time_bins4], per_ROI(this_ROI).results.all_dFFpredl4_xyop+y_shift*ii_plot-0.3*y_shift,'-r','LineWidth',1)

    %Lane 1
    for ii_tr=1:length(per_ROI(this_ROI).results.all_trl1)
        this_x=per_ROI(this_ROI).results.all_trl1(ii_tr);
        plot([this_x this_x], [-y_shift y_shift*(ii_plot+1)],'-','Color',[0.6 0.6 0.6])
    end

    %Lane 4
    for ii_tr=1:length(per_ROI(this_ROI).results.all_trl4)
        this_x=no_time_bins1+100+per_ROI(this_ROI).results.all_trl4(ii_tr);
        plot([this_x this_x], [-y_shift y_shift*(ii_plot+1)],'-','Color',[0.6 0.6 0.6])
    end

    ylim([-y_shift y_shift*(ii_plot+1)])

    if no_time_bins1==0
        xlim([-100 100+no_time_bins4])
        no_time_bins1=-100;
    end
    %Show Rhos

    %Rop
    % if sum(this_ROI==iiROI_Rops_to_plot)==1
        text(no_time_bins1+10,0,['op ' num2str(per_ROI(this_ROI).results.Rop)], 'Color','k')
    % else
    %     text(no_time_bins1+10,0,['op ' num2str(per_ROI(this_ROI).results.Rop)], 'Color',[0.6 0.6 0.6])
    % end

    %Rxy
    % if sum(this_ROI==iiROI_Rxys_to_plot)==1
        text(no_time_bins1+10,y_shift,['xy ' num2str(per_ROI(this_ROI).results.Rxy)], 'Color','k')
    % else
    %     text(no_time_bins1+10,y_shift,['xy ' num2str(per_ROI(this_ROI).results.Rxy)], 'Color',[0.6 0.6 0.6])
    % end

    %Rxyop
    % if sum(this_ROI==iiROI_Rxyops_to_plot)==1
        text(no_time_bins1+10,2*y_shift,['xyop ' num2str(per_ROI(this_ROI).results.Rxyop)], 'Color','k')
    % else
    %     text(no_time_bins1+10,2*y_shift,['xyop ' num2str(per_ROI(this_ROI).results.Rxyop)], 'Color',[0.6 0.6 0.6])
    % end

    %Plot space activity maps as in Moser https://doi.org/10.1126/science.1114037
    max_this_dFF_activity=max(this_dFF_activity(:));
    max_this_dFFl4_activity=max(this_dFFl4_activity(:));
    max_this_dFFl1_activity=max(this_dFFl1_activity(:));

    % %Entire arena
    % subplot(2,4,5)
    % max_activity=max_this_dFF_activity;
    % delta_ac=max_activity/255;
    % this_masked_dFF_activity=this_dFF_activity;
    % this_masked_dFF_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
    % drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFF_activity)
    % colormap fire
    % this_cmap=colormap;
    % this_cmap(1,:)=[0.2 0.2 0.2];
    % colormap(this_cmap)
    % clim([-1.5*delta_ac max_activity])
    % shading interp
    % set(gca, 'YDir', 'reverse');
    % 
    % yticks(0:48:480)
    % xticks(0:50:500)
    % xlabel('x (mm)')
    % ylabel('y (mm)')
    % title(['All trials, SSI= ' num2str(SSI(this_ROI))])

    colormap fire
    this_cmap=colormap;
    this_cmap(1,:)=[0.3 0.3 0.3];

    if max_this_dFFl1_activity>max_this_dFFl4_activity
        %Lane 1
        subplot(2,3,4)
        max_activity=max_this_dFFl1_activity;
        delta_ac=max_activity/255;
        this_masked_dFFl1_activity=this_dFFl1_activity;
        this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
        drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
        colormap(this_cmap)
        clim([-1.5*delta_ac max_activity])
        shading interp
        set(gca, 'YDir', 'reverse');

        yticks(0:48:480)
        xticks(0:50:500)
        xlabel('x (mm)')
        ylabel('y (mm)')
        title(['Lane 1, SSI= ' num2str(SSIl1(this_ROI))])

        %Lane 4 normalized to lane 1
        subplot(2,3,5)
        max_activity=max_this_dFFl1_activity;
        this_masked_dFFl4_activity=this_dFFl4_activity;
        this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
        drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
        colormap(this_cmap)
        clim([-1.5*delta_ac max_activity])
        shading interp
        set(gca, 'YDir', 'reverse');

        yticks(0:48:480)
        xticks(0:50:500)
        xlabel('x (mm)')
        ylabel('y (mm)')
        title(['Lane 4, SSI= ' num2str(SSIl4(this_ROI))])

        %Lane 4 normalized to lane 4
        subplot(2,3,6)
        max_activity=max_this_dFFl4_activity;
        this_masked_dFFl4_activity=this_dFFl4_activity;
        this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
        drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
        colormap(this_cmap)
        clim([-1.5*delta_ac max_activity])
        shading interp
        set(gca, 'YDir', 'reverse');

        yticks(0:48:480)
        xticks(0:50:500)
        xlabel('x (mm)')
        ylabel('y (mm)')
        title(['Lane 4'])

    else


        %Lane 4 normalized to lane 4
        subplot(2,3,4)
        max_activity=max_this_dFFl4_activity;
        delta_ac=max_activity/255;
        this_masked_dFFl4_activity=this_dFFl4_activity;
        this_masked_dFFl4_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
        drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl4_activity)
        colormap(this_cmap)
        clim([-1.5*delta_ac max_activity])
        shading interp
        set(gca, 'YDir', 'reverse');

        yticks(0:48:480)
        xticks(0:50:500)
        xlabel('x (mm)')
        ylabel('y (mm)')
       title(['Lane 4, SSI= ' num2str(SSIl4(this_ROI))])

        %Lane 1 normalized to lane 4
        subplot(2,3,5)
        max_activity=max_this_dFFl4_activity;
        this_masked_dFFl1_activity=this_dFFl1_activity;
        this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
        drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
        colormap(this_cmap)
        clim([-1.5*delta_ac max_activity])
        shading interp
        set(gca, 'YDir', 'reverse');

        yticks(0:48:480)
        xticks(0:50:500)
        xlabel('x (mm)')
        ylabel('y (mm)')
        title(['Lane 1, SSI= ' num2str(SSIl1(this_ROI))])

       %Lane 1 normalized to lane 1
        subplot(2,3,6)
        max_activity=max_this_dFFl1_activity;
        this_masked_dFFl1_activity=this_dFFl1_activity;
        this_masked_dFFl1_activity(this_dFF_activity_n==0)=-0.9*delta_ac;
        drg_pcolor(repmat(x,length(y),1)',repmat(y,length(x),1),this_masked_dFFl1_activity)
        colormap(this_cmap)
        clim([-1.5*delta_ac max_activity])
        shading interp
        set(gca, 'YDir', 'reverse');

        yticks(0:48:480)
        xticks(0:50:500)
        xlabel('x (mm)')
        ylabel('y (mm)')
        title(['Lane 1'])
    end


    sgt_legend=['dFF map ROI No ' num2str(this_ROI) ' rho ' ...
        num2str(spatial_rhol1l4(this_ROI)) ' dcm ' num2str(delta_center_of_mass(this_ROI))...
        ' SSIl1= ' num2str(SSIl1(this_ROI)) ' SSIl4= ' num2str(SSIl4(this_ROI))];

    if p_value_lane_trial(ii_ROI)<drsFDRpval(p_value_lane_trial)
        sgt_legend=[sgt_legend ' lane sig']
    end

    if (p_value_xy(ii_ROI)<drsFDRpval(p_value_xy))
        sgt_legend=[sgt_legend ' space sig']
    end

    if (p_value_xyl1(ii_ROI)<drsFDRpval(p_value_xyl1))
        sgt_legend=[sgt_legend ' space l1 sig']
    end

    if (p_value_xyl4(ii_ROI)<drsFDRpval(p_value_xyl4))
        sgt_legend=[sgt_legend ' space l4 sig']
    end

    sgtitle(sgt_legend)

    if (abs(spatial_rhol1l4(this_ROI))<=0.15)||(spatial_rhol1l4(this_ROI)>=0.5)||(delta_center_of_mass(this_ROI)<=50)||(delta_center_of_mass(this_ROI)>=350)
        pffft=1;
    end

    if (per_ROI(this_ROI).results.Rop>=0.4)||(per_ROI(this_ROI).results.Rxy>=0.4)||(per_ROI(this_ROI).results.Rxyop>=0.4)
            pffft=1;
    end

    %Now show the place cells
    if (abs(spatial_rhol1l4(this_ROI))>=0.4)&(delta_center_of_mass(this_ROI)<=100)...
            &(p_value_xy(ii_ROI)<drsFDRpval(p_value_xy))&(p_value_xyl1(ii_ROI)<drsFDRpval(p_value_xyl1))&(p_value_xyl4(ii_ROI)<drsFDRpval(p_value_xyl4))
        no_place_cells=no_place_cells+1;
        place_cells=[place_cells this_ROI];
        pffft=1;
    end

    %Now show the odor cells
    if (abs(spatial_rhol1l4(this_ROI))<=0.1)&(delta_center_of_mass(this_ROI)>=300)...
            &(p_value_lane_trial(ii_ROI)<drsFDRpval(p_value_lane_trial))
        no_lane_trial_cells=no_lane_trial_cells+1;
        lane_trial_cells=[lane_trial_cells this_ROI];
        pffft=1;
    end
    
    pffft=1;
end
    
%Let's see where the place cells and odor cells fall with respect to dFF fits

%Plot a cumulative histogram for Rxys
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on

%place cells
[f,x] = drg_ecdf(Rxy);
plot(x,f,'Color',[0 0 0],'LineWidth',3)

[f,x] = drg_ecdf(Rxy(place_cells));
plot(x,f,'Color',[1 0 0],'LineWidth',3)

[f,x] = drg_ecdf(Rxy(lane_trial_cells));
plot(x,f,'Color',[0 0 1],'LineWidth',3)


title('Cumulative histogram for Rxy (red place, blue odor')
ylabel('Probability')
xlabel('Rxy')


%Plot a cumulative histogram for Rxyops
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on

%place cells
[f,x] = drg_ecdf(Rxyop);
plot(x,f,'Color',[0 0 0],'LineWidth',3)

[f,x] = drg_ecdf(Rxyop(place_cells));
plot(x,f,'Color',[1 0 0],'LineWidth',3)

[f,x] = drg_ecdf(Rxyop(lane_trial_cells));
plot(x,f,'Color',[0 0 1],'LineWidth',3)


title('Cumulative histogram for Rxyop (red place, blue odor')
ylabel('Probability')
xlabel('Rxyop')

%Plot a cumulative histogram for Rops
figNo=figNo+1;
try
    close(figNo)
catch
end

hFig = figure(figNo);

set(hFig, 'units','normalized','position',[.05 .1 .3 .3])
hold on

%place cells
[f,x] = drg_ecdf(Rop);
plot(x,f,'Color',[0 0 0],'LineWidth',3)

[f,x] = drg_ecdf(Rop(place_cells));
plot(x,f,'Color',[1 0 0],'LineWidth',3)

[f,x] = drg_ecdf(Rop(lane_trial_cells));
plot(x,f,'Color',[0 0 1],'LineWidth',3)


title('Cumulative histogram for Rxyop (red place, blue odor')
ylabel('Probability')
xlabel('Rop')
pffft=1;