%drgMini_analyze_batch_DecodeCorrXYConcSingleROIsWithPrevious
close all
clear all

figureNo=0;

fileID = fopen('/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/agent_stats.txt','w');



%Now plot Hit rate for kappa
kappa_input_file='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Fabio_results_final_agent_navigation/Results_experiments/kappa_results.mat';
load(kappa_input_file)

all_hit_rates_kappa=hit_rates_kappa_all;

sh_kappa_input_file='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Fabio_results_final_agent_navigation/Results_control/kappa_results.mat';
load(sh_kappa_input_file)

sh_all_hit_rates_kappa=hit_rates_kappa_all;

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[0:5:100];
rand_offset=0.5;

glm_hitRate=[];
glm_hitRate_ii=0;

id_hitRate_ii=0;
input_hitRate_data=[];

 
these_groups=[1 5]; %2 cm and 1 cm all trials


these_means_hitRates=[];
these_means_hitRates_sh=[];
bar_offsets=[];

ii_kappa_ordered=[2 1 3 4 5];
kappa_ordered=[1 0 2 5 10];

for kappa_ii=1:length(kappa_values)
    %First get the hitRates for each number of ROIs
    these_hitRates=[];
    these_hitRates=all_hit_rates_kappa(kappa_ii,:);
    these_means_hitRates=[these_means_hitRates mean(these_hitRates)];
    

    %Then get the hitRates_sh for each number of ROIs
    these_hitRate_sh=[];
    these_hitRate_sh=sh_all_hit_rates_kappa(kappa_ii,:);
    these_means_hitRates_sh=[these_means_hitRates_sh mean(these_hitRate_sh)];
    
    bar_offsets=[bar_offsets bar_offset];

    %plot bar
    bar(bar_offset,mean(these_hitRates),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])

    %Violin plot
    try
        [mean_out, CIout,violin_x]=drgViolinPoint(these_hitRates...
            ,edges,bar_offset,rand_offset,'k','k',5);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;



    %plot bar
    bar(bar_offset,mean(these_hitRate_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])

    %Violin plot
    try
        [mean_out, CIout,violin_x]=drgViolinPoint(these_hitRate_sh...
            ,edges,bar_offset,rand_offset,'k','k',5);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    glm_hitRate.data(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=these_hitRates;
    glm_hitRate.shuffled(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=0*ones(1,length(these_hitRates));
    glm_hitRate.kappa(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=ii_kappa_ordered(kappa_ii)*ones(1,length(these_hitRates));
    glm_hitRate_ii=glm_hitRate_ii+length(these_hitRates);

    id_hitRate_ii=id_hitRate_ii+1;
    input_hitRate_data(id_hitRate_ii).data=these_hitRates;
    input_hitRate_data(id_hitRate_ii).description=['Original k' num2str(kappa_ordered(kappa_ii))];


    glm_hitRate.data(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRate_sh))=these_hitRate_sh;
    glm_hitRate.shuffled(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRate_sh))=1*ones(1,length(these_hitRate_sh));
    glm_hitRate.kappa(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=ii_kappa_ordered(kappa_ii)*ones(1,length(these_hitRates));
    glm_hitRate_ii=glm_hitRate_ii+length(these_hitRate_sh);

    id_hitRate_ii=id_hitRate_ii+1;
    input_hitRate_data(id_hitRate_ii).data=these_hitRate_sh;
    input_hitRate_data(id_hitRate_ii).description=['Shuffled k' num2str(kappa_ordered(kappa_ii))];


end

plot(bar_offsets,these_means_hitRates,'-','Color',[0/255 158/255 115/255],'LineWidth',2,'MarkerFaceColor',[0/255 158/255 115/255])
plot(bar_offsets+1,these_means_hitRates_sh,'-','Color',[240/255 228/255 66/255],'LineWidth',2,'MarkerFaceColor',[240/255 228/255 66/255])

xlim([-0.5 bar_offset-0.5])
ylim([0 60])
xticks(bar_offsets+0.5);           % place ticks at 0, 2, 5, 10 [web:7]
xticklabels(string(kappa_values));  % use those numbers as labels [web:1]
  % tick labels from your cell array
xlabel('Kappa')
ylabel('Hit rate')
title('Temporal offset')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for temporal offset\n']);
fprintf(1, ['\nkappa values: 1:1, 2:0, 3:2, 4:5, 5:10\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for temporal offset\n']);
fprintf(fileID, ['\nkappa values: 1:1, 2:0, 3:2, 4:5, 5:10\n']);



tbl = table(glm_hitRate.data',glm_hitRate.shuffled',glm_hitRate.kappa',...
    'VariableNames',{'r2ES','shuffled','kappa'});
mdl = fitglm(tbl,'r2ES~shuffled+kappa'...
    ,'CategoricalVars',[2,3])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for temporal offset\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for temporal offset\n']);

 
[output_data] = drgMutiRanksumorTtest(input_hitRate_data, fileID,0);


%Now plot Hit rate for speed
speed_input_file='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Fabio_results_final_agent_navigation/Results_experiments/speed_results.mat';
load(speed_input_file)

all_hit_rates_speed=hit_rates_speed_all;

sh_speed_input_file='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Fabio_results_final_agent_navigation/Results_control/speed_results.mat';
load(sh_speed_input_file)

sh_all_hit_rates_speed=hit_rates_speed_all;

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[0:5:100];
rand_offset=0.5;

glm_hitRate=[];
glm_hitRate_ii=0;

id_hitRate_ii=0;
input_hitRate_data=[];

 
these_groups=[1 5]; %2 cm and 1 cm all trials


these_means_hitRates=[];
these_means_hitRates_sh=[];
bar_offsets=[];


for speed_ii=1:length(speeds_cm)
    %First get the hitRates for each number of ROIs
    these_hitRates=[];
    these_hitRates=all_hit_rates_speed(speed_ii,:);
    these_means_hitRates=[these_means_hitRates mean(these_hitRates)];
    

    %Then get the hitRates_sh for each number of ROIs
    these_hitRate_sh=[];
    these_hitRate_sh=sh_all_hit_rates_speed(speed_ii,:);
    these_means_hitRates_sh=[these_means_hitRates_sh mean(these_hitRate_sh)];
    
    bar_offsets=[bar_offsets bar_offset];

    %plot bar
    bar(bar_offset,mean(these_hitRates),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])

    %Violin plot
    try
        [mean_out, CIout,violin_x]=drgViolinPoint(these_hitRates...
            ,edges,bar_offset,rand_offset,'k','k',5);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;



    %plot bar
    bar(bar_offset,mean(these_hitRate_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])

    %Violin plot
    try
        [mean_out, CIout,violin_x]=drgViolinPoint(these_hitRate_sh...
            ,edges,bar_offset,rand_offset,'k','k',5);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    glm_hitRate.data(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=these_hitRates;
    glm_hitRate.shuffled(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=0*ones(1,length(these_hitRates));
    glm_hitRate.speed(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=speeds_cm(speed_ii)*ones(1,length(these_hitRates));
    glm_hitRate_ii=glm_hitRate_ii+length(these_hitRates);

    id_hitRate_ii=id_hitRate_ii+1;
    input_hitRate_data(id_hitRate_ii).data=these_hitRates;
    input_hitRate_data(id_hitRate_ii).description=['Original speed ' num2str(speeds_cm(speed_ii))];


    glm_hitRate.data(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRate_sh))=these_hitRate_sh;
    glm_hitRate.shuffled(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRate_sh))=1*ones(1,length(these_hitRate_sh));
    glm_hitRate.speed(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=speeds_cm(speed_ii)*ones(1,length(these_hitRates));
    glm_hitRate_ii=glm_hitRate_ii+length(these_hitRate_sh);

    id_hitRate_ii=id_hitRate_ii+1;
    input_hitRate_data(id_hitRate_ii).data=these_hitRate_sh;
    input_hitRate_data(id_hitRate_ii).description=['Shuffled speed ' num2str(speeds_cm(speed_ii))];


end

plot(bar_offsets,these_means_hitRates,'-','Color',[0/255 158/255 115/255],'LineWidth',2,'MarkerFaceColor',[0/255 158/255 115/255])
plot(bar_offsets+1,these_means_hitRates_sh,'-','Color',[240/255 228/255 66/255],'LineWidth',2,'MarkerFaceColor',[240/255 228/255 66/255])

xlim([-0.5 bar_offset-0.5])
ylim([0 80])
xticks(bar_offsets+0.5);           % place ticks at 0, 2, 5, 10 [web:7]
xticklabels(string(speeds_cm));  % use those numbers as labels [web:1]
  % tick labels from your cell array
xlabel('Speed')
ylabel('Hit rate')
title('Speed')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for speed\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for speed\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);



tbl = table(glm_hitRate.data',glm_hitRate.shuffled',glm_hitRate.speed',...
    'VariableNames',{'r2ES','shuffled','speed'});
mdl = fitglm(tbl,'r2ES~shuffled+speed+shuffled*speed'...
    ,'CategoricalVars',[2])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for speed\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for speed\n']);

 
[output_data] = drgMutiRanksumorTtest(input_hitRate_data, fileID,0);



%Now plot Hit rate for noise
noise_input_file='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Fabio_results_final_agent_navigation/Results_experiments/noise_results.mat';
load(noise_input_file)

all_hit_rates_noise=hit_rates_noise_all;

sh_noise_input_file='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Fabio_results_final_agent_navigation/Results_control/noise_results.mat';
load(sh_noise_input_file)

sh_all_hit_rates_noise=hit_rates_noise_all;

figureNo = figureNo + 1;
try
    close(figureNo)
catch
end
hFig=figure(figureNo);
hold on

ax=gca;ax.LineWidth=3;
set(hFig, 'units','normalized','position',[.2 .2 .3 .3])

bar_offset=0;

edges=[0:5:100];
rand_offset=0.5;

glm_hitRate=[];
glm_hitRate_ii=0;

id_hitRate_ii=0;
input_hitRate_data=[];

 
these_groups=[1 5]; %2 cm and 1 cm all trials


these_means_hitRates=[];
these_means_hitRates_sh=[];
bar_offsets=[];


for noise_ii=1:length(noise_levels)
    %First get the hitRates for each number of ROIs
    these_hitRates=[];
    these_hitRates=all_hit_rates_noise(noise_ii,:);
    these_means_hitRates=[these_means_hitRates mean(these_hitRates)];
    

    %Then get the hitRates_sh for each number of ROIs
    these_hitRate_sh=[];
    these_hitRate_sh=sh_all_hit_rates_noise(noise_ii,:);
    these_means_hitRates_sh=[these_means_hitRates_sh mean(these_hitRate_sh)];
    
    bar_offsets=[bar_offsets bar_offset];

    %plot bar
    bar(bar_offset,mean(these_hitRates),'LineWidth', 3,'EdgeColor','none','FaceColor',[0/255 158/255 115/255])

    %Violin plot
    try
        [mean_out, CIout,violin_x]=drgViolinPoint(these_hitRates...
            ,edges,bar_offset,rand_offset,'k','k',5);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;



    %plot bar
    bar(bar_offset,mean(these_hitRate_sh),'LineWidth', 3,'EdgeColor','none','FaceColor',[240/255 228/255 66/255])

    %Violin plot
    try
        [mean_out, CIout,violin_x]=drgViolinPoint(these_hitRate_sh...
            ,edges,bar_offset,rand_offset,'k','k',5);
    catch
        pffft=1;
    end
    bar_offset=bar_offset+1;

    glm_hitRate.data(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=these_hitRates;
    glm_hitRate.shuffled(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=0*ones(1,length(these_hitRates));
    glm_hitRate.noise(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=noise_levels(noise_ii)*ones(1,length(these_hitRates));
    glm_hitRate_ii=glm_hitRate_ii+length(these_hitRates);

    id_hitRate_ii=id_hitRate_ii+1;
    input_hitRate_data(id_hitRate_ii).data=these_hitRates;
    input_hitRate_data(id_hitRate_ii).description=['Original noise ' num2str(noise_levels(noise_ii))];


    glm_hitRate.data(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRate_sh))=these_hitRate_sh;
    glm_hitRate.shuffled(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRate_sh))=1*ones(1,length(these_hitRate_sh));
    glm_hitRate.noise(glm_hitRate_ii+1:glm_hitRate_ii+length(these_hitRates))=noise_levels(noise_ii)*ones(1,length(these_hitRates));
    glm_hitRate_ii=glm_hitRate_ii+length(these_hitRate_sh);

    id_hitRate_ii=id_hitRate_ii+1;
    input_hitRate_data(id_hitRate_ii).data=these_hitRate_sh;
    input_hitRate_data(id_hitRate_ii).description=['Shuffled noise ' num2str(noise_levels(noise_ii))];


end

plot(bar_offsets,these_means_hitRates,'-','Color',[0/255 158/255 115/255],'LineWidth',2,'MarkerFaceColor',[0/255 158/255 115/255])
plot(bar_offsets+1,these_means_hitRates_sh,'-','Color',[240/255 228/255 66/255],'LineWidth',2,'MarkerFaceColor',[240/255 228/255 66/255])

ylim([0 60])
xlim([-0.5 bar_offset-0.5])
xticks(bar_offsets+0.5);           % place ticks at 0, 2, 5, 10 [web:7]
xticklabels(string(noise_levels));  % use those numbers as labels [web:1]
  % tick labels from your cell array
xlabel('Noise')
ylabel('Hit rate')
title('Noise')

%Perform the glm  for prediction of odor all, hit, miss, shuffled
fprintf(1, ['\nglm for Fig ' num2str(figureNo) ' for noise\n']);
% fprintf(1, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);
fprintf(fileID, ['\nglm for Fig ' num2str(figureNo) ' for noise\n']);
% fprintf(fileID, ['\n1 Hit, 2 All, 3 Miss, 4 Shuffled\n']);



tbl = table(glm_hitRate.data',glm_hitRate.shuffled',glm_hitRate.noise',...
    'VariableNames',{'r2ES','shuffled','noise'});
mdl = fitglm(tbl,'r2ES~shuffled+noise+shuffled*noise'...
    ,'CategoricalVars',[2])

txt = evalc('mdl');
txt=regexp(txt,'<strong>','split');
txt=cell2mat(txt);
txt=regexp(txt,'</strong>','split');
txt=cell2mat(txt);

fprintf(fileID,'%s\n', txt);


%Do the ranksum/t-test
fprintf(1, ['\n\nRanksum or t-test p values for noise\n'])
fprintf(fileID, ['\n\nRanksum or t-test p values for noise\n']);

 
[output_data] = drgMutiRanksumorTtest(input_hitRate_data, fileID,0);

fclose(fileID);

pffft=1;