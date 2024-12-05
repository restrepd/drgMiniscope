% drgMini_read_undistortion
%This code reshapes the output of the mp4 file to a rectangular odor arena
% The file is structured like this:
% orig_x,orig_y,corr_x,corr_y
clear all
close all

%These are the files to be processed, we need the mp4 and the arena files
%and the crop values from the yaml metadata file
%MATLAB is not able to read the mp4 file. Use Flie convert/stream in VLC,
%it saves with m4v suffix, change the suffix to mp4

%mp4 path
this_mp4_path='/data2/SFTP/PreProcessed/20221028_FCM22_1/';

% mp4_file='20220713_arena01_session001_topCam-0000conv.mp4' %Note: I had to open with VLC and file convert/stream saved as mp4
mp4_file='20221028_arena01_session004_topCam-0000_conv.mp4';

this_arena_path=this_mp4_path;

% arena_file='20220713_FCM6withodor_odorarena_L1andL4_syn.mat';
arena_file='20221028_FCM22_1withodor_odorarena_L1andL4_sync.mat';

%Enter the crop values from the metadata yaml file
%These are the first and thrid crop values in the file
%e.g. for file 20220804_arena01_session001_metadata.yaml
% Crop:
% - - 67
%   - 280
%   - 3
%   - 250
%
% cropx=67;
% cropy=3;

cropx=67;
cropy=3;

%If the mouse is at the corner obstructing the click change this to a
%different value
numImagesToRead=100;

%Odor arena size
y_length=480; %mm
x_length=500; %mm

is_sphgpu=1;

if is_sphgpu==1
    this_path='/data/SFTP/PreProcessedDR/Undistortion_odorarena/';
else
    this_path='/Users/restrepd/Documents/Projects/SFTP/Fabio_OdorArena_GoodData/Undistortion_odorarena/';
end
correction_filename='coordinate_corrections.txt';
correction_table=readtable([this_path correction_filename]);





%Plot the edges of the camera imaged odor arena
x_orig=correction_table.Var1;
y_orig=correction_table.Var2;

original_corners=[60 65; 935 67; 925 1029; 70 1031];


border_x=[];
border_y=[];
figure(1)
hold on
these_x_orig=unique(x_orig);
for ii=1:length(these_x_orig)
    ii_x=these_x_orig(ii);
    [this_y_orig_min min_ii]=min(y_orig(x_orig==ii_x));
    border_x=[border_x ii_x];
    border_y=[border_y this_y_orig_min];
    plot(ii_x,this_y_orig_min,'.')
    [this_y_orig_max max_ii]=max(y_orig(x_orig==ii_x));
    border_x=[border_x ii_x];
    border_y=[border_y this_y_orig_max];
    plot(ii_x,this_y_orig_max,'.')
end

these_y_orig=unique(y_orig);
for ii=1:length(these_y_orig)
    ii_y=these_y_orig(ii);
    [this_x_orig_min min_ii]=min(x_orig(y_orig==ii_y));
    border_x=[border_x this_x_orig_min];
    border_y=[border_y ii_y];
    plot(this_x_orig_min,ii_y,'.')
    [this_x_orig_max max_ii]=max(x_orig(y_orig==ii_y));
    border_x=[border_x this_x_orig_max];
    border_y=[border_y ii_y];
    plot(this_x_orig_max,ii_y,'.')
end

for ii=1:size(original_corners,1)
    plot(original_corners(ii,1),original_corners(ii,2),'ok')
end

set(gca, 'YDir', 'reverse');
title('Odor arena boundaries, camera view')
xlabel('x')
ylabel('y')



%Now plot the edges of the physical coordinates of the odor arena
x_corr=correction_table.Var3;
x_corr=x_length*(x_corr-min(x_corr))/(max(x_corr)-min(x_corr));

y_corr=correction_table.Var4;
y_corr=y_length*(y_corr-min(y_corr))/(max(y_corr)-min(y_corr));

figure(2)
hold on
these_x_corr=unique(x_corr);
for ii=1:length(these_x_corr)
    ii_x=these_x_corr(ii);
    [this_y_corr_min min_ii]=min(y_corr(x_corr==ii_x));
    plot(ii_x,this_y_corr_min,'.')
    [this_y_corr_max max_ii]=max(y_corr(x_corr==ii_x));
    plot(ii_x,this_y_corr_max,'.')
end

these_y_corr=unique(y_corr);
for ii=1:length(these_y_corr)
    ii_y=these_y_corr(ii);
    [this_x_corr_min min_ii]=min(x_corr(y_corr==ii_y));
    plot(this_x_corr_min,ii_y,'.')
    [this_x_corr_max max_ii]=max(x_corr(y_corr==ii_y));
    plot(this_x_corr_max,ii_y,'.')
end

set(gca, 'YDir', 'reverse');
title('Odor arena boundaries, flat rendering')
ylim([-50 530])
xlim([-50 550])
xlabel('x')
ylabel('y')

%Read a frame of the mp4 file
figure(3)
vidObj.CurrentTime = 1;
vidObj = VideoReader([this_mp4_path mp4_file]);
for iif=1:numImagesToRead
    vidFrame = readFrame(vidObj);
end
imshow(vidFrame)
hFig = figure(3);
set(hFig, 'units','normalized','position',[.1 .1 .8 .8])

title('Click on the four corners of the arena (start top left, clockwise)');
hold on

% Step 2: Use ginput to capture four points
[x_in, y_in] = ginput(4);

% Step 3: Plot the captured points
plot(x_in, y_in, 'ro', 'MarkerSize', 10, 'LineWidth', 2);

% hold on
% these_x_orig=unique(x_orig);
% for ii=1:length(these_x_orig)
%     ii_x=these_x_orig(ii);
%     [this_y_orig_min min_ii]=min(y_orig(x_orig==ii_x));
%     plot(ii_x,this_y_orig_min,'.')
%     [this_y_orig_max max_ii]=max(y_orig(x_orig==ii_x));
%     plot(ii_x,this_y_orig_max,'.')
% end
% 
% these_y_orig=unique(y_orig);
% for ii=1:length(these_y_orig)
%     ii_y=these_y_orig(ii);
%     [this_x_orig_min min_ii]=min(x_orig(y_orig==ii_y));
%     plot(this_x_orig_min,ii_y,'.')
%     [this_x_orig_max max_ii]=max(x_orig(y_orig==ii_y));
%     plot(this_x_orig_max,ii_y,'.')
% end




%Save the frame as a tiff file
imwrite(vidFrame, [this_mp4_path mp4_file(1:end-4) '.tiff']);

%Now downsample the odor arena
x_orig_dn=[];
y_orig_dn=[];
x_corr_dn=[];
y_corr_dn=[];
ii_dn=0;
delta_length=2;
for y=0:delta_length:y_length
    for x=0:delta_length:x_length
        [this_sq_dis, this_min_ii]=min((x_corr-x).^2+(y_corr-y).^2);
        ii_dn=ii_dn+1;
        ii_final=this_min_ii(1);
        x_corr_dn(ii_dn)=x;
        y_corr_dn(ii_dn)=y;
        x_orig_dn(ii_dn)=x_orig(ii_final);
        y_orig_dn(ii_dn)=y_orig(ii_final);
    end
end

%plot downsampled camera view coordinates
figure(4)
hold on
plot(x_orig_dn,y_orig_dn,'.')
set(gca, 'YDir', 'reverse');
title('Odor arena coordinates, camera view')
xlabel('x')
ylabel('y')



%Now plot downsampled physical coordinates of the odor arena
figure(5)
hold on
plot(x_corr_dn,y_corr_dn,'.')

set(gca, 'YDir', 'reverse');
title('Odor arena coordinates, flat rendering')
ylim([-50 530])
xlim([-50 550])
xlabel('x')
ylabel('y')

%Now perform an isometric shift of the pattern of x_orig_dn y_orig_dn
%Code was from perplexity



% %Enter corners of the odor arena mesured in the movie
% x_top_left=77;
% y_top_left=58;
% 
% x_top_right=1039;
% y_top_right=65;
% 
% x_bottom_right=1074;
% y_bottom_right=970;
% 
% x_bottom_left=56;
% y_bottom_left=975;

%Corners of the odor arena measured by the user
x_top_left=x_in(1);
y_top_left=y_in(1);

x_top_right=x_in(2);
y_top_right=y_in(2);

x_bottom_right=x_in(3);
y_bottom_right=y_in(3);

x_bottom_left=x_in(4);
y_bottom_left=y_in(4);




%Here are the corners of the odor arena for this session
these_corners = [x_top_left y_top_left;...
    x_top_right y_top_right;...
    x_bottom_right y_bottom_right;...
    x_bottom_left y_bottom_left];

% Calculate the transformation matrix
T = fitgeotform2d(original_corners, these_corners, 'projective');

% Define your original points
original_points=[x_orig_dn' y_orig_dn']; 

% Apply the transformation to your points
these_points = transformPointsForward(T, original_points);

x_orig_dn_sh=these_points(:,1)';
y_orig_dn_sh=these_points(:,2)';


%plot downsampled shifted camera view coordinates
figure(6)
hold on
plot(x_orig_dn_sh,y_orig_dn_sh,'.')
set(gca, 'YDir', 'reverse');
title('Odor arena coordinates, camera view, shifted')
xlabel('x')
ylabel('y')

figure(7)
vidObj.CurrentTime = 1;
vidObj = VideoReader([this_mp4_path mp4_file]);
vidFrame = readFrame(vidObj);
imshow(vidFrame)

hold on
plot(x_orig_dn_sh,y_orig_dn_sh,'.y')
set(gca, 'YDir', 'reverse');
title('Odor arena coordinates, camera view, shifted')
xlabel('x')
ylabel('y')

figure(8)

imshow(vidFrame)

hold on
this_delta=2*delta_length;

min_x_orig_dn_sh=min(x_orig_dn_sh);
max_x_orig_dn_sh=max(x_orig_dn_sh);

for x=min_x_orig_dn_sh-this_delta/2:this_delta:max_x_orig_dn_sh+this_delta/2
    [this_y_orig_dn_sh_min min_ii]=min(y_orig_dn_sh((x_orig_dn_sh>=x)&(x_orig_dn_sh<x+this_delta)));
    if ~isempty(this_y_orig_dn_sh_min)
        plot(x,this_y_orig_dn_sh_min,'.')
    end
    [this_y_orig_dn_sh_max min_ii]=max(y_orig_dn_sh((x_orig_dn_sh>=x)&(x_orig_dn_sh<x+this_delta)));
    if ~isempty(this_y_orig_dn_sh_max)
        plot(x,this_y_orig_dn_sh_max,'.')
    end
end

min_y_orig_dn_sh=min(y_orig_dn_sh);
max_y_orig_dn_sh=max(y_orig_dn_sh);

for y=min_y_orig_dn_sh-this_delta/2:this_delta:max_y_orig_dn_sh+this_delta/2
    [this_x_orig_dn_sh_min min_ii]=min(x_orig_dn_sh((y_orig_dn_sh>=y)&(y_orig_dn_sh<y+this_delta)));
    if ~isempty(this_x_orig_dn_sh_min)
        plot(this_x_orig_dn_sh_min, y,'.')
    end
    [this_x_orig_dn_sh_max min_ii]=max(x_orig_dn_sh((y_orig_dn_sh>=y)&(y_orig_dn_sh<y+this_delta)));
    if ~isempty(this_x_orig_dn_sh_max)
        plot(this_x_orig_dn_sh_max,y,'.')
    end
end


title('Odor arena outline, camera view, shifted')
xlabel('x')
ylabel('y')

%Now translate the arena file x and y coordinates to mm
load([this_arena_path arena_file])
arena.xsync_orig=4*(arena.xsync-cropx);
arena.ysync_orig=4*(arena.ysync-cropy);
 
figure(9)
hold on

plot(arena.xsync(end),arena.ysync(end),'ok')
% plot(arena.xsync(end-these_no_points:end),arena.ysync(end-these_no_points:end),'-b')
plot(arena.xsync_orig,arena.ysync_orig,'-b')
set(gca, 'YDir', 'reverse');
% title(['End of trajectory for ' arena_file])
title(['Trajectory before morphing '])
ylabel('y')
xlabel('x')
legend('end of trajectory','trajectory')

for ii=1:length(arena.xsync)
    this_x=arena.xsync_orig(ii);
    this_y=arena.ysync_orig(ii);
    sq_distances=(y_orig_dn_sh-this_y).^2 + (x_orig_dn_sh-this_x).^2;
    [this_min_sq_dist this_min_ii]=min(sq_distances);
    arena.xsync(ii)=x_corr_dn(this_min_ii);
    arena.ysync(ii)=y_corr_dn(this_min_ii);
end
 

figure(10)
hold on

plot(arena.xsync(end),arena.ysync(end),'ok')
% plot(arena.xsync(end-these_no_points:end),arena.ysync(end-these_no_points:end),'-b')
plot(arena.xsync,arena.ysync,'-b')
set(gca, 'YDir', 'reverse');
% title(['End of trajectory for ' arena_file])
title(['Trajectory after morphing '])
ylabel('y')
xlabel('x')
legend('end of trajectory','trajectory')

save([this_arena_path arena_file(1:end-4) '_mm.mat'],'arena','-v7.3')
pffft=1;

