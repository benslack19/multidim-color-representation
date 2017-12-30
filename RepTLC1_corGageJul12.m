% load BV_waveplay.mat
% save BV_waveplay.mat

% Purpose: Color code each frame (ex. 1st frame is blue, 2nd is red, etc.
% of a movie and merge them to create a 2D wave map.

% The only thing that should be imported is the raw time sequence - convert
% to 8-bit tif, and save in workspace


% Input file, # frames

[x1,y1] = size(stack);  % Input file name without .tif
frame_range = 2:12;   % Input frame range to look at
t=length(frame_range);
AugStrike = uint8(zeros(x1,y1,1,length(frame_range))); % preallocate 4-D array; 2nd and 3rd values are picture resolution

for frame=1:t;
	[AugStrike(:,:,:,frame),map] = imread('stack.tif', frame_range(frame));  % Input file name with .tif Note: files must be in workspace
end

AugStrike1 = im2double(AugStrike(:,:,:,:)); % Changed variable type from unit8 to double of my original time series

% get frame range for average
frame_range2 = frame_range;  
frame_range2 = 1:3;
t2=length(frame_range2);
AugStrike5 = uint8(zeros(x1,y1,1,length(frame_range2))); % preallocate 4-D array; 2nd and 3rd values are picture resolution

for frame=1:t2;
	[AugStrike5(:,:,:,frame),map] = imread('stack.tif', (frame + max(frame_range2)-t2));  % Input file name with .tif Note: files must be in workspace
end

AugStrike2 = im2double(AugStrike5(:,:,:,:)); % Changed variable type from unit8 to double of my original time series

%% Subtract average before assigning color
% Generate a composite (addition of all pixels)

c=1:t2; % determine frames for baseline/average, choose t for all frames
d=length(c);
a=1; b=d;
x=0;
AugStrikePreMerge = zeros(x1,y1);
while b >= a
    x = AugStrike2(:,:,:,c(a));
    AugStrikePreMerge = x + AugStrikePreMerge;
    a=a+1;
end

% figure, imshow(AugStrikePreMerge)    % The composite image without
% baseline subtraction

% Get the average
AvgFact = 1;
AvgPreMerge = AvgFact*(AugStrikePreMerge/d);    % The averaged-across-time image/values

% Try standardizing?    
% figure              % Compare composite and average image
% subplot(1,2,1), imshow(AugStrikePreMerge)  
% subplot(1,2,2), imshow(AvgPreMerge)

    % Then subtract the AvgPreMerge from every image in the time series
AugStrikePre = zeros(x1,y1,1,t);
for tq=1:t;
AugStrikePre(:,:,:,tq) = AugStrike1(:,:,:,tq)-AvgPreMerge;    % I get "Warning: Input arguments must be scalar."
end

%% Assign color

cmap = colormap(jet);
mycolormap = cmap;
tq=1:t;
colorValues = linspace(min(tq),max(tq),length(mycolormap(:,1)));

AugStrikePreEd = zeros(x1,y1,3,t);
for tq = 1:t;                                   % 6 for the number of images in the time series
    q = find((abs(colorValues - tq)) == min(abs(colorValues - tq)));   % I'm assigning one color for each timepoint
    if length(q) == 1;    
    k = mycolormap(q,:);
    else length(q) > 1;   %  This is for the midpoint of the jet colormap
    j = mycolormap(q,:);    
    k = mean(j); %
    end %  This is the color assignment from the jet colormap
        for z = 1:3;
        AugStrikePreEd(:,:,z,tq) = (AugStrikePre(:,:,1,tq))*k(z); % I multiplied the one starting intensity value by the 3 values of the assigned color       
        end
end

%{
figure        % I'm showing each image of the time series withOUT its assigned color
for t = 1:size(AugStrikePre(:,:,:),3);
subplot(2,3,t), imshow(AugStrikePre(:,:,:,t))      
end

figure        % I'm showing each image of the time series with its assigned color
for t = 1:size(AugStrikePreEd(:,:,:,:),4);
subplot(1,size(AugStrikePreEd(:,:,:,:),4),t), imshow(AugStrikePreEd(:,:,:,t))      
end
%}

%% THIS COMPOSITE IMAGE WITH BASELINE SUBTRACTION LOOKS GOOD 

% Generate a composite image, average subtracted from each frame

c=1:t;
a=1; b=t;
x=0;
AugStrikeMergePost = zeros(x1,y1,3);
while b >= a
    x = AugStrikePreEd(:,:,:,c(a));
    AugStrikeMergePost = x + AugStrikeMergePost;
    a=a+1;
end

figure
imshow(AugStrikeMergePost)

% imwrite(AugStrikeMergePost, 'TLC1.tif','tif')  % save tif file to workspace folder

%}



%% SHOW EACH IMAGE WITH ASSIGNED COLOR WITHOUT BASELINE SUBTRACTION
% LOOKS BETTER

%assign color

cmap = colormap(jet);
mycolormap = cmap;
tq=1:t;
colorValues = linspace(min(tq),max(tq),length(mycolormap(:,1)));

AugStrikePreEd = zeros(x1,y1,3,t);
for tq = 1:t;                                   % 6 for the number of images in the time series
    q = find((abs(colorValues - tq)) == min(abs(colorValues - tq)));   % I'm assigning one color for each timepoint
    if length(q) == 1;    
    k = mycolormap(q,:);
    else length(q) > 1;   %  This is for the midpoint of the jet colormap
    j = mycolormap(q,:);    
    k = mean(j); %
    end %  This is the color assignment from the jet colormap
        for z = 1:3;
        AugStrikePreEd(:,:,z,tq) = (AugStrike1(:,:,1,tq))*k(z); % I multiplied the one starting intensity value by the 3 values of the assigned color       
        end
end



 % Composite without subtracting average, usually ugly
c=1:t;
a=1; b=t;
x=0;
AugStrikeMergePostCol = zeros(x1,y1,3);
while b >= a
    x = AugStrikePreEd(:,:,:,c(a));
    AugStrikeMergePostCol = x + AugStrikeMergePostCol;
    a=a+1;
end

figure
imshow(AugStrikeMergePostCol)





figure        % I'm showing each image of the time series with its assigned color
for ts = 1:t;
subplot(t/3,t/3,ts), imshow(AugStrikePreEd(:,:,:,ts))      
end


%% Make a movie
%{
Exporting Audio/Video DataMATLAB includes several functions that you can use to export audio or
video data from the MATLAB workspace. These functions write audio data to
a file using specific file formats. The following sections describe

%}

mov = immovie(AugStrikePreEd);
movie2avi(mov,'TLC.avi','compression', 'none')

aviobj = avifile('AugStrikePreEd','fps',1)

for k=1:16
    h = plot(fft(eye(k+16)));
    set(h,'EraseMode','xor');
    axis equal;
    frame = getframe(gca);
    aviobj = addframe(aviobj,frame);
end

%}