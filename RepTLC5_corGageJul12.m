% load RepTLC1_space.mat
% save RepTLC1_space.mat

%% Try to take the frame with the biggest intensity that also exceeds SD

% Purpose: Color code each frame (ex. 1st frame is blue, 2nd is red, etc.
% of a movie and merge them to create a 2D wave map.

% This will try to use the standard deviation as a criterion for assigning
% color.  If a pixel intensity goes above one SD, then assign a color at
% the first timepoint that it goes above one SD.  If one pixel is less than
% 1 SD during the whole time-series, it will be assigned a zero.

% The code will also be simplified by removing some comments, etc.

% The only thing that should be imported is the raw time sequence - convert
% to 8-bit tif, and save in workspace

% Input file, # frames

[x1,y1] = size(stack);  % Input file name without .tif
frame_range = 240:2:290;   % Input frame range to look at
t=length(frame_range);
MultDArray8 = uint8(zeros(x1,y1,1,length(frame_range))); % preallocate 4-D array; 2nd and 3rd values are picture resolution

for frame=1:t;
	[MultDArray8(:,:,:,frame),map] = imread('stack.tif', frame_range(frame));  % Input file name with .tif Note: files must be in workspace
end

MultDArrayDbl = im2double(MultDArray8(:,:,:,:)); % Changed variable type from unit8 to double of my original time series

frame_range2 = 1:1;   % Input frame range to look at
ts=length(frame_range2);
MultDArray8 = uint8(zeros(x1,y1,1,length(frame_range2))); % preallocate 4-D array; 2nd and 3rd values are picture resolution

% get frame range for baseline
for frame=1:ts;
	[MultDArray8(:,:,:,frame),map] = imread('stack.tif', frame_range2(frame));  % Input file name with .tif Note: files must be in workspace
end

MultDArrayDbl2 = im2double(MultDArray8(:,:,:,:)); % Changed variable type from unit8 to double of my original time series

%% Determine average and SD
% Generate a composite (addition of all pixels)

c=1:ts;     % determine frames for baseline/average, choose t for all frames
AvgFact = 1;
SDlevel = 1;
Thresh = 1.2;
d=length(c);
a=1; b=d;
x=0;
MultDArrayComp = zeros(x1,y1);
while b >= a
    x = MultDArrayDbl2(:,:,:,c(a));
    MultDArrayComp = x + MultDArrayComp;
    a=a+1;
end

MultDArrayMod = MultDArrayDbl2(:,:,:,1:d);  % MD matrix for average

% figure, imshow(MultDArrayComp)    % The composite image without
% baseline subtraction, usually ugly

% Get the average - old way
ArrayAvg = (AvgFact) * (MultDArrayComp/d);    % The averaged-across-time image/values

% Get the average, all timeframes - mean function, new way, much easier
% AvgPreMerge2 = mean(MultDArrayDbl,4);

% Get the SD the SD of time-constrained matrix
ArrayModSD = SDlevel * std(MultDArrayMod,1,4);

% Get the SD - much easier way
% MultDArraySD = SDlevel*std(MultDArrayDbl,1,4);
% SD, dividing by n, not n-1 THIS WORKS, MUCH SIMPLER, GET DIMENSION I WANT

% Generate a new array of the threshold, defined as avg+SD
% ArrayThresh = ArrayAvg + ArrayModSD;
ArrayThresh = (Thresh) * (MultDArrayComp/d) + 0.2;  % this way just set an average of 20% above

% determine which pixels exceed threshold
MultDArrayFindSDs = zeros(size(MultDArrayDbl,1),size(MultDArrayDbl,2),1,t);
for ts=1:t;
    MultDArrayHold = zeros(size(MultDArrayDbl,1),size(MultDArrayDbl,2));
    k = find(MultDArrayDbl(:,:,1,ts) > ArrayThresh);
    MultDArrayHold(k) = 1;
    MultDArrayFindSDs(:,:,:,ts) = MultDArrayHold;
end

%% find timeframes for each pixel where it exceeds threshold, save in MultDArrayThresh

MultDArrayThresh = zeros(size(MultDArrayDbl,1),size(MultDArrayDbl,2),1,t);
for ts=1:t;
    MultDArrayHold = MultDArrayDbl(:,:,:,ts);
    k = find(MultDArrayDbl(:,:,1,ts) < ArrayThresh);    % what it was: (not equal) k = find(MultDArrayDbl(:,:,1,ts) <= ArrayThresh);
    MultDArrayHold(k) = 0;
    MultDArrayThresh(:,:,:,ts) = MultDArrayHold;
end

%{
figure        % I'm showing each image of the time series after determining threshold
for ts = 1:t;
subplot(t/3,t/3,ts), imshow(MultDArrayThresh(:,:,:,ts))      
end
%}

%% find the first frame that exceeds threshold, and assign intensity,
% otherwise set as zero

MultDArrayFirstF = zeros(size(MultDArrayDbl,1),size(MultDArrayDbl,2),1,t);
for a=1:x1
    for b=1:y1
    k=find(MultDArrayFindSDs(a,b,1,:)==1,1,'first');     
    MultDArrayFirstF(a,b,1,k) = MultDArrayDbl(a,b,1,k);
    end
end

%{
figure        % I'm showing each image of the time series after determining threshold
for ts = 1:t;
subplot(t/3,t/3,ts), imshow(MultDArrayFirstF(:,:,:,ts))      
end
%}

%% find the frame that exceeds threshold, and has the biggest intensity,
% otherwise set as zero
MultDArrayMaxF = zeros(size(MultDArrayDbl,1),size(MultDArrayDbl,2),1,t);
for a=1:x1
    for b=1:y1
    [c,i]=max(MultDArrayDbl(a,b,1,:));  % find the frame of the max intensity
    k=find(MultDArrayFindSDs(a,b,1,:)== 1); % find the frames which are above threshold
    l=find(k==i); % find where both conditions are met
    m=k(l);
    MultDArrayMaxF(a,b,1,m) = MultDArrayDbl(a,b,1,m);
    end
end

%{
figure        % I'm showing each image of the time series after determining threshold
for ts = 1:t;
subplot(t/3,t/3,ts), imshow(MultDArrayMaxF(:,:,:,ts))      
end
%}
%% choose method before assigning color
%MultDArrayPre1 = zeros(x1,y1,1,t);
MultDArrayPre2 = zeros(x1,y1,1,t);
MultDArrayPre3 = zeros(x1,y1,1,t);
%MultDArrayPre1 = MultDArrayThresh;
MultDArrayPre2 = MultDArrayFirstF;
MultDArrayPre3 = MultDArrayMaxF;
%% subtract average from SD method, usually bad
%{
MultDArrayPre = zeros(x1,y1,1,t);
for tq=1:t;
MultDArrayPre(:,:,:,tq) = MultDArrayMaxF(:,:,:,tq)-AvgPreMerge;
end
%}
%% Assign color

cmap = colormap(jet);
mycolormap = cmap;
tq=1:t;
colorValues = linspace(min(tq),max(tq),length(mycolormap(:,1)));
%{
MultDArrayColor1 = zeros(x1,y1,3,t);
for tq = 1:t;                                   % 6 for the number of images in the time series
    q = find((abs(colorValues - tq)) == min(abs(colorValues - tq)));   % I'm assigning one color for each timepoint
    if length(q) == 1;    
    k = mycolormap(q,:);
    else length(q) > 1;   %  This is for the midpoint of the jet colormap
    j = mycolormap(q,:);    
    k = mean(j); %
    end %  This is the color assignment from the jet colormap
        for z = 1:3;
        MultDArrayColor1(:,:,z,tq) = (MultDArrayPre1(:,:,1,tq))*k(z); % I multiplied the one starting intensity value by the 3 values of the assigned color       
        end
end

%}

MultDArrayColor2 = zeros(x1,y1,3,t);
for tq = 1:t;                                   % 6 for the number of images in the time series
    q = find((abs(colorValues - tq)) == min(abs(colorValues - tq)));   % I'm assigning one color for each timepoint
    if length(q) == 1;    
    k = mycolormap(q,:);
    else length(q) > 1;   %  This is for the midpoint of the jet colormap
    j = mycolormap(q,:);    
    k = mean(j); %
    end %  This is the color assignment from the jet colormap                            %  This is the color assignment from the jet colormap
        for z = 1:3;
        MultDArrayColor2(:,:,z,tq) = (MultDArrayPre2(:,:,1,tq))*k(z); % I multiplied the one starting intensity value by the 3 values of the assigned color       
        end
end

MultDArrayColor3 = zeros(x1,y1,3,t);
for tq = 1:t;                                   % 6 for the number of images in the time series
    q = find((abs(colorValues - tq)) == min(abs(colorValues - tq)));   % I'm assigning one color for each timepoint
    if length(q) == 1;    
    k = mycolormap(q,:);
    else length(q) > 1;   %  This is for the midpoint of the jet colormap
    j = mycolormap(q,:);    
    k = mean(j); %
    end %  This is the color assignment from the jet colormap                            %  This is the color assignment from the jet colormap
        for z = 1:3;
        MultDArrayColor3(:,:,z,tq) = (MultDArrayPre3(:,:,1,tq))*k(z); % I multiplied the one starting intensity value by the 3 values of the assigned color       
        end
end

%{
figure        % I'm showing each image of the time series withOUT its assigned color
for t = 1:size(MultDArrayPre(:,:,:),3);
subplot(2,3,t), imshow(MultDArrayPre(:,:,:,t))      
end

figure        % I'm showing each image of the time series with its assigned color
for t = 1:size(MultDArrayColor2(:,:,:,:),4);
subplot(1,size(MultDArrayColor2(:,:,:,:),4),t), imshow(MultDArrayColor2(:,:,:,t))      
end
%}

%% THIS COMPOSITE IMAGE WITH BASELINE SUBTRACTION LOOKS GOOD 

%{
% Generate a composite image, average subtracted from each frame
c=1:t;
a=1; b=t;
x=0;
ArrayTLC1 = zeros(x1,y1,3);
while b >= a
    x = MultDArrayColor1(:,:,:,c(a));
    ArrayTLC1 = x + ArrayTLC1;
    a=a+1;
end

%}

c=1:t;
a=1; b=t;
x=0;
ArrayTLC2 = zeros(x1,y1,3);
while b >= a
    x = MultDArrayColor2(:,:,:,c(a));
    ArrayTLC2 = x + ArrayTLC2;
    a=a+1;
end

c=1:t;
a=1; b=t;
x=0;
ArrayTLC3 = zeros(x1,y1,3);
while b >= a
    x = MultDArrayColor3(:,:,:,c(a));
    ArrayTLC3 = x + ArrayTLC3;
    a=a+1;
end

figure
%subplot (2,2,1),imshow(ArrayTLC1)
subplot (1,2,1),imshow(ArrayTLC2)
subplot (1,2,2),imshow(ArrayTLC3)

% figure, imshow(ArrayTLC2)
% figure, imshow(ArrayTLC3)


%% SHOW EACH IMAGE WITH ASSIGNED COLOR WITHOUT BASELINE SUBTRACTION - NOTE
%% SEE TLC1 TO FIX THIS (BEN, 10/8/09)
MultDArrayMov = zeros(x1,y1,3,t);   
% Preallocated a new 4D matrix-I have 6 images in the time series..

for tq = 1:t;                                   % 6 for the number of images in the time series
    q = find((abs(colorValues - tq)) == min(abs(colorValues - tq)));   % I'm assigning one color for each timepoint
    if length(q) == 1;    
    k = mycolormap(q,:);
    else length(q) > 1;   %  This is for the midpoint of the jet colormap
    j = mycolormap(q,:);    
    k = mean(j); %
    end %  This is the color assignment from the jet colormap
        for z = 1:3;
        MultDArrayMov(:,:,z,tq) = (MultDArrayThresh(:,:,1,tq))*k(z); % I multiplied the one starting intensity value by the 3 values of the assigned color       
        end
end


% Make a movie
mov = immovie(MultDArrayMov);
movie2avi(mov,'TLCmovie.avi','compression', 'none')

aviobj = avifile('MultDArrayMov','fps',1);

%{
for k=1:9
    h = plot(fft(eye(k+16)));
    set(h,'EraseMode','xor');
    axis equal;
    frame = getframe(gca);
    aviobj = addframe(aviobj,frame);
end
%}

