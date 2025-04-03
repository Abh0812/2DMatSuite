%% Hall Bar Device Generation Using Markers
% This script loads an image of an SOI chip with pre‐printed alignment markers.
% The markers are used to define the chip’s coordinate system. The user then:
% 1) clicks on three marker centers (e.g., top‐right, bottom‐right, bottom‐left),
%    and enters the marker letters.
% 2) outlines the flake (by drawing a polygon).
% 3) clicks on two points along the flake edge to define the current channel.
% 4) The script generates two current electrodes (at the channel ends) and four 
%    voltage electrodes (placed perpendicular to the channel at 1/4 and 3/4 positions).
% 5) Finally, all electrode coordinates (translated using the markers) are saved 
%    to a text file for use in GDS layout.
%
% Note: The marker positions and letter-to-value mapping are assumed to be as in your
% original code.

close all; clear; clc;

%% Marker & Pad Information (from original code)
pads_centers = [ -2750,  500;
                 -3500, 2200;
                 -2500, 3500;
                 -500,  3500;
                  500,  3500;
                 2500, 3500;
                 3500, 2200;
                 2750,  500;
                 2750, -500;
                 3500, -2200;
                 2500, -3500;
                  500, -3500;
                 -500, -3500;
                 -2500,-3500;
                 -3500,-2200;
                 -2750, -500];
pad_size = 500;
lettervals = dictionary(["A" "B" "C" "D" "E" "F" "G" "H" "I" "J"], [0 1 2 3 4 5 6 7 8 9]);

%% Initialize containers
cornerpositions = {};
nn = 0;
electrodes = {};

%% Step 1: Load Image & Get Marker Positions
% Open image file
[filename, folder] = uigetfile('*.*', 'Select chip image');
if isequal(filename,0)
    error('No image selected.');
end
[~, baseFilename, ~] = fileparts(filename);
addpath(folder);
ChipImage = imread(fullfile(folder, filename));
figure; imshow(ChipImage);
title('Zoom in if necessary, then press Enter');
zoom on;
pause;  % Allow zooming

% Prompt user to click on three alignment markers.
% For example: (1) TOP RIGHT, (2) BOTTOM RIGHT, (3) BOTTOM LEFT.
title({'Click on the (1) TOP RIGHT cross center,', 'then the (2) BOTTOM RIGHT cross center,', 'then the (3) BOTTOM LEFT cross center.'});
disp('Click on the (1) TOP RIGHT, (2) BOTTOM RIGHT, and (3) BOTTOM LEFT marker centers.');
numMarkers = 3;
xmark = zeros(numMarkers,1); ymark = zeros(numMarkers,1);
for i = 1:numMarkers
    [xmark(i), ymark(i)] = ginput(1);
    plot(xmark(i), ymark(i), 'r+', 'MarkerSize', 10, 'LineWidth', 2);
    disp(['Marker ' num2str(i) ': (' num2str(xmark(i)) ', ' num2str(ymark(i)) ')']);
end

% Identify the field using marker letters.
markNum = baseFilename; % default is the file name, if desired
userQ = ['Letters are currently "' markNum '", type in the letters of the last marker you clicked:'];
changemarkNum = inputdlg(userQ, 'Input Dialog', [1 50], {markNum});
if ischar(changemarkNum{1}) && length(changemarkNum{1}) == 4
    markNum = changemarkNum{1};
end

% Define translation: use the third marker (BOTTOM LEFT) as the origin.
Translation = [xmark(3); ymark(3)];
% For consistency with your original code:
xmark(1) = xmark(3);
ymark(2) = ymark(3);
% Translate markers so that origin is at Marker 3.
xmark = xmark - Translation(1);
ymark = ymark - Translation(2);
% Invert y-axis (as done originally)
ymark = ymark * (-1);

%% Step 2: Get Flake Outline Using Polygon (Same as Original)
% Reset axes limits
xlim([1 size(ChipImage, 2)]);
ylim([1 size(ChipImage, 1)]);
title('Click on the corners of the nanostructure (flake) with the mouse');
hPolygon = drawpolygon('LineWidth', 1, 'Color', 'g');  % Draw flake boundary
vertices = hPolygon.Position;   % Nx2 matrix in image coordinates
og_cornerdata = vertices;
plot(og_cornerdata(:,1), og_cornerdata(:,2),'r+','MarkerSize', 20);
for i = 1:size(og_cornerdata,1)
    text(og_cornerdata(i,1)+7, og_cornerdata(i,2), num2str(i), 'Color', 'red', 'FontSize', 12);
end
% Translate flake coordinates using the marker translation:
cornerdata = og_cornerdata' - Translation;  % subtract translation from x-coords
cornerdata(2,:) = cornerdata(2,:) * (-1);
cornerpositions{1} = cornerdata;

%% Step 3: Define Current Electrode (Channel) Axis from Flake Boundary
title('Click on two points along the flake boundary for the CURRENT electrode axis');
[current_x, current_y] = ginput(2);
current_axis = [current_x, current_y];  % Two points along the flake edge
hold on;
plot(current_axis(:,1), current_axis(:,2), 'm+', 'MarkerSize', 10, 'LineWidth', 2);
disp('Current electrode axis defined.');

%% Step 4: Generate Current Electrodes (Source/Drain)
% Prompt for current electrode dimensions (in pixels)
prompt = {'Enter current electrode length (pixels):', 'Enter current electrode width (pixels):'};
dlgtitle = 'Current Electrode Dimensions';
dims = [1 50];
definput = {'80','20'};
currDims = inputdlg(prompt, dlgtitle, dims, definput);
currLength = str2double(currDims{1});
currWidth  = str2double(currDims{2});

% Compute channel vector and unit vector
chVec = current_axis(2,:) - current_axis(1,:);
L_channel = norm(chVec);
u = chVec / L_channel;  % Unit vector along channel
p = [-u(2), u(1)];      % Perpendicular unit vector

% Create a rectangle (local coordinates) for current electrodes
rect_curr = [ currLength/2,  currWidth/2;
              currLength/2, -currWidth/2;
             -currLength/2, -currWidth/2;
             -currLength/2,  currWidth/2 ];
         
% Rotation matrix for current electrodes (aligned with channel)
theta = atan2(u(2), u(1));
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];

% Current electrode 1: centered at the first endpoint of channel
center1 = current_axis(1,:);
electrode1 = (R * rect_curr')' + center1;
% Current electrode 2: centered at the second endpoint of channel
center2 = current_axis(2,:);
electrode2 = (R * rect_curr')' + center2;

%% Step 5: Generate Voltage Electrodes (Probes)
% Voltage electrodes will be placed perpendicular to the channel.
% They are placed at positions 1/4 and 3/4 along the channel.
prompt = {'Enter voltage electrode length (pixels):', 'Enter voltage electrode width (pixels):', 'Enter lateral offset from channel (pixels):'};
dlgtitle = 'Voltage Electrode Dimensions & Offset';
definput = {'40','15','30'};
voltDims = inputdlg(prompt, dlgtitle, dims, definput);
voltLength = str2double(voltDims{1});
voltWidth  = str2double(voltDims{2});
offsetDist = str2double(voltDims{3});

% Define positions along the channel:
pos1 = current_axis(1,:) + u * (L_channel/4);
pos2 = current_axis(1,:) + u * (3*L_channel/4);

% Local rectangle for voltage electrodes (centered at 0)
rect_volt = [ voltLength/2,  voltWidth/2;
              voltLength/2, -voltWidth/2;
             -voltLength/2, -voltWidth/2;
             -voltLength/2,  voltWidth/2 ];
         
% Voltage electrodes are oriented perpendicular to the channel.
theta_volt = atan2(p(2), p(1));
R_volt = [cos(theta_volt) -sin(theta_volt); sin(theta_volt) cos(theta_volt)];

% For position pos1, create two electrodes (one on each side)
volt1_center = pos1 + p * offsetDist;
volt2_center = pos1 - p * offsetDist;
voltage_electrode1 = (R_volt * rect_volt')' + volt1_center;
voltage_electrode2 = (R_volt * rect_volt')' + volt2_center;

% For position pos2, similarly:
volt3_center = pos2 + p * offsetDist;
volt4_center = pos2 - p * offsetDist;
voltage_electrode3 = (R_volt * rect_volt')' + volt3_center;
voltage_electrode4 = (R_volt * rect_volt')' + volt4_center;

%% Step 6: Translate All Coordinates Using Markers
% The flake and electrode coordinates are now translated relative to Marker 3.
% (If desired, additional scaling or conversion to microns can be applied.)
%
% Here, we translate current and voltage electrode coordinates:
electrodes.current = {electrode1, electrode2};
electrodes.voltage = {voltage_electrode1, voltage_electrode2, voltage_electrode3, voltage_electrode4};

% (Optional) If your chip coordinate conversion uses a change-of-coordinates
% matrix (e.g., from pixel to microns using marker references), apply it here.
% For example, if "changecoordmat" is defined from your marker analysis:
% electrode_chip = changecoordmat * electrode_pixel + Center;

%% Step 7: Plot the Hall Bar Configuration on the Chip Image
figure; imshow(ChipImage); hold on;
plot(electrode1(:,1), electrode1(:,2), 'b-', 'LineWidth', 2);
plot(electrode2(:,1), electrode2(:,2), 'b-', 'LineWidth', 2);
plot(voltage_electrode1(:,1), voltage_electrode1(:,2), 'r-', 'LineWidth', 2);
plot(voltage_electrode2(:,1), voltage_electrode2(:,2), 'r-', 'LineWidth', 2);
plot(voltage_electrode3(:,1), voltage_electrode3(:,2), 'r-', 'LineWidth', 2);
plot(voltage_electrode4(:,1), voltage_electrode4(:,2), 'r-', 'LineWidth', 2);
title('Hall Bar Device Configuration');
legend('Current Electrode 1','Current Electrode 2','Voltage Electrodes');

%% Step 8: Save the Coordinates to a Text File
[savefilename, pathname] = uiputfile('*.txt', 'Save Hall Bar Coordinates As', baseFilename);
if ~isequal(savefilename,0)
    fileid = fopen(fullfile(pathname, savefilename), 'w');
    fprintf(fileid, 'Hall Bar Electrode Coordinates (Pixel Units):\n\n');
    
    fprintf(fileid, 'Current Electrodes (Source/Drain):\n');
    for i = 1:2
        electrode = electrodes.current{i};
        fprintf(fileid, 'Electrode %d:\n', i);
        fprintf(fileid, '%f, %f\n', electrode');
        fprintf(fileid, '\n');
    end
    
    fprintf(fileid, 'Voltage Electrodes (Probes):\n');
    for i = 1:4
        electrode = electrodes.voltage{i};
        fprintf(fileid, 'Electrode %d:\n', i+2);
        fprintf(fileid, '%f, %f\n', electrode');
        fprintf(fileid, '\n');
    end
    fclose(fileid);
    disp(['Hall Bar coordinates saved to: ' fullfile(pathname, savefilename)]);
else
    disp('File save canceled.');
end
