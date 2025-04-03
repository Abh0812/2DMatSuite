%findFlakeCoord.m loads images of the marker chip and allows the user to
%identify the coordinates of nanostructures relative to the alignment marks
%on the chip.

close all

pads_centers = [-2750, 500; -3500, 2200; -2500, 3500; -500, 3500; 500, 3500;
    2500, 3500; 3500, 2200; 2750, 500; 2750, -500; 3500, -2200; 2500, -3500;
    500, -3500; -500, -3500; -2500, -3500; -3500, -2200; -2750, -500];
pad_size = 500;

lettervals = dictionary(["A" "B" "C" "D" "E" "F" "G" "H" "I" "J"], ...
    [0 1 2 3 4 5 6 7 8 9]);

cornerpositions = {};
nn = 0;
electrodes = {};
%Open image file

    %Open image file
    % filename = uigetfile('*.jpg*');
    [filename, folder] = uigetfile('*.*')
    [~, baseFilename, ~] = fileparts(filename)
    % add to path
    addpath(folder);
    % Read image file
    ChipImage = imread(filename);

    % Display image
    imshow(ChipImage);
    hold on;

    % Enable zoom
    zoom on; % Allow user to zoom
    
    % Wait until the user presses Enter after zooming
    title(["ZOOM IN a region w/ 4 crosses, press Enter to continue."]);
    pause; % Wait for user to press Enter

    % Get positions of alignment marks from user
    title(["Click on the (1) TOP RIGHT cross center," "then the (2) BOTTOM RIGHT cross center," "then the (3) BOTTOM LEFT cross center."]);
    disp(["Click on the (1) TOP RIGHT cross center,"; "then the (2) BOTTOM RIGHT cross center,"; "then the (3) BOTTOM LEFT cross center."]);

    % [xmark,ymark] = ginput(4);

    % Number of points to click
    numPoints = 3;
    
    % Preallocate arrays for coordinates
    xmark = zeros(numPoints, 1);
    ymark = zeros(numPoints, 1);
    
    % Loop to get points and plot them immediately
    for i = 1:numPoints
        [xmark(i), ymark(i)] = ginput(1); % Get one point from ginput
        
        % Plot the marker at the clicked point
        plot(xmark(i), ymark(i), 'r+', 'MarkerSize', 10, 'LineWidth', 2); % Red circle markers
        
        % Display the coordinates of the current point
        disp(['Point ' num2str(i) ': (' num2str(xmark(i)) ', ' num2str(ymark(i)) ')']);
    end

    %Identify field
    markNum = baseFilename;
    userQ = ['Letters are currently "' markNum '", type in the letters of the last cross you clicked:'];
    changemarkNum = inputdlg(userQ, 'Input Dialog', [1 50], {markNum});

    if ischar(changemarkNum{1}) && length(changemarkNum{1}) == 4
        markNum = changemarkNum{1};
    else
    end
    %Determine translation of the alignment mark to 0,0
    Translation = [xmark(3);ymark(3)];
    
    xmark(1)=xmark(3);
    ymark(2)=ymark(3);

    %Translate alignment mark coordinates such that origin corresponds to Mark
    %4 position
    xmark = xmark-Translation(1);
    ymark = ymark-Translation(2);
    ymark = ymark*(-1);

    %Get positions of nanostructures from user
    
    % Reset the axes to the original limits
    xlim([1 size(ChipImage, 2)]); % Reset x-axis limits based on image size
    ylim([1 size(ChipImage, 1)]); % Reset y-axis limits based on image size
    
       title('Click on the corners of the nanostructure with the mouse');
       % Prompt the user to draw an outline around the flake
        hPolygon = drawpolygon('LineWidth', 1, 'Color', 'g');  % Let the user draw the boundary
        vertices = hPolygon.Position; % This returns an Nx2 matrix where N is the number of vertices

       cornerdata = vertices;
       og_cornerdata = cornerdata;

       plot(og_cornerdata(:,1), og_cornerdata(:,2),'r+','MarkerSize', 20)
       for i = 1:length(og_cornerdata(:,1))
        % Adjust the text position slightly to avoid overlap with the dot
        text(og_cornerdata(i,1) + 7, og_cornerdata(i,2), num2str(i), 'Color', 'red', 'FontSize', 12);
       end

       %Translate positions such that origin corresponds to Mark4 position
       cornerdata = cornerdata' - Translation;
       cornerdata(2,:) = cornerdata(2,:)*(-1);
       cornerpositions{1} = cornerdata;
       % jj = jj + 1;
    


    %electrode axis
    axis_num = 0;
    % axis = inputdlg("Enter the 2 indices of the coordinates for the axis. " + ...
    %        "For example, if it was the first point clicked, enter 1. Do not separate inputs. Enter 0 if none: ");

        % Prompt user to click on polygon vertices
    title('Click on the 2 polygon vertices for the 1st electrode.');
    
    % Number of vertices to select
    numSelections = 2;
    selectedIndices = zeros(numSelections, 1);

    for i = 1:numSelections
        [xClick, yClick] = ginput(1); % Get one point from ginput
        
        % Find the closest vertex to the clicked point
        distances = sqrt((vertices(:, 1) - xClick).^2 + (vertices(:, 2) - yClick).^2);
        [~, closestIndex] = min(distances); % Get index of the closest vertex
        
        % Store the index
        selectedIndices(i) = closestIndex;
        
        % Mark the selected vertex
        plot(vertices(closestIndex, 1), vertices(closestIndex, 2), 'g+', 'MarkerSize', 10, 'LineWidth', 2); % Green cross marker
        
        % Display the selected index
        disp(['Selected vertex index: ' num2str(closestIndex)]);
    end

    % while axis ~= 0
        axis_num = axis_num +1;
        %corner data for electrodes
        y1 = og_cornerdata(selectedIndices(1),2);
        x1 = og_cornerdata(selectedIndices(1),1);

        y2 = og_cornerdata(selectedIndices(2),2);
        x2 = og_cornerdata(selectedIndices(2),1);

        len1 = ((x1-x2)^2+(y1-y2)^2)^0.5;
        len2 = len1+40;
        wid = 10;

        % Compute the vector from the points
        vector = [x1 y1]- [x2 y2];
        % Compute the angle of the vector
        angle = atan2(vector(2), vector(1)) + pi/2;
        % Create the rotation matrix using the computed angle
        R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
        % Rotate the 2x2 matrix A by the angle of the vector x
        
        electrode_corners = [wid/2, len2/2; wid/2,-len2/2; -wid/2,-len2/2; -wid/2,len2/2];
        % Apply the rotation matrix to each point
        electrode_corners = (R * electrode_corners')';  % Transpose, multiply, then transpose back
        electrode_corners = electrode_corners + [0.5*(x1+x2), 0.5*(y1+y2)];
        
        plot(electrode_corners(:,1),  electrode_corners(:,2),'b-','MarkerSize', 20)
        og_electrode_corners= electrode_corners;
        %Translate positions such that origin corresponds to Mark4 position
        electrode_corners = electrode_corners' - Translation;
        electrode_corners(2,:) = electrode_corners(2,:)*(-1);
        electrodes{axis_num} = electrode_corners;
        % axis = input("Enter 0: ");
    % end

    %parallel electrode
    parallel_electrodes = {};
    for ll = 1:(length(electrodes))
        title("Click on where the center of the parallel electrode should be:")
        [center2_x, center2_y] = ginput(1);
        parallel_electrode = og_electrode_corners - [0.5*(x1+x2), 0.5*(y1+y2)] + [center2_x, center2_y];

        plot(parallel_electrode(:,1),  parallel_electrode(:,2),'b-','MarkerSize', 20)
        %Translate positions such that origin corresponds to Mark4 position
        parallel_electrode = parallel_electrode' - Translation;
        parallel_electrode(2,:) = parallel_electrode(2,:)*(-1);
        parallel_electrodes{ll} = parallel_electrode;
        %parallel_electrodes{ll} = mat * [[wid/2;len2/2] [wid/2;-len2/2] [-wid/2;-len2/2] [-wid/2;len2/2]] + [center2_x; center2_y];
    end

    %Find change of coordinates matrix to translate corner positions from image
    %coordinates to positions (in microns) relative to Mark4
    basisvec1 = [xmark(2);ymark(2)];
    basisvec2 = [xmark(1);ymark(1)];
    basismat = [basisvec1 basisvec2];
    changecoordmat = [200,0;0,200]*basismat^(-1);

    %Determine absolute positions of the marker in microns
    Center = [-2900+2000*lettervals(markNum(3))+200*lettervals(markNum(4));
        2900-2000*lettervals(markNum(1))-200*lettervals(markNum(2))];

    %Translate coordinates of corners to AutoCAD coordinates
    for kk = 1:(length(cornerpositions)-nn)
        tempmat = cornerpositions{kk+nn};
        tempmat = changecoordmat*tempmat;
        tempmat = tempmat + Center;
        cornerpositions{kk+nn} = tempmat;
    end

    for kk = 1:(length(electrodes)-nn)
        tempmat = electrodes{kk+nn};
        tempmat = changecoordmat*tempmat;
        tempmat = tempmat + Center;
        electrodes{kk+nn} = tempmat;

        parallel_electrodes{kk+nn}= changecoordmat*parallel_electrodes{kk+nn}+Center;
    end
    nn = length(cornerpositions);


%% Choose the closest pad for routing
special = cornerpositions{1}(:,1)';

[padIndex, padCoord] = get_closest_pad(special, pads_centers);

% disp_mes = "Closest pad is at ["+ string(padCoord(1))+"," +string(padCoord(2))+ "], pad #" + padIndex+", press 0 to change, otherwise press enter.";
% padChange = inputdlg(disp_mes);
% 
% if padChange == {0}
%     padIndex = input('Type desired pad index. <a href="https://github.com/morganblevins/scanning-photocurrent-microscope/tree/main/nanofab-electrode-design">Reference PPMS electrode map on GitHub</a>.');
%     padCoord = pads_centers(padIndex,:)
% end

%%
%Write corner position data to file

    % disp('Create a file for saving points.  Use extension .txt');
    % savefilename = uiputfile;
    fileName = markNum;
    [savefilename, pathname] = uiputfile('*.txt', 'Save Coordinates As', markNum);

    fileid = fopen(savefilename,'w');

    fprintf(fileid,'preferred pad\n');

    tempmat = padCoord;
    fprintf(fileid,'%f,%f\n',tempmat);
    fprintf(fileid,'\n');

    fprintf(fileid,'nanostructures\n');

    for ll = 1:length(cornerpositions)
       tempmat = cornerpositions{ll};
       fprintf(fileid,'%f,%f\n',tempmat);
       fprintf(fileid,'\n');
    end

    fprintf(fileid,'electrodes\n');

    for ll = 1:length(electrodes)
       tempmat = electrodes{ll};
       fprintf(fileid,'%f,%f\n',tempmat);
       fprintf(fileid,'\n');
    end

    fprintf(fileid, 'parallel electrodes\n');
    for ll = 1:length(electrodes)
       tempmat = parallel_electrodes{ll};
       fprintf(fileid,'%f,%f\n',tempmat);
       fprintf(fileid,'\n');
    end

    % fprintf(fileid, 'cross\n');
    % for ll = 1:length(xmark)
    %    tempmat = [xmark(ll); ymark(ll)]+Center;
    %    fprintf(fileid,'%f,%f\n',tempmat);
    %    % fprintf(fileid,'\n');
    % end

    fclose(fileid);


%%

% FUNCTIONS

function [padIndex, padCoord] = get_closest_pad(special, pads_centers)
    % Function to return the index and coordinates of the pad closest to a special coordinate

    % Initialize minimum distance to a large value
    minDistance = inf;
    padIndex = -1;
    padCoord = [];

    % Iterate through each pad center
    for i = 1:length(pads_centers(:,1))
        % Calculate the distance between the current pad center and the special coordinate
        currentDistance = distance(pads_centers(i,:), special);
        
        % Update the minimum distance and corresponding pad index and coordinates if a closer pad is found
        if currentDistance < minDistance
            minDistance = currentDistance;
            padIndex = i;
            padCoord = pads_centers(i,:);
        end
    end
end

function d = distance(coord1, coord2)
    % Function to calculate the distance between two coordinates
    d = sqrt(sum((coord1 - coord2).^2));
end



