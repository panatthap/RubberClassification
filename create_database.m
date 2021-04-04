clc
close all
clear all
%%
format short
for i = 1:10
    ii = num2str(i);
    strr = strcat('test_data/', ii,'.jpg')
    Img = imread(strr);
    Img = imresize(Img, [512 512]);
    I = rgb2gray(Img);
    cc = rgb2hsv(Img);
%     imshow(cc), impixelinfo
    %% spilt HSV color
    hh = cc(:,:,1);
    ss = cc(:,:,2);
    vv = cc(:,:,3);
    [r1 c1] = size(hh);
    %% segmentation with HSV color space
    for i1 = 1:r1
        for j1 = 1:c1
            if hh(i1,j1) >= 0.10 && hh(i1,j1) <= 0.30 &&...
                    ss(i1,j1) >= 0.50 && ss(i1,j1) <= 0.80 &&...
                    vv(i1,j1) >= 0.20 && vv(i1,j1) <= 0.80
                out(i1,j1) = 255;
            else
                out(i1,j1) = 0;
            end
        end
    end
     %% Post-processing
    out = im2bw(out); % Binary image
    out = bwareaopen(out, 100); % Remove object < 100 out from image
    out = imfill(out, 'holes'); % fill holes
    converted = uint8(out) .* I;
%     imshow(out),impixelinfo
    GLCM2 = graycomatrix(converted,'Offset',[-1 0;0 1]);
    stats = GLCM_Features1(GLCM2,0);
    dataGLCM(i,:) = [stats.autoc stats.contr stats.corrm stats.corrp stats.cprom...
        stats.cshad stats.dissi stats.energ stats.entro stats.homom...
        stats.homop stats.maxpr stats.sosvh stats.savgh stats.svarh...
        stats.senth stats.dvarh stats.denth stats.inf1h stats.inf2h...
        stats.indnc stats.idmnc stats.contr stats.contr stats.contr];
    dataGLCMss = [stats.autoc stats.contr stats.corrm stats.corrp stats.cprom...
        stats.cshad stats.dissi stats.energ stats.entro stats.homom...
        stats.homop stats.maxpr stats.sosvh stats.savgh stats.svarh...
        stats.senth stats.dvarh stats.denth stats.inf1h stats.inf2h...
        stats.indnc stats.idmnc stats.contr stats.contr stats.contr];
   %% Finding Edge
    edged = edge(converted, 'canny');
%     imshow(edged),impixelinfo
    data_edged(i,:) = std2(edged);
    data_edgedss = std2(edged);
    %% Finding gabor
    gaborArray = gabor([4 8],[0 90]); % g = gabor(wavelength,orientation) creates a Gabor filter with the specified wavelength (in pixels/cycle) and orientation (in degrees). 
    gaborMag = imgaborfilt(I,gaborArray);
    data_gaborMag(i,:) = std2(gaborMag);
    data_gaborMagss = std2(gaborMag);
    %%%
    [labeledImage numberOfBlobs] = bwlabel(out, 8);
    blobMeasurements = regionprops(labeledImage, 'Centroid', 'Orientation');
    xCenter = blobMeasurements(1).Centroid(1);
    yCenter = blobMeasurements(1).Centroid(2);
    imshow(I);
    axis on;
    hold on;
    % Plot the centroid.
    plot(xCenter, yCenter, 'b+', 'MarkerSize', 10, 'LineWidth', 3);
    hold on;
    boundaries = bwboundaries(out);
    for k = 1 : length(boundaries)
      thisBoundary = boundaries{k};
      plot(thisBoundary(:,2), thisBoundary(:,1), 'b', 'LineWidth', 2);
      numberOfBoundaryPoints = length(thisBoundary);
      angle = 0: 10 : 360;
      for a = 1 : length(angle)
          xb = thisBoundary(a,2);
          yb = thisBoundary(a,1);
          angles(a,:) = atand((yb-yCenter) / (xb-xCenter));
          distances(a,:) = sqrt((xb-xCenter).^2+(yb-yCenter).^2);
      end
    end
    stats = regionprops('table',out,'Orientation','Perimeter',...
    'MajorAxisLength','MinorAxisLength');
    Orientationss = mean(stats.Orientation);
    Perimeterss = mean(stats.Perimeter);
    Orientation(i,:) = mean(stats.Orientation);
    Perimeter(i,:) = mean(stats.Perimeter);

end
format short
database = [data_gaborMag data_edged dataGLCM Orientation Perimeter];
save('database.mat','database')