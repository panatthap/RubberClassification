clc
close all
clear all
%%
format short
load database.mat
load t.mat y
for i = 1:127
    ii = num2str(i);
    strr = strcat('test_data/', ii,'.jpg');
    Img = imread(strr);
    Img = imresize(Img, [512 512]);
    I = rgb2gray(Img);
    %% convert to HSV color
    cc = rgb2hsv(Img);
    %% spilt HSV color
    hh = cc(:,:,1);
    ss = cc(:,:,2);
    vv = cc(:,:,3);
%     figure, imshow(cc),impixelinfo
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
%     figure, imshow(out),impixelinfo
    converted = uint8(out) .* I;
%     figure, imshow(converted),impixelinfo
    %% Finding GLCM
    % glcms = graycomatrix(I,Name,Value) returns one or more gray-level co-occurrence matrices, ...
    % depending on the values of the optional name-value pair arguments.
    
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
%     figure, imshow(edged),impixelinfo
    data_edged(i,:) = std2(edged);
    data_edgedss = std2(edged);
    %% Finding gabor
    gaborArray = gabor([4 8],[0 90]); % g = gabor(wavelength,orientation) creates a Gabor filter with the specified wavelength (in pixels/cycle) and orientation (in degrees). 
    gaborMag = imgaborfilt(I,gaborArray);
    data_gaborMag(i,:) = std2(gaborMag);
    data_gaborMagss = std2(gaborMag);
    %%%
    [labeledImage numberOfBlobs] = bwlabel(out, 8); % L = bwlabel(BW,conn) returns a label matrix, where conn specifies the connectivity.
    blobMeasurements = regionprops(labeledImage, 'Centroid', 'Orientation');
    xCenter = blobMeasurements(1).Centroid(1);
    yCenter = blobMeasurements(1).Centroid(2);
    %subplot(1,2,2);
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
    %angless = std2(angles);
    %distancess = std2(distances);
    %%
    stats = regionprops('table',out,'Orientation','Perimeter',...
    'MajorAxisLength','MinorAxisLength');
    Orientationss = mean(stats.Orientation)
    Perimeterss = mean(stats.Perimeter);
%     Orientations(i,:) = mean(Orientationss);
%     Perimeters(i,:) = mean(Perimeterss);
    %%
    format short
    test_data = [data_gaborMagss data_edgedss dataGLCMss Orientationss Perimeterss];
    [rr cc] = size(p);
    %t=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1];
    t = y';
    [pn,ps] = mapminmax(p');
    [tn,ts] = mapminmax(t);
    net=newff(pn,tn,[52 30 10 30 10 2],{'tansig' 'tansig' 'tansig' 'tansig' 'tansig' 'purelin'},'trainrp');

    net.trainParam.show = 1000;
    net.trainParam.epochs = 5000;
    net.trainParam.goal = 1e-5;

    net = train(net, pn, tn);
    aaa = sim(net, test_data');
    y2(i,:) = round(mapminmax('reverse',aaa,ts));
    
end
error = abs(y2-t');
[errorr errorc] = size(error);
total = errorr;
yy = find(error==0);
[rrr ccc] = size(yy);
format short
percenge_acuracy = ((total-rrr)/total)*100;
formatSpec = 'percenge_acuracy is %f percenge \n';
fprintf(formatSpec,percenge_acuracy)