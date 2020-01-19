%%
% clear all;

% basePath = "~/Dropbox/xps_matlab/jsonData/";
basePath = "~/Dropbox/xps/eQ/build/";

% LOAD FILE
% [file path] = uigetfile('~/Dropbox/xps/eQ/build/*.json', 'Open File');
[file path] = uigetfile(sprintf("%s*.json", basePath), 'Open File');
thisExtension = '.json';


titrationCurve = [];
meanReadout= [];
stdReadout = [];



%==========================================================================
%  SEQUENCE ALL .JSON FILES IN FOLDER:
%==========================================================================
% for whichFile = 1 : fileCounter
    
    % LOAD FILE, READ DIRECTLY JSON FILE VIA MATLAB INTERFACE
%     dataIn = fileread(myFiles(whichFile));
    dataIn = fileread(fullfile(path,file));
    fprintf('decoding JSON file...%s\n', file);
    jsonData =  jsondecode(dataIn);
    
    
    % EXTRACT META DATA, FRAMED DATA, AND CELL DIVISION DATA STRUCTURES:
    lengthCellDataField = length(jsonData.frames(1).cells(1).d);
    maxValues = zeros(1, lengthCellDataField);
    minValues = zeros(1, lengthCellDataField);

    % CELL DATA ARRAY KEY:
    CENTER_X = 1;
    CENTER_Y = 2;
    ANGLE = 3;
    LENGTH = 4;
    COMP = 5;
    HSL = 6;

    simNumber   = jsonData.simNumber
    simdt       = jsonData.parameters.dt;
    timeStamp = num2str(jsonData.timeSinceEpoch);
    trapHeight  = jsonData.parameters.simulationTrapHeightMicrons;
    trapWidth   = jsonData.parameters.simulationTrapWidthMicrons;

    % DETERMINE NUMBER OF CELLS TO TRACK:
    numFrames = length(jsonData.frames);
    numCells = zeros(numFrames,1);
    maxCellNumber = 0;
    for i = 1:numFrames
        numCells(i) = length(jsonData.frames(i).cells);
        for j=1:numCells(i)
            thisCellNumber = jsonData.frames(i).cells(j).i;
            if(thisCellNumber > maxCellNumber)
                maxCellNumber = thisCellNumber;
            end
        end
    end
    % DEFINE THE CELLS DATA STRUCTURE: (indexed by cellID, arrays of frame
    % data)
    cells(1).divisions = [];
    cells(1).length = [];
    cells(1).angle = [];
    cells(1).xy= [];
    cells(1).parent = [];
    cells(1).dlist = [];
    cells(maxCellNumber) = cells(1);

    % DEFINE THE FRAMES DATA STRUCTURE: (data across all frames)
    frameData.strainRatio = zeros(numFrames,1);
    frameData.avg_comp = zeros(numFrames,1);
    frameData.avg_length = zeros(numFrames,1);
    frameData.qScalar = zeros(numFrames,1);
    frameData.forceX = zeros(numFrames,1);
    frameData.forceY = zeros(numFrames,1);

    % FOR X,Y GRID AND Y-AVERAGED ARRAYS PER FRAME
        global qGrid;
        qGrid =  zeros(trapHeight+1, trapWidth+1, numFrames);
        xFvec =  zeros(trapWidth+1, numFrames);
        xQvec =  zeros(trapWidth+1, numFrames);


    %pretty-print progress counters:    
    displayCounter=1;
    displaySlice = numFrames/10.0; %display progress in 10% intervals
    display('Beginning Frame iteration...');

%==========================================================================
    % LOOP PER FRAME:  jsonData.frames(i)
%==========================================================================
    for i= 1:numFrames

        %prettry-print progress
%         if(i > displayCounter*displaySlice)
%             thisMessage = sprintf('Completed %d percent of frames...', 10*displayCounter);
%             display(sprintf('Completed %d percent of frames...', 10*displayCounter));
%             displayCounter = displayCounter + 1;
%         end
% 
        % HANDLE DIVISION EVENTS        % HANDLE DIVISION EVENTS

        lengthDivArray = size(jsonData.frames(i).divisions);
        
        for d = 1:lengthDivArray(1)
           thisDivision = jsonData.frames(i).divisions(d,:);
               timeStamp    = thisDivision(1);
               parent       = thisDivision(2);
               plength      = thisDivision(3);
               daughter     = thisDivision(4);
               dlength      = thisDivision(5);
           %parent list updated if not a seed cell:    
           if(parent > 0)
                thisDivision = [timeStamp; plength];
                cells(parent).divisions = [cells(parent).divisions, thisDivision];
                cells(parent).dlist = [cells(parent).dlist, daughter];
           end
           %daughter list created with birth entry:
            thisDivision = [timeStamp; dlength];
            cells(daughter).divisions = thisDivision;%birth entry
            cells(daughter).dlist = [];
        end

        % HANDLE CELL DATA:
        thisFrameCells = size(jsonData.frames(i).cells);
        numCells(i) = length(jsonData.frames(i).cells);
            %zero accumulators for averaging across all cells in frame:
            cos2phi=0;  sin2phi=0; fx=0; fy=0;
            xF=zeros(1,trapWidth+1); xFi=zeros(1,trapWidth+1);
            xQ=zeros(1,trapWidth+1); xQi=zeros(1,trapWidth+1);

%==========================================================================
                % LOOP PER CELL::           jsonData.frames(i).cells(j)
%==========================================================================
        for j=1:numCells(i)
            thisCell = jsonData.frames(i).cells(j).i;
            thisStrain = jsonData.frames(i).cells(j).p;
                x = jsonData.frames(i).cells(j).d(1); xi = round(x) + 1;
                y = jsonData.frames(i).cells(j).d(2); yj = round(y) + 1;
                phi = mod(jsonData.frames(i).cells(j).d(3), 2*pi);
                l = jsonData.frames(i).cells(j).d(4);
                comp = jsonData.frames(i).cells(j).d(5);

            %update accumulators:
            frameData.strainRatio(i)    = frameData.strainRatio(i) + thisStrain;    
            frameData.avg_comp(i)       = frameData.avg_comp(i) + comp;
            frameData.avg_length(i)     = frameData.avg_length(i) + l;
            cos2phi                     = cos2phi + cos(2.0*phi);
            sin2phi                     = sin2phi + sin(2.0*phi);
            fx                          = fx + power(comp*cos(phi), 2.0);  
            fy                          = fy + power(comp*sin(phi), 2.0);

            %update cell data arrays (per frame)
            cells(thisCell).length      = [cells(thisCell).length, l];
            cells(thisCell).angle       = [cells(thisCell).angle, phi];
            cells(thisCell).xy          = [cells(thisCell).xy, [x; y] ];

            qGrid(yj,xi,i) = cos(2*phi) + 1;
            %laterl force averaged over column
                    xF(xi) = xF(xi) + power(comp*cos(phi), 2.0);  
                    %increment index counter:
                    xFi(xi) = xFi(xi) + 1;

                    %iterate the data array and check new max/min
                    for dvalue = 1:lengthCellDataField
                        thisvalue = jsonData.frames(i).cells(j).d(dvalue);
                        if thisvalue > maxValues(dvalue)
                            maxValues(dvalue) = thisvalue;
                        end
                        if thisvalue < minValues(dvalue)
                            minValues(dvalue) = thisvalue;
                        end
                    end
        end%end cells in this frame

                % COMPUTE AVERAGES PER FRAME:
                frameData.strainRatio(i)    = frameData.strainRatio(i)/numCells(i);
                frameData.avg_comp(i)       = frameData.avg_comp(i)/numCells(i);
                frameData.avg_length(i)     = frameData.avg_length(i)/numCells(i);
                cos2phi                     = cos2phi/numCells(i); 
                sin2phi                     = sin2phi/numCells(i); 
                frameData.qScalar(i)        = sqrt(cos2phi*cos2phi + sin2phi*sin2phi);
                frameData.forceX(i)         = fx/numCells(i);
                frameData.forceY(i)         = fy/numCells(i);
                %average columns (but only if counter > 0)
                    ix = find(xFi>0); 
                    xF(ix) = xF(ix) ./ xFi(ix);
                    xFvec(:,i) = xF;

    end%end frames
%==========================================================================
    % AVERAGE DATA OVER LAST n FRAMES:
%==========================================================================
    
%     colFrameWidth=10;%#frames-1; 0 ==> 1 frame
    colFrameWidth=14;%#frames-1; 0 ==> 1 frame
    readoutX = zeros(1,trapWidth);
    rXcount = zeros(1,trapWidth);
    
    for i=numFrames-colFrameWidth:numFrames
        readoutData = zeros(1,numCells(i));
        meanKH=0; counter=0;  meanH = 0;
        % LOOP OVER ALL CELLS IN THIS FRAME:
        for j=1:numCells(i)
            thisCell = jsonData.frames(i).cells(j).i;
            thisStrain = jsonData.frames(i).cells(j).p;

            %LIMIT TO MIDDLE OF TRAP HERE:
%             ymid = trapHeight/2;
%             y = jsonData.frames(i).cells(j).d(2);
%             if(abs(y - ymid) > 2)
%                 continue;
%             end
            %END LIMITING REGION
            
            counter = counter + 1;

            whichData = 6;%HSL
            rThresh = 0.5;
%             whichData = 6;   %HSL
%             rThresh = 1e4;


            readout = jsonData.frames(i).cells(j).d(whichData);
            readoutData(counter) = readout;

            meanH = meanH + jsonData.frames(i).cells(j).d(6);
%             meanKH = meanKH + jsonData.frames(i).cells(j).d(8);
            
            
            x = round(jsonData.frames(i).cells(j).d(1));
            if(x==0) x=1; end
            
            readoutX(x) = readoutX(x) + readout;
            rXcount(x) = rXcount(x) + 1;
            
        end
    end %last n frames

    readoutData = readoutData(1:counter);
    
    %average columns (but only if counter > 0)
        ix = find(rXcount>0); 
        readoutX(ix) = readoutX(ix) ./ rXcount(ix);
    
    figure;
    plot(1:trapWidth, readoutX);
    
%     iptgValues = [];
%     flatResponses = [];
    
%     iptg = jsonData.parameters.MODULUS_IPTG;
%     iptgValues = [iptgValues iptg];
%     flatResponses = [flatResponses readoutX'];

    return;
    
    
    % DISPLAY HISTOGRAM OF READOUT DATA:
    figure(200);
    numBins = 100;
    h = histogram(readoutData, numBins);    
    binsOverPointFive = find(h.BinEdges(1:(end-1)) >= rThresh, h.NumBins);
    activated = sum(h.BinCounts(binsOverPointFive))/sum(h.BinCounts)
    
    
%     display(sprintf("iptg = %i", jsonData.parameters.MODULUS_IPTG));
%     display(sprintf("meanKH = %i", meanKH/counter));
%     display(sprintf("meanH = %i", meanH/counter));

%     
%     readoutValues{whichFile} =  readoutData;
%     readoutDistributions{whichFile} = {h.BinEdges(1:(end-1)), h.BinCounts};
%     iptgValues(whichFile) = jsonData.parameters.MODULUS_IPTG;
% 
%     
%     axis([0 1 -inf inf]);
%     ax = gca;
%     FS=20;
%     ax.FontSize = FS; 
%     xlabel('READOUT'); 
%     ylabel('Histogram');
%     iptgString = sprintf("IPTG SIGNAL = %d \\muM", jsonData.parameters.MODULUS_IPTG);
%     title({'ACTIVATION READOUT', ...
%         iptgString, sprintf("\\phi = %d", activated)}, ...
%         'FontSize', FS);
%     
%     %APPEND THE TITRATION DATA TO SAVE TO FILE:
%     titrationCurve = [titrationCurve; [jsonData.parameters.MODULUS_IPTG, activated]];
%     meanReadout = [meanReadout; [jsonData.parameters.MODULUS_IPTG, mean(readoutData)]];
%     stdReadout = [stdReadout; [jsonData.parameters.MODULUS_IPTG, std(readoutData)]];
%     
% end %for whichFile = 1 : fileCounter

%==========================================================================
    % PLOT TITRATION CURVE FOR ALL FILES' DATA:
%==========================================================================
%%
figure;
semilogx(titrationCurve(:,1), titrationCurve(:,2), 's');
    title({'ACTIVATION READOUT', ...
        sprintf("%s", jsonData.parameters.MODULUS_option)}, ...
        'FontSize', 20);
axis([-Inf Inf 0 1]);
xlabel('IPTG \muM', 'FontSize', 20); 
ylabel('Readout Fraction Activated', 'FontSize', 20);

figure;
semilogx(meanReadout(:,1), meanReadout(:,2), 's');
    title({'ACTIVATION READOUT', ...
        sprintf("%s", jsonData.parameters.MODULUS_option)}, ...
        'FontSize', 20);
axis([-Inf Inf 0 1]);
xlabel('IPTG \muM', 'FontSize', 20); 
ylabel('Mean Readout', 'FontSize', 20);  
hold on;
errorbar(meanReadout(:,1), meanReadout(:,2), stdReadout(:,2), 's');
hold off;



% [dsort isort] = sort(iptgValues);
% figure;
% for rplot=[1, 15, length(iptgValues)]   
%     thisData = readoutDistributions{rplot};
%     bar(thisData{1}, thisData{2});  hold on;
% end
% hold off;


%==========================================================================
%   SAVE DATA TO .MAT FILE:
%==========================================================================
% saveFileName = sprintf("%sjson%s.mat", path, jsonData.parameters.MODULUS_option);
% saveFileName = sprintf("json%s.mat", jsonData.parameters.MODULUS_option);
% save (saveFileName, 'titrationCurve');

oldPath = cd(path);
cd ..
saveFileName = sprintf("readouts%s.mat", jsonData.parameters.MODULUS_option);
save (saveFileName, 'iptgValues',  'readoutValues');

cd(oldPath);
%==========================================================================
    return;
%==========================================================================

%==========================================================================

            
%==========================================================================
% LOAD ALL .mat FILES IN FOLDER AS TITRATION CURVES:
%==========================================================================
%%
thisPath = '~/Dropbox/xps_matlab/jsonData/';
% thisPath = '~/Dropbox/xps_matlab/jsonData/var0_100x800scale4/';
% thisPath = '~/Dropbox/xps_matlab/jsonData/var0_4fb2_50x800s4/';

[file selpath] = uigetfile(sprintf("%s*.mat", thisPath), 'Open File');

% selpath = uigetdir(thisPath)
dinfo = dir(sprintf("%s/*.mat", selpath))
nfiles = length(dinfo);

titrations =[];
legendNames = [];

inputConcs = cell(1,nfiles);
meanReadouts = cell(1,nfiles);
stdReadouts = cell(1,nfiles);

for i = 1 : nfiles  %index number of topologies... 
  filename = fullfile(dinfo(i).folder, dinfo(i).name);
  [filepath,name,ext] = fileparts(filename);
  load (filename);
  
  if strcmp(name, 'readouts-D-F') xi = 1; 
  elseif strcmp(name, 'readouts-D+F') xi=2; 
  elseif strcmp(name, 'readouts+D-F') xi=3; 
  else xi=4; end
  
  xmeans = zeros(1,length(readoutValues));
  xstds = zeros(1,length(readoutValues));
  
  for j=1:length(readoutValues)
      xmeans(j) = mean(readoutValues{j});
      xstds(j) = std(readoutValues{j});
  end
  
  inputConcs{xi} = iptgValues;
  meanReadouts{xi} = xmeans;
  stdReadouts{xi} = xstds;
%   titrationCurve = sortrows(titrationCurve);%sort by iptg induction value
%   titrations = [titrations titrationCurve];
%   legendNames = [legendNames; name];
end


figure;
myorange = 1/255 * [204,85,0]


subplot(1,2,1);

    semilogx(inputConcs{1}, meanReadouts{1}, 'ks');  hold on;
    semilogx(inputConcs{2}, meanReadouts{2}, 'bo'); 
    errorbar(inputConcs{1}, meanReadouts{1}, stdReadouts{1}, 'ks');
    errorbar(inputConcs{2}, meanReadouts{2}, stdReadouts{2}, 'bo');
    hold off;
ax = gca;
ax.FontSize = 16;     
title({'READOUT TITRATIONS', 'NO SIGNALING'}, 'FontSize', 20);
axis([-Inf Inf 0 1]);
xlabel('IPTG \muM', 'FontSize', 20); 
ylabel('READOUT VALUE', 'FontSize', 20);  
legend({'NO FEEDBACK', 'WITH FEEDBACK'}, 'Location', 'southeast', 'FontSize', 10);  


subplot(1,2,2);

    semilogx(inputConcs{3}, meanReadouts{3}, 'bs');  hold on;
    semilogx(inputConcs{4}, meanReadouts{4}, 'o', 'Color', myorange); 
    errorbar(inputConcs{3}, meanReadouts{3}, stdReadouts{3}, 'bs');
    errorbar(inputConcs{4}, meanReadouts{4}, stdReadouts{4}, 'o', 'Color', myorange);
    hold off;

ax = gca;
ax.FontSize = 16; 
title({'READOUT TITRATIONS', 'HSL SIGNALING'}, 'FontSize', 20);
axis([-Inf Inf 0 1]);
xlabel('IPTG \muM', 'FontSize', 20); 
ylabel('READOUT VALUE', 'FontSize', 20);  
legend({'NO FEEDBACK', 'WITH FEEDBACK'}, 'Location', 'southeast', 'FontSize', 10);  

return;


% figure;
% %skip every other column (repeats input conc.)
% plot(titrations(:,1), titrations(:, 2:2:end), 's-', 'MarkerSize',5, 'LineWidth', 3);
% legend(legendNames, 'Location', 'southeast', 'FontSize', 15);
% axis([0 100 -Inf Inf]);
%     title(selpath, ...
%         'FontSize', 10, 'Interpreter', 'none');
% 
% figure;
% semilogx(titrations(:,1), titrations(:, 2:2:end), 's-', 'MarkerSize',5, 'LineWidth', 3);
% legend(legendNames, 'Location', 'southeast', 'FontSize', 15);
% axis([0 100 -Inf Inf]);
%     title(selpath, ...
%         'FontSize', 10, 'Interpreter', 'none');
% 
% return;

%==========================================================================
% FOR PRINTING IPTG INDUCTION VALUES FOR USE IN .CPP PROGRAM
%==========================================================================
%%
xa =power(10, linspace(0.5, 2, 20));
xs = string(xa);
xfs = xs(1);
for i=2:length(xs)
    xfs = sprintf("%s, %s", xfs, xs(i));
end

xfs
%==========================================================================
% CREATE DIVISION TIME VS. BIRTH LENGTH SCATTER PLOT
%==========================================================================
%%

%  now we want to plot tb vs. vb (where tb is "from vb birth time")
% only consider one division length from each division, thus iterate the
% cell list: divisions array
simDivisions = [];
for i=1:length(cells)
    [h w] = size(cells(i).divisions);
    for j=1:w-1
        vb = cells(i).divisions(2,j);
%         if(cells(i).divisions(1,j+1) > 500)
        if(true)
            tb = cells(i).divisions(1,j+1) - cells(i).divisions(1,j);
            simDivisions = [simDivisions [vb; tb; cells(i).divisions(1,j); i]];
        end
    end
end

figure(20); clf;

    plot(simDivisions(1,:), simDivisions(2,:), 's'); hold on;

xlabel('V_b BIRTH LENGTH (um)');  ylabel('t_d TIME TO DIVISION (min)');

numBins = 100;
[N,edges, bin] = histcounts(simDivisions(1,:), numBins);
mTime = zeros(1,length(N));
sTime = zeros(1,length(N));
for ibin = 1:length(N)
    theseIndices = find(bin==ibin);
    mTime(ibin) = mean(simDivisions(2,theseIndices));
    sTime(ibin) = std(simDivisions(2,theseIndices));
end
% errorbar(edges(1:end-1), mTime, sTime);

    plot(edges(2:end), mTime, 'r', 'LineWidth', 3); hold off;

meanbl = mean(simDivisions(1,:));
meandt = mean(simDivisions(2,:));
blstring = sprintf('MEAN BIRTH LENGTH = %.2f \\mum', meanbl);
dtstring = sprintf('MEAN DIVISION TIME = %.2f min', meandt);
axis([meanbl*0.9 meanbl*1.1 -inf inf]);
legend({'(V_b, t_d) scatter', 'mean t_d'}, 'Location', 'northeast');
title({'DIVISION TIME VS. BIRTH LENGTH', ...
    blstring, dtstring}, ...
    'FontSize', 20);


%%

%==========================================================================
%  CREATE QSCALARE, COMPRESSION AND FORCE PLOTS:
%==========================================================================
%%
FS = 30;
figure(11); clf;

    plot(frameData.qScalar, 'LineWidth', 3); hold off;

axis([-Inf Inf 0 1.1]);
ax = gca;
ax.FontSize = FS; 
xlabel('Time (mins)'); 
ylabel('Whole-Trap q-Scalar');
% legend({'Strain Ratio','q-Scalar'}, 'Location', 'northeast');
title({'Q ORDER PARAMETER VS. TIME', ...
    'MEAN DIVISION LENGTH = 2.8 \mum'}, ...
    'FontSize', FS);

figure(10); clf;
% subplot(1,2,1);

yyaxis left;
    plot(frameData.strainRatio, 'LineWidth', 5); hold on;

axis([-Inf Inf 0 1]);
ax = gca;
ax.FontSize = FS; 
xlabel('Time (mins)','FontSize', FS);  ylabel('Strain Ratio','FontSize', FS);

yyaxis right
    plot(frameData.qScalar, 'LineWidth', 3); hold off;

axis([-Inf Inf 0 1]);
ylabel('Whole-Trap q-Scalar');
legend({'Strain Ratio','q-Scalar'}, 'Location', 'northeast');
title({'STRAIN RATIO, Q-SCALAR', ...
%     'MEAN DIVISION LENGTH = 4.2 \mum'}, ...
    'l_d = 4.2--> 2.8 \mum'}, ...
    'FontSize', FS);
subplot(1,2,2);
yyaxis left;
% plot(frameData.avg_comp); hold on;
% plot([abs(frameData.forceX) abs(frameData.forceY)], 'LineWidth', 5); 
plot(abs(frameData.forceX), 'k-', 'LineWidth', 5); hold on;
plot(abs(frameData.forceY), 'b-', 'LineWidth', 5);  hold off;
axis([-Inf Inf 0 Inf]);
xlabel('Time (mins)');  ylabel('Mean-square Force (a.u.)');
yyaxis right
plot(frameData.avg_comp, 'LineWidth', 3); hold off;
axis([-Inf Inf 0 Inf]);
ylabel('Mean Cell-spring Compression');
legend({'[ForceX]^2', '[ForceY]^2','Comp'}, 'Location', 'east');
title({'MEAN SPRING FORCE', ...
    'l_d = 4.2-->2.8 \mum (A@600min)'}, ...
    'FontSize', 20);

%%

% figure(30); clf;
% for i=1:numFrames
%     imagesc(qGrid(:,:,i)); colorbar;
%      M(i) = getframe; 
% end

%  CREATE QSCALAR KYMOGRAPH

figure(40);
meanCos2phi = mean(qGrid,1);
imagesc(squeeze(meanCos2phi)'); colorbar;
xlabel('Trap Position X');  ylabel('Time (mins)');
title('Kymograph: cos(2phi): Time(mins) vs. Trap X position', 'FontSize', 30);

% return;
%==========================================================================


%% 
%==========================================================================
%  LEGI COMPUTATIONS:
%==========================================================================
% 
xfolder = sprintf('%s%d-%d_tiffImages', path, jsonData.timeSinceEpoch, jsonData.simNumber);
mkdir(xfolder);

% threshold = jsonData.parameters.hslThresh;
% dso_hslThresh = jsonData.parameters.dso_hslThresh;
% threshold = .25;


tic
mf = figure(50);  clear v; clear M;
% display("Creating parpool object and generating tif files...");
% poolobj = parpool(5);

% colmeansR = zeros(1,trapWidth+1);
% colmeansA = zeros(1,trapWidth+1);
% colcounter = zeros(1,trapWidth+1);

% colFrameWidth=20;
colFrameWidth=0;

for i=1:numFrames
% for i=numFrames-colFrameWidth:numFrames
% for i=1:100
% for i=790:815
    if(0 == mod(i,10)) fprintf("\ri=%d", i); end
% parfor i=1:numFrames
%     clf;
    readoutData = zeros(1,numCells(i));
    meanKH=0; KHcounter=0;  meanH = 0;
    for j=1:numCells(i)
        thisCell = jsonData.frames(i).cells(j).i;
        thisStrain = jsonData.frames(i).cells(j).p;
        
        whichData = 6;  %legi_A
%         whichData = 7;   %legi_R
        selectedStrain = 0;%1=repressore, 0=activator
%         threshold = maxValues(whichData);
%         threshold = 1600;%hard-coded from -D-F
        
%         whichData = HSL;  zc
%         selectedStrain = 0;%1=repressore, 0=activator

        readout = jsonData.frames(i).cells(j).d(whichData);
        readoutData(j) = readout;
        
%         meanKH = meanKH + jsonData.frames(i).cells(j).d(8);
%         KHcounter = KHcounter + 1;
        
        meanH = meanH + jsonData.frames(i).cells(j).d(6);
%         continue;
        
        
        
%         hsl = jsonData.frames(i).cells(j).d(whichData);
        
%         %legi objective 2:
%         hsl = jsonData.frames(i).cells(j).d(7);
%         legiScale = 20/log(2);
%         hsl = hsl/legiScale - 0.9;
%         threshold = 0.2;
%         selectedStrain = 0;%1=repressore, 0=activator
%         %end legi
        
        threshold = 1500;%change for each sim max HSL in cell
        hsl1 = jsonData.frames(i).cells(j).d(6);
        hsl2 = jsonData.frames(i).cells(j).d(7);
        lacI = jsonData.frames(i).cells(j).d(8);
        
        if 0==thisStrain
            gain = hsl1/threshold;
            if(gain>1) gain=1; end
            cscale = uint8(255*gain);
            cellPaint = [0 cscale 255];  %'c'; 
        else
            gain = hsl2/threshold;
            if(gain>1) gain=1; end
            cscale = uint8(255*gain);
            cellPaint = [cscale 255 0];  %'c'; 
        end
        
        
        x = jsonData.frames(i).cells(j).d(1); xi = round(x) + 1;
        y = jsonData.frames(i).cells(j).d(2); yj = round(y) + 1;
        phi = mod(jsonData.frames(i).cells(j).d(3), 2*pi);
        l = jsonData.frames(i).cells(j).d(4);
        comp = jsonData.frames(i).cells(j).d(5);

        cpoly = makeCell(x,y, phi,l);
        fill(cpoly(1,:), cpoly(2,:), cellPaint); hold on;    
        
        %legi 
%         if (i<=numFrames) && (i > numFrames-colFrameWidth)
%             legiA = jsonData.frames(i).cells(j).d(8);
%             legiR = jsonData.frames(i).cells(j).d(7);
%             colmeansA(xi) = colmeansA(xi) + legiA;
%             colmeansR(xi) = colmeansR(xi) + legiR;
%             colcounter(xi) = colcounter(xi) + 1;
%         end
    end
    hold off;
    axis equal;
    dxy = 5;
    axis([-dxy trapWidth+dxy -dxy trapHeight+dxy]);
    figScale = 10;
    set(gcf,'Position',[100 100 0.6*figScale*trapWidth figScale*trapHeight])

    drawnow;  %will slow down the code
    M(i) = getframe;
    xfile = sprintf('%s\\snapshot_%d.tif', xfolder, i);
    saveas(gcf,xfile);
end
toc




return;
%==========================================================================
%==========================================================================

%%

%LEGI computation;
for col = 1:length(colmeansA)
   if(colcounter(col) ~= 0)
       colmeansA(col) = colmeansA(col)/colcounter(col);
       colmeansR(col) = colmeansR(col)/colcounter(col);
   end
end
figure;  fontSize=30;
legiScale = 20/log(2);
yyaxis left;
plot(colmeansA, 'LineWidth', 3); hold on;
ax = gca;
ax.FontSize = fontSize; 
yyaxis right;
% plot(colmeansR/legiScale, 'LineWidth', 3); hold off;  
plot((smooth(colmeansR/legiScale)), 'LineWidth', 3); hold off;  
legend('a', 'r','FontSize', fontSize, 'Location', 'southeast');
axis([0 105 0.9 1.1]);
xlabel('Trap Position X \mum', 'FontSize', fontSize);  
ylabel('Intracellular Protein (a.u.)', 'FontSize', fontSize);
title({'c_0 = 0.7, H~ = 1, leftWall', 'Intracellular A and R', '(averaged over y-axis of trap)'}, 'FontSize', fontSize);

% delete(poolobj)
display("parpool object deleted.  Done");
% % 
% filename=strcat(num2str(jsonData.timeSinceEpoch), '.avi');
% v = VideoWriter(filename);
% open(v);
% writeVideo(v,M);
% close(v);


%==========================================================================
function cpoly = makeCell(x,y,phi,l) 
% creates a two semicircles from two poles' halves: pA(B)=poleA(B), c=center,t=theta, h=hemicircle
    dtheta = pi/8;
    cc = [x;y];%cell center
    pAc = cc + (l-1)/2 * [cos(phi); sin(phi)];
    pBc = cc - (l-1)/2 * [cos(phi); sin(phi)];
    pAt = (-pi/2+phi):dtheta:(pi/2+phi);
    pBt = (pi/2+phi):dtheta:(3*pi/2+phi);
    pAh = pAc + 0.5*[cos(pAt); sin(pAt)];
    pBh = pBc + 0.5*[cos(pBt); sin(pBt)];
    cpoly = [pAh pBh];
end
%==========================================================================


