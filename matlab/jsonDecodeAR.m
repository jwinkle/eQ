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
    
    fprintf("Frames scanned; maximum cell ID=%d\n", maxCellNumber);
    
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

    %use force data for each strain type A,B:
    frameData.qScalarA = zeros(numFrames,1);
    frameData.qScalarB = zeros(numFrames,1);
    frameData.forceXA = zeros(numFrames,1);
    frameData.forceYA = zeros(numFrames,1);
    frameData.forceXB = zeros(numFrames,1);
    frameData.forceYB = zeros(numFrames,1);

    % FOR X,Y GRID AND Y-AVERAGED ARRAYS PER FRAME
        global qGrid;
        qGrid =  zeros(trapHeight+1, trapWidth+1, numFrames);
        xFvec =  zeros(trapWidth+1, numFrames);
        xQvec =  zeros(trapWidth+1, numFrames);


    %pretty-print progress counters:    
    displayCounter=1;
    displaySlice = numFrames/10.0; %display progress in 10% intervals
    fprintf('Beginning Frame iteration...\n');

%==========================================================================
    % LOOP PER FRAME:  jsonData.frames(i)
%==========================================================================
    for i= 1:numFrames

        %prettry-print progress
        if(i > displayCounter*displaySlice)
            fprintf('Completed %d percent of frames...\n', round(100*i/numFrames));
            displayCounter = displayCounter + 1;
        end

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
            cos2phiA=0;  sin2phiA=0; fxA=0; fyA=0;
            cos2phiB=0;  sin2phiB=0; fxB=0; fyB=0;
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
            if(0 == thisStrain)
                cos2phiA                     = cos2phiA + cos(2.0*phi);
                sin2phiA                     = sin2phiA + sin(2.0*phi);
                fxA                          = fxA + power(comp, 2.0)*cos(phi);  
                fyA                          = fyA + power(comp, 2.0)*sin(phi);
            else
                cos2phiB                     = cos2phiB + cos(2.0*phi);
                sin2phiB                     = sin2phiB + sin(2.0*phi);
                fxB                          = fxB + power(comp, 2.0)*cos(phi);  
                fyB                          = fyB + power(comp, 2.0)*sin(phi);
            end

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
                numA = numCells(i) - frameData.strainRatio(i);
                numB = frameData.strainRatio(i);

                frameData.strainRatio(i)    = frameData.strainRatio(i)/numCells(i);
                frameData.avg_comp(i)       = frameData.avg_comp(i)/numCells(i);
                frameData.avg_length(i)     = frameData.avg_length(i)/numCells(i);
                
                cos2phiA                     = cos2phiA/numA; 
                sin2phiA                     = sin2phiA/numA; 
                frameData.qScalarA(i)        = sqrt(cos2phiA*cos2phiA + sin2phiA*sin2phiA);
                frameData.forceXA(i)         = fxA;
                frameData.forceYA(i)         = fyA;
                    cos2phiB                     = cos2phiB/numB; 
                    sin2phiB                     = sin2phiB/numB; 
                    frameData.qScalarB(i)        = sqrt(cos2phiB*cos2phiB + sin2phiB*sin2phiB);
                    frameData.forceXB(i)         = fxB;
                    frameData.forceYB(i)         = fyB;
                %average columns (but only if counter > 0)
                    ix = find(xFi>0); 
                    xF(ix) = xF(ix) ./ xFi(ix);
                    xFvec(:,i) = xF;

    end%end frames
    fprintf('done.\n');
    
%==========================================================================
    % AVERAGE DATA OVER LAST n FRAMES:
%==========================================================================
%  <deleted 1/20/20

    return;
    
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

smdataA = frameData.qScalarA;
smdataB = frameData.qScalarB;
% smdataA = smooth(frameData.qScalarA);
% smdataB = smooth(frameData.qScalarB);

    plot(smdataA, 'b', 'LineWidth', 3); hold on;
    plot(smdataB, 'g', 'LineWidth', 3); 
    plot(frameData.strainRatio, 'k', 'LineWidth', 5); hold off;

axis([-Inf Inf 0 1.1]);
ax = gca;
ax.FontSize = FS; 
xlabel('Time (mins)'); 
ylabel('Q, STRAIN FRACTION');
legend({'Q(STRAIN A)', 'Q(STRAIN B)', 'FRACTION B'}, 'Location', 'southeast');
title({'STRAIN FRACTION OSCILLATIONS', 'and ORDER PARAMETER Q'}, ...
    'FontSize', FS);

figure(10); clf;
% subplot(1,2,1);


yyaxis right
%     plot(frameData.qScalar, 'LineWidth', 3); hold off;
%     plot([frameData.qScalarA; frameData.qScalarB], 'LineWidth', 3); hold off;
% 
% axis([-Inf Inf 0 1]);
% ylabel('Whole-Trap q-Scalar');
% legend({'Strain Ratio','QA', 'QB'}, 'Location', 'northeast');
% title({'STRAIN RATIO, Q-SCALAR', ...
% %     'MEAN DIVISION LENGTH = 4.2 \mum'}, ...
%     'l_d = 4.2--> 2.8 \mum'}, ...
%     'FontSize', FS);
% 
% 
% subplot(1,2,2);
% yyaxis left;
% plot(frameData.avg_comp); hold on;
% plot([abs(frameData.forceX) abs(frameData.forceY)], 'LineWidth', 5); 
smspan = 20;
smxdataA = smooth(abs(frameData.forceXA), smspan);
smxdataB = smooth(abs(frameData.forceXB), smspan);
maxsdata = max([smxdataA' smxdataB']);

plot(smxdataA, 'b-', 'LineWidth', 5); hold on;
plot(smxdataB, 'g-', 'LineWidth', 5); hold on;
% plot(abs(frameData.forceY), 'b-', 'LineWidth', 5);  hold off;
axis([-Inf Inf 0 2*maxsdata]);
xlabel('Time (mins)');  
ylabel('Mean X component Force (a.u.)');
% yyaxis right
% plot(frameData.avg_comp, 'LineWidth', 3); hold off;
% axis([-Inf Inf 0 Inf]);
% ylabel('Mean Cell-spring Compression');
% legend({'[ForceX]^2', '[ForceY]^2','Comp'}, 'Location', 'east');

yyaxis left;
    plot(frameData.strainRatio, 'k', 'LineWidth', 5); hold on;

axis([-Inf Inf 0 1]);
ax = gca;
ax.FontSize = FS; 
xlabel('Time (mins)','FontSize', FS);  ylabel('Strain Ratio','FontSize', FS);

legend({'[ForceX]', '[ForceY]'}, 'Location', 'east');
title('MEAN SPRING FORCE', ...
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

    %pretty-print progress counters:    
    displayCounter=1;
    displaySlice = numFrames/10.0; %display progress in 10% intervals
    fprintf('Beginning Frame iteration for movie file...\n');

for i=1:numFrames
% for i=numFrames-colFrameWidth:numFrames
% for i=1:100
% for i=790:815
% parfor i=1:numFrames

        %prettry-print progress
        if(i > displayCounter*displaySlice)
            fprintf('Completed %d percent of frames...\n', round(100*i/numFrames));
            displayCounter = displayCounter + 1;
        end
%     if(0 == mod(i,10)) fprintf("\ri=%d", i); end

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


