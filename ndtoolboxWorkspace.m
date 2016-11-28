clear;
close all;

% AM Last Modified 05/07/2016
% Mostly based on the runND.m function, which is called in demoND.m
% function [machine, outputMisc, outputConf, outputROC, outputData] = runND(traindataNorOri,...
%     validdataNorOri, testdataNorOri, validdataAbOri, testdataAbOri, NDtype, paramsPreset, paramsPlot)
%
% Please first download the NDtoolbox from the CHI lab page http://www.robots.ox.ac.uk/~davidc/publications_NDtool.php
% Also, run demoND.m to understand what it does.
% Then load the data that was created in spectrawriteupAMmodifed.m script
% (uncomment the corresponding lines)

optbyAUC = false; % whether optimise parameter values using my partial AUC method.

% Choose which figures to plot.
paramsPlot.plotoriftr = true; % plot both original and scaled features; I have only coded for 2D and 4D features.    
paramsPlot.plotprojftr = true; % plot projected features on 2D space; if nftrs <=2, my code will automatically set plotprojftr = false.
paramsPlot.plotROC = true; % set plotROC = true if you want an AUC value.
paramsPlot.plotROC_Pr = false; % should set it to "false", because I decided not to use CalROC_Pr for my ROC output. See my note "ROC2011.txt"
paramsPlot.plotoutput = true; % plot ND output (i.e, novelty scores).
paramsPlot.plotcontour = true; % plot 2D contour of output.
paramsPlot.ploterr = true; % plot errors on validation data, used by MinErr_thr.
if ~isfield(paramsPlot, 'plotoriftr')  
    paramsPlot.plotoriftr = false;
end

if ~isfield(paramsPlot, 'plotprojftr')  
    paramsPlot.plotprojftr = false;
end

if ~isfield(paramsPlot, 'plotROC')  
    paramsPlot.plotROC = false;
end

if ~isfield(paramsPlot, 'plotROC_Pr')  
    paramsPlot.plotROC_Pr = false;
end

if ~isfield(paramsPlot, 'plotoutput')  
    paramsPlot.plotoutput = false;
end

if ~isfield(paramsPlot, 'plotcontour')  
    paramsPlot.plotcontour = false;
end

if ~isfield(paramsPlot, 'ploterr')  
    paramsPlot.ploterr = false;
end

plotoriftr = paramsPlot.plotoriftr;
plotprojftr = paramsPlot.plotprojftr;
plotROC = paramsPlot.plotROC;
plotROC_Pr = paramsPlot.plotROC_Pr;
plotoutput = paramsPlot.plotoutput;
plotcontour = paramsPlot.plotcontour;
ploterr = paramsPlot.ploterr;
    
paramsPreset =[];

% Load data
% whichdata = 'banana_n200200s1';
% data = load(whichdata);
% Pump 18
% whichdata = 'pump18RecPartial12';
% load('C:\Users\engs1602\research\data\pumpWiimoteND\pump18RecPartial12');
% Pump 32
whichdata = 'WDT32RecZ_HpfMaDelNdsHamm256FlZPFreq56_2_64dB';
load('C:\Users\engs1602\research\meetings\smallGroup\20161104ManandharAnalysisCodeUpdate\data\WDT32RecZ_HpfMaDelNdsHamm256FlZPFreq56_2_64dB');
fprintf('\nLoading data set %s...\n', whichdata);
alldataOri = data.x; % numdata by numftrs
classlabels = data.y; % numdata by 1, class labels = 1, 2
isnor = classlabels == 1; % regard class 1 as normal.
isab = ~isnor;

[traindataNorOri, testdataNorOri, validdataNorOri, validdataAbOri, testdataAbOri] = splitData(alldataOri, isab);

%% Scale data, according to traindataNor. 
[testdataNor, traindataNor] = scaleData(testdataNorOri, traindataNorOri);
[validdataNor] = scaleData(validdataNorOri, traindataNorOri);
[validdataAb] = scaleData(validdataAbOri, traindataNorOri);
[testdataAb] = scaleData(testdataAbOri, traindataNorOri);

%% Re-pack all the data together, for the purpose of plotting.
alldataOri = [traindataNorOri; validdataNorOri; testdataNorOri; validdataAbOri; testdataAbOri]; % alldataOri are original data before scaling.
alldata = [traindataNor; validdataNor; testdataNor; validdataAb; testdataAb]; % alldata are scaled data
isab = [false(size(traindataNor,1),1); false(size(validdataNor, 1),1); false(size(testdataNor,1),1); true(size(validdataAb,1),1); true(size(testdataAb,1),1)];

% NDtype = 'parzen';
NDtype = 'nn';
% NDtype = 'svmtax';
% NDtype = 'svmsch';
% NDtype = 'kde';
% NDtype = 'pca';
%% Set parameters of the ND method
[paramsND] = setNDParams(NDtype); % set parameters; for svm, paramsND.nu = 0.1; paramsND.d = 1.5;

% Use preset parameter values if any.
if ~isempty(paramsPreset) && strncmpi(NDtype, 'svm', 3)
    if isfield(paramsPreset, 'd')
        paramsND.d = paramsPreset.d;
    end
    
    if isfield(paramsPreset, 'nu')
        paramsND.nu = paramsPreset.nu;
    end
end

%% Use validation data (both normal and abnormal) to find optimal param values by AUC. See NDarti.
% paramsROC are used by MaxAUC_SVMdnu, CalROC and CalROC_Pr.
if optbyAUC || plotROC || plotROC_Pr
    paramsROC.nthresh = 100; % number of thresholds for calculating ROC
    % AM Modified
    paramsROC.whichAUC = 'whole'; % = 'whole', or 'partial'
%     paramsROC.whichAUC = 'partial'; % = 'whole', or 'partial'
%     paramsROC.Pfpmin = 0; % only used for paramsROC.whichAUC = 'partial';
%     paramsROC.Pfpmax = 0.3; % only used for paramsROC.whichAUC = 'partial';
end

if optbyAUC
    if strncmpi(NDtype, 'svm', 3)        
        paramsAUC.NDtype = NDtype;
        paramsAUC.dRange  = [1, 4];        
        paramsAUC.nuRange = [0.1, 0.5];
        [maxAUC, optd, optnu] = MaxAUC_SVMdnu(traindataNor, validdataNor, validdataAb, paramsAUC, paramsROC);
        % Use the opimal values of (d, nu) by AUC.
        paramsND.d = optd;
        paramsND.nu = optnu;
        fprintf('Have chosen optimal parameter values by AUC: d = %.2f, nu = %.2f at AUC = %.2f using validation data.\n', optd, optnu, maxAUC);       
    else
        fprintf('Message from "%s": you have chosen to optimise parameters using AUC by optbyAUC = true,\n', mfilename);
        fprintf('\t but I have only programmed this for NDtype = "svmSch" and "svmTax", and your NDtype = %s.\n', NDtype);        
        fprintf('\t Therefore, default parameters (set by function "SetNDParams") have been used.\n');
        fprintf('\n');
    end
end

% Display values of (d, nu) if svm
if strncmpi(NDtype, 'svm', 3)
    fprintf('SVM: have chosen (d, nu) = (%.2f, %.2f)\n', paramsND.d, paramsND.nu);
end

% % Train the machine.
% fprintf('Train the machine...\n');
% machine = train_parzen(traindataNor, paramsND);
% % Use the trained machine to get output
% fprintf('Use the trained machine to get output...\n');
% outputTrain = out_parzen(traindataNor, machine);
% outputTest_nor = out_parzen(testdataNor, machine); % ntest_nor by 1.
% outputTest_ab  = out_parzen(testdataAb,  machine); % ntest_nor by 1.
% outputData.outputTrain = outputTrain;
% outputData.outputTest_nor = outputTest_nor;
% outputData.outputTest_ab = outputTest_ab;
%% Run the chosen novelty detection method, using test data.
trainfunc = ['train_', NDtype];
outfunc = ['out_', NDtype];
% Train the machine.
fprintf('Train the machine...\n');
machine = feval(trainfunc, traindataNor, paramsND); % machine is a structure.
% Use the trained machine to get output
fprintf('Use the trained machine to get output...\n');
if strcmp(NDtype, 'nn') || strcmp(NDtype, 'som')
    outputTrain = feval(outfunc, traindataNor, machine, 1); % ntrain by 1.
else
    outputTrain = feval(outfunc, traindataNor, machine); % ntrain by 1.
end
outputTest_nor = feval(outfunc, testdataNor, machine); % ntest_nor by 1.
outputTest_ab  = feval(outfunc, testdataAb,  machine); % ntest_nor by 1.
% output = [outputTrain; outputTest_nor; outputTest_ab]; % output is 500 x 1.
outputData.outputTrain = outputTrain;
outputData.outputTest_nor = outputTest_nor;
outputData.outputTest_ab = outputTest_ab;

figure;
subplot(2,1,1);plot(cat(1,outputTrain,outputTest_nor,outputTest_ab));grid on;
subplot(2,1,2);plot(cat(1,true(size(traindataNor,1),1),true(size(testdataNor,1),1),false(size(testdataAb,1),1)));grid on;

%% Set threshold for each NDtype, using validation data if possible
switch lower(NDtype)
    case 'svmsch'
        thresh = 0;  % 0 is the boundary by support vectors.
    case 'svmtax'
        thresh = machine.R; % R is the radius.
    case {'gmm', 'parzen', 'kde', 'gpoc'} % use validation data to set threshold.
        % thresh = 0.01;  % pdf > 0.01 is deemed normal; an empirical value.        
        [minErr, thresh] = minErr_thr(machine, validdataNor, validdataAb, NDtype, ploterr);
        fprintf('Have chosen optimal threshold = %.5f for NDtype = %s at min Err = %.3f, using validation data.\n', thresh, NDtype, minErr);        
    case{'dist', 'nn', 'kmeans', 'pca', 'kpca'}
        thresh = max(outputTrain); % threshold is the maximum output of training normal data.
    case{'som'}
        thresh = mean(outputTrain); % threshold is the maximum output of training normal data
    otherwise
        error('Unknown NDtype for setting ND threshold - please add code.');        
end
fprintf('Threshold = %0.3f for NDtype = %s.\n', thresh, NDtype);
outputMisc.thresh = thresh;

%% Assign class label on the output of test data: 0 for normal, 1 for abnormal. (0-1 coding)
[predTrain] = assignCls(NDtype, outputTrain, thresh);
[predTest_nor] = assignCls(NDtype, outputTest_nor, thresh); % predTest_nor is same size as outputTest_nor.
[predTest_ab] = assignCls(NDtype, outputTest_ab, thresh);
predTest = [predTest_nor; predTest_ab];

%% Compute confusion matrix of test data.
tarTest_nor = zeros(size(outputTest_nor, 1), 1);
tarTest_ab = ones(size(outputTest_ab, 1), 1);
tarTest = [tarTest_nor; tarTest_ab];
[conf, rate] = confmat(predTest, tarTest); % predTest and tarTest are 0-1 coding.
fprintf('\n');
disp('Confusion matrix using test data is:');
disp(conf);
accuracy = (conf(1,1) + conf(2,2)) / sum(conf(:)); % accuracy = rate(1)
sensitivity = conf(2,2)/(conf(2,1)+conf(2,2));
specificity = conf(1,1)/(conf(1,1)+conf(1,2));
fprate = conf(1,2)/(conf(1,1)+conf(1,2)); % false positive rate
fnrate = conf(2,1)/(conf(2,1)+conf(2,2)); % false negative rate
fprintf('[accuracy, sensitivity, specificity] = [%0.2f, %0.2f, %0.2f]%%\n', accuracy*100, sensitivity*100, specificity*100);
fprintf('[FP rate, FN rate] = [%0.2f, %0.2f]%%\n', fprate*100, fnrate*100);
fprintf('accuracy = total correct / total data = %d / %d = %0.2f%%\n', rate(2), sum(conf(:)), accuracy*100);
fprintf('\n');
outputConf.conf = conf;
outputConf.accuracy = accuracy;
outputConf.sensitivity = sensitivity;
outputConf.specificity = specificity;
outputConf.fprate = fprate;
outputConf.fnrate = fnrate;

%% Set class symbols for plotting
normalSymbol = 'g.';
abnormalSymbol = 'r.';
normalSymbolB = 'b.';
abnormalSymbolB = 'm.';
pointsize = 14;

%% Sanity check: to see whether one needs to plot projected features. Do this before plotoriftr.
nftrs = size(traindataNor, 2);
if plotprojftr
    if nftrs <= 2
        plotprojftr = false;
        plotoriftr = true;
        fprintf('Message from "%s": you have chosen to plot the projected 2D features by setting plotprojftr = true,\n', mfilename);
        fprintf('\t but the number of features = %d which is <= 2; therefore no need to reduce dimensions!\n', nftrs);        
        fprintf('\t I will plot the original features for you.\n');
        fprintf('\n');
    end
end

%% Take a look at the original and scaled data
if plotoriftr    % plot original features before scaling.
    switch nftrs
        case 2
            figure('name', 'ori');
            % Plot original features before scaling.
            subplot(121);
            hold on;
            plot(alldataOri(~isab, 1), alldataOri(~isab, 2), normalSymbol, 'markersize', pointsize);
            plot(alldataOri(isab, 1), alldataOri(isab, 2), abnormalSymbol, 'markersize', pointsize);
            title('Original 2D features');
            hold off;
            legend('nor', 'ab', 'location', 'best');
            % Plot scaled features, by simply replacing "alldataOri" with "alldata".
            subplot(122); % same legend as subplot(121), so I have omitted it here.
            hold on;
            plot(alldata(~isab, 1), alldata(~isab, 2), normalSymbol, 'markersize', pointsize);
            plot(alldata(isab, 1), alldata(isab, 2), abnormalSymbol, 'markersize', pointsize);
            title('Scaled 2D features');
            hold off;            
        case 4
            figure('name', 'ori');
            % Plot original features before scaling.
            subplot(121)
            hold on;
            plot(alldataOri(~isab, 1), alldataOri(~isab, 2), normalSymbol, 'markersize', pointsize);
            plot(alldataOri(isab, 1), alldataOri(isab, 2), abnormalSymbol, 'markersize', pointsize);
            plot(alldataOri(~isab, 3), alldataOri(~isab, 4), normalSymbolB, 'markersize', pointsize);
            plot(alldataOri(isab, 3), alldataOri(isab, 4), abnormalSymbolB, 'markersize', pointsize);
            title('Original 4D features');
            hold off;
            legend('(ftr1,ftr2), nor', '(ftr1,ftr2), ab', '(ftr3,ftr4), nor', '(ftr3,ftr4), ab', 'location', 'best');
            % Plot scaled features, by simply replacing "alldataOri" with "alldata".
            subplot(122); % same legend as subplot(121), so I have omitted it here.
            hold on;
            plot(alldata(~isab, 1), alldata(~isab, 2), normalSymbol, 'markersize', pointsize);
            plot(alldata(isab, 1), alldata(isab, 2), abnormalSymbol, 'markersize', pointsize);
            plot(alldata(~isab, 3), alldata(~isab, 4), normalSymbolB, 'markersize', pointsize);
            plot(alldata(isab, 3), alldata(isab, 4), abnormalSymbolB, 'markersize', pointsize);
            title('Scaled 4D features');
        otherwise
            fprintf('Message from "%s": You have chosen to plot the original features.\n', mfilename);
            fprintf('\t Your number of features = %d, but I have only coded for 2D and 4D...\n', nftrs);
            fprintf('\n');
    end
end

%% Take a look at the projected 2D data.
% The dimestion reduction here is only for the purpose of visualisation. I have not used the projected features in my ND code.
if plotprojftr % plot projected features in 2D space. I have ensured that nftrs > 2 in code above.
    if nftrs <= 2 % Just to be sure again that nftrs > 2.        
        fprintf('Message from "%s": You have chosen to plot the projected 2D features,\n', mfilename);
        fprintf('\t but the number of features = %d which is <= 2; therefore no need to reduce dimensions!\n', nftrs);
        fprintf('\n');
    else % now we are absolutely sure that nftrs > 2.        
        whichproject = 'pca'; % pca is a netlab function, see Netlab p.230
        % whichproject = 'neuroscale'; % Netlab p.264. demns1        
        switch lower(whichproject)
            case 'pca'
                [pcacoef, pcavec] = pca(alldata); % call Netlab function pca.
                mu = mean(alldata);
                projdata = (alldata - ones(size(alldata, 1), 1)*mu) * pcavec(:, 1:2);
            case 'neuroscale' % code is from Netlab's demns1.
                ncentres = 10;
                input_dim = size(alldata, 2);
                output_dim = 2; % project onto 2D space.
                net = rbf(input_dim, ncentres, output_dim, 'tps', 'neuroscale');
                % First row controls shadow targets, second row controls rbfsetbf
                options(1, :) = foptions;
                options(2, :) = foptions;
                options(1, 1) = 1;
                options(1, 2) = 1e-2;
                options(1, 3) = 1e-2;
                options(1, 6) = 1;    % Switch on PCA initialisation
                options(1, 14) = 60;
                options(2, 1) = -1;   % Switch off all warnings
                options(2, 5) = 1;
                options(2, 14) = 10;
                net2 = rbftrain(net, options, alldata);
                projdata = rbffwd(net2, alldata);
        end                
        figure('name', 'proj');
        hold on;
        plot(projdata(~isab,1), projdata(~isab, 2), normalSymbol, 'MarkerSize', pointsize);
        plot(projdata(isab,1), projdata(isab, 2), abnormalSymbol, 'MarkerSize', pointsize);
        hold off;
        legend('nor', 'ab', 'location', 'best');
    end
end

%% Calculate ROC curve: Pfp, Ptp and fscore are nthresh x 1, valAUC is a scalar.
% The "CalROC" method is the preferred method to the "CalROC_Pr" method.
if plotROC    
    fprintf('Message from "%s": you have set plotROC = true. I will plot ROC curve produced by function "CalROC".\n', mfilename);    
    fprintf('\n');
    
    % Calculate ROC curve: Pfp, Ptp and fscore are nthresh x 1, valAUC is a scalar.
    [Pfp, Ptp, fscore, valAUC] = calROC(machine, testdataNor, testdataAb, paramsROC);
    outputROC.Pfp = Pfp;
    outputROC.Ptp = Ptp;
    outputROC.fscore = fscore;
    outputROC.valAUC = valAUC;
    
    if isempty(valAUC) % function CalROC or CalROC_Pr did not get a valid AUC value.
        fprintf('Warning from "%s": you have chosen plotROC = true, but function "CalROC" did not get a valid AUC value!\n', mfilename);
        fprintf('\t Therefore no ROC curve is plotted.\n');
        fprintf('\n');
    else
        fprintf('Output by "CalROC": You have chosen %s AUC method, %sAUC = %0.2f using test data.\n', paramsROC.whichAUC, paramsROC.whichAUC, valAUC);
        fprintf('\n');
        figure('name', 'ROC');
        titlestr = sprintf('ROC, AUC=%.2f', valAUC);
        plot(Pfp(:,1), Ptp(:,1), '.-'); axis([0 1 0 1]); title(titlestr); xlabel('Pfp'); ylabel('Ptp');
    end
end

%% Plot ND score and contour.
if plotoutput
    figure('name', 'out');
    hold on;
    plot(outputTrain, normalSymbol, 'MarkerSize', pointsize);
    plot(outputTest_nor, normalSymbolB, 'MarkerSize', pointsize);
    plot(outputTest_ab, abnormalSymbol, 'MarkerSize', pointsize);
    xlimits = get(gca,'XLim');
    line([xlimits(1) xlimits(2)], [thresh thresh],'Color','k','LineStyle','-.'); 
    hold off;
    % Set my own ylim, because Matlab sometimes leave too much margin on y-axis.
    ymin = min([outputTrain;outputTest_nor;outputTest_ab]);
    ymax = max([outputTrain;outputTest_nor;outputTest_ab]);
    yrange = ymax - ymin;
    ymargin = 0.05 * yrange; % 0.05 is a good margin. Change 0.05 to 0.1 for wider margin on y-axis.
    ylim([ymin-ymargin, ymax+ymargin]);
    % Add legend and labels.
    legend('train nor', 'test nor', 'test ab', 'thresh', 'location', 'best');    
    xlabel('samples');
    ylabel('ND output');
    title(NDtype);
end

if plotcontour && size(alldataOri,2)==2  % plot 2D contour of output.
    switch lower(NDtype)
        case 'parzen'
            [fh] = plotParzen(machine, testdataNor, testdataAb);        
        case {'gmm', 'kde'}
            paramsCont.thresh = thresh;
            plotGMM2Dcont(machine, testdataNor, testdataAb, paramsCont);
            title(gca, 'Contour of gmm output on test data');
        case {'svmsch', 'svmtax'}
            plotSVM2Dcont(machine, testdataNor, testdataAb, [], []);
        case {'nn', 'kmeans', 'som', 'pca','kpca'}
            plot2DDecisionBound(NDtype, machine, testdataNor, testdataAb, outfunc, thresh);
        case {'gpoc'}
            paramsCont.thresh = thresh;
            plotGP2Dcont(machine, testdataNor, testdataAb, validdataNor, validdataAb, paramsCont);
    end
end

% paramsPreset = [];
% 
% paramsPlot.plotoriftr = true;   % plot both original and scaled features; I have only coded for 2D and 4D features.
% paramsPlot.plotprojftr = true;  % plot projected features on 2D space; if nftrs <=2, my code will automatically set plotprojftr = false.
% paramsPlot.plotROC = true;      % set plotROC = true if you want an AUC value.
% paramsPlot.plotROC_Pr = false;  % should set it to "false", because I decided not to use CalROC_Pr for my ROC output. See my note "ROC2011.txt"
% paramsPlot.plotoutput = true;   % plot ND output (i.e, novelty scores).
% paramsPlot.plotcontour = true;  % plot 2D contour of output.
% 
% NDtype = {'parzen'};
% % NDtype = {'svmTax'};
% % NDtype = {'pca'};
% 
% load('C:\Users\engs1602\research\data\pumpWiimoteND\pump32RecPartial');
% alldataOri = data.x; % numdata by numftrs
% classlabels = data.y; % numdata by 1, class labels = 1, 2
% 
% %% Separate normal data from abnormal data (if any).
% % Produce quantitative results (2 classes are presented)
% isnor = classlabels == 1; % regard class 1 as normal.
% isab = ~isnor;
% normaldataOri = alldataOri(isnor);
% abnormaldataOri = alldataOri(~isnor);
% [traindataNorOri, testdataNorOri, validdataNorOri, validdataAbOri, testdataAbOri] = splitData(alldataOri, isab);
% 
% [machine, outputMisc, outputConf, outputROC, outputData] = runND(traindataNorOri,...
%     validdataNorOri, testdataNorOri, validdataAbOri, testdataAbOri, lower(NDtype{1}), paramsPreset, paramsPlot);        