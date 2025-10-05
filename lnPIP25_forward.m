function lnpip25 = lnPIP25_forward(sic, index, bayes, plotOpt)
% lnPIP25_forward_model  Forward PSM for ln(PIP25) given SIC.
% Written by matt osman (mo549@cam.ac.uk), Oct 2025
%  adapted from original python code by Chung Yan (Crystal) Fu (cyf25@cam.ac.uk)
%
% Inputs:
% sic -- A vector of sea-ice concentrations in fractional units [0..1]
%        (1 x N) or (N x 1). Values 0 and 1 are internally nudged by 1e-4.
%
% index -- 'dino' or 'bras' (selects which calibration/posterior to use)
% bayes -- string denoting Bayesian posterior file to use for the calibration.
%        If empty, the default posterior file for the chosen index is used.
%        Accepted forms are ... 
%           - struct with fields: b0 (M×1), b1 (M×1), sd (M×1)
%           - string/char: path to a 3-col text file [b0 b1 phi]
%           - [] or not supplied: load default file next to this function
% plotOpt -- if true, plots histogram (scalar SIC) or series with HDIs (vector SIC)
%        (optional, defaults to false)
% 
% Output:
% lnpip25 : A set of 1000 ln(PIP25) estimates for each SIC. (N x 1000)
%
% Notes:
%   ln(PIP25) is modeled as:
%       mu = (-log(1/SIC - 1) - b0) ./ b1
%       ln(PIP25) ~ Normal(mu, sd)
%   where sd = sqrt(phi). We draw 1000 posterior predictive samples per SIC
% 
% EXAMPLE; 
% lnpip25 = lnPIP25_forward([0.3], 'dino', [], true);
%
% Dependencies: 
%   if you're plotting, you need cbrewer! 
%   available here: https://uk.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2
%   MBO NOTE: this appears to have been changed to cbrewer2, not cbrewer at
%   some point ... i'm antiquated so i am hard-calling to the 'cbrewer'
%   subdirectory -- you may need to update line96 if you're replacing with cbrewer2
% 
% For details, see:
% Fu, C. Y., Osman, M. B., & Aquino-López, M. A. (2025). Bayesian calibration 
%     for the Arctic sea ice biomarker IP25. Paleoceanography and Paleoclimatology, 
%     40, e2024PA005048. https://doi.org/10.1029/2024PA005048
% -------------------------------------------------------------------------

if nargin < 3; bayes = [];end
if nargin < 4, plotOpt = false; end
if strcmpi(index, 'gatos'); lnpip25 = gatos(); return; end

% vectorise sic + validate
sic = sic(:);
if any(sic < 0 | sic > 1)
    error('SIC values must be within [0,1].');
end

% nudge 0/1 to avoid infinities
sic(sic==0) = 1e-4;
sic(sic==1) = 1 - 1e-4;

% load calibration posterior (b0,b1,sd)
[b0, b1, sd] = load_bayes(index, bayes);

% choose number of predictive draws to return
ndraws = 1000;
% If the posterior has fewer than ndraw draws, resample with replacement
M = numel(b0);
if M < ndraws
    idx = randi(M, ndraws, 1);
    b0s = b0(idx);  b1s = b1(idx);  sds = sd(idx);
else
    % otherwise take first ndraw draws (or sample if you prefer)
    b0s = b0(1:ndraws);  b1s = b1(1:ndraws);  sds = sd(1:ndraws);
end

% compute predictive mean for each (SIC, draw) pair
%     mu = (-log(1/SIC - 1) - b0) ./ b1
%     lnPIP25 ~ N(mu, sd)
N = numel(sic);
mu = (-log(1./sic - 1)) - b0s.'; % (N x ndraw) minus b0 per draw
mu = mu ./ b1s.'; % divide by b1 per draw

% draw posterior predictive noise
eps = randn(N, ndraws) .* sds.';  % add sd per draw

% return posterior predictive samples of ln(PIP25)
lnpip25 = mu + eps;

% =======  optional plotting capabilities ======= 

if ~plotOpt, return; end

% default HDI masses
hdiMass = [0.15, 0.35, 0.55, 0.75, 0.95, 0.9999]; hdiMass = fliplr(hdiMass); 
hdiLabel = ["BaySIC","15% HDI","35% HDI","55% HDI","75% HDI","95% HDI",">95% HDI"]; hdiLabel = fliplr(hdiLabel); 
lnPIP_grid = linspace(-12, 0, 1000);
% grab colours
cd cbrewer/
    warning('off','all'); 
        if strcmp(index,'dino')
        colors = cbrewer('seq','Reds',length(hdiMass)+1); 
        elseif strcmp(index,'bras')
        colors = cbrewer('seq','Blues',length(hdiMass)+1); 
        end
        colors(colors<0) = 0; colors(colors>1) = 1; % colors = flipud(colors); 
    warning('on','all')
cd ../

if N == 1 % plots scalar SIC: histogram + HDIs + MAP ----
   
    x = lnpip25(:);
    [pdfVals, xi] = ksdensity(x, lnPIP_grid, 'Bandwidth', 0.15);
    pdfVals = pdfVals / trapz(xi, pdfVals);
    % compute HDIs
    [lo, hi, bounds] = computeHDI(pdfVals, xi, hdiMass);
    figure('Color','w','Position',[100 100 400 300]); hold on
    for i = 1:numel(bounds)
        % dummy: 
        f(i) = fill([xi(bounds{i}(1):bounds{i}(2)), fliplr(xi(bounds{i}(1):bounds{i}(2)))], ...
                 [pdfVals(bounds{i}(1):bounds{i}(2)), zeros(1, numel(bounds{i}(1):bounds{i}(2)))], ...
                 colors(i,:), 'FaceAlpha',1, 'EdgeColor','none');
    end
    [~, imax] = max(pdfVals);
    f(i+1) = plot([xi(imax), xi(imax)], [0, max(pdfVals)], 'color', colors(i+1,:), 'LineWidth',1.0);
    text(xi(imax)-0.15, max(pdfVals), sprintf('%.2f', xi(imax)), 'HorizontalAlignment','right');
    xlabel(sprintf('ln(PIP_{25}) (%s)', lower(index)));
    ylabel('Probability');
    xlim([-12 0]); ylim([0 max(pdfVals)*1.1]);
    set(gca, 'YTick', []);
    box on
    l1 = legend(gca,f,hdiLabel,'Box','off','Orientation','vertical','fontsize',10); l1.Position = [0.1875    0.5500    0.2200    0.3150]; 

else % vector SIC.. plots series HDI bands + MAP line 
    
    mapVals = nan(N,1);
    loVals  = nan(N,numel(hdiMass));
    hiVals  = nan(N,numel(hdiMass));
    for i = 1:N
        [pdfVals, xi] = ksdensity(lnpip25(i,:), lnPIP_grid, 'Bandwidth', 0.15);
        pdfVals = pdfVals / trapz(xi,pdfVals);
        [~, imax] = max(pdfVals); mapVals(i) = xi(imax);
        [lo, hi] = computeHDI(pdfVals, xi, hdiMass);
        loVals(i,:) = lo(:);
        hiVals(i,:) = hi(:);
    end

    figure('Color','w','Position',[100 100 max(500,40*N) 350]); hold on
    for i = 1:numel(hdiMass)
        xfill = [1:N, fliplr(1:N)];
        yfill = [loVals(:,i)', fliplr(hiVals(:,i)')];
        fill(xfill, yfill, colors(i,:), 'FaceAlpha', 1, 'EdgeColor','none');
    end
    plot(1:N, mapVals, 'color', colors(i+1,:), 'LineWidth',1.2);
    xlabel('input SIC index'); ylabel(sprintf('ln(PIP_{25}) (%s)', lower(index)));
    box on; xlim([1 N]); ylim([-12 0]);

end

end % end of main function

% -------------------------- little helpers ------------------------------
function [b0, b1, sd] = load_bayes(index, bayes)
% locate and load a BaySIC posterior file for the given index
% If 'bayes' input file exists --> use it directly (file or struct), oherwise:
% - find all .txt files in this directory
% - pick the one containing the index ('dino' or 'bras') in its name
% - if multiple matches, prompt user to choose interactively

    index = lower(string(index));

    % if the bayes struct is already provided as a struct :) 
    if isstruct(bayes)
        req = {'b0','b1','sd'};
        if ~all(isfield(bayes, req))
            error('bayes struct must contain fields: b0, b1, sd.');
        end
        b0 = bayes.b0(:);
        b1 = bayes.b1(:);
        sd = bayes.sd(:);
        return
    end

    % if the filename has been manually provided ... 
    if (ischar(bayes) || isstring(bayes)) && ~isempty(bayes)
        M = readmatrix(string(bayes), fileType="text");
        if size(M,2) < 3
            error('posterior file must have 3 columns: b0, b1, phi.');
        end
        b0 = M(:,1); b1 = M(:,2); sd = sqrt(M(:,3));
        return
    end

    % if no input found, throw an error!
    [thisDir, ~, ~] = fileparts(mfilename('fullpath'));
    txtFiles = dir(fullfile(thisDir, '*.txt'));
    if isempty(txtFiles)
        error('no .txt calibration files found in %s. Please run "download_calib_files.m".', thisDir);
    end

    % find matches containing index name
    names = {txtFiles.name};
    matches = contains(lower(names), index);

    if ~any(matches)
        error('no .txt file found containing ''%s'' in its name.', index);
    end

    idxMatch = find(matches);
    if numel(idxMatch) > 1
        % prompt user to choose
        prompt = sprintf('Multiple "%s" posterior files found. Choose one : ', index);
        listStr = names(matches);
        [sel, ok] = listdlg( ...
            'PromptString', prompt, ...
            'SelectionMode', 'single', ...
            'ListString', listStr, ...
            'ListSize', [400, 100]);   % <-- wider (600 px), shallower (150 px)        
        if ~ok
            error('file selection cancelled');
        end
        chosenFile = fullfile(thisDir, names{idxMatch(sel)});
    else
        chosenFile = fullfile(thisDir, names{idxMatch});
    end

    % finally, load selected file!
    M = readmatrix(chosenFile, fileType="text");
    if size(M,2) < 3
        error('Posterior file must have 3 columns: b0, b1, phi.');
    end
    b0 = M(:,1);
    b1 = M(:,2);
    sd = sqrt(M(:,3)); % remaining = coretop SIC posterior draws
end

function lnpip25 = gatos() % meow
    figure('Color','w','Position',[200 200 300 200]); axis off
    text(0.5,0.5, ...
        ["  ∧,,,∧         "; ...
         " (• ˕ •)        "; ...
         "----- U  U ------------"; ...
         "  HAVE A NICE DAY!  "; ...
         "------------------------"], ...
        'FontName','Courier','FontSize',16,'HorizontalAlignment','center','VerticalAlignment','middle');
    lnpip25 = "crystal <3 cats";
end

function [lo, hi, bounds] = computeHDI(pdf, x, masses) % for plotting purposes only! 
    pdf = pdf / trapz(x, pdf);
    sorted = sort(pdf, 'descend');
    cumsumProb = cumsum(sorted) / sum(sorted);
    lo = nan(numel(masses),1);
    hi = nan(numel(masses),1);
    bounds = cell(numel(masses),1);
    for i = 1:numel(masses)
        thr = sorted(find(cumsumProb >= masses(i), 1));
        mask = pdf >= thr;
        idxRange = find(mask);
        lo(i) = x(idxRange(1));
        hi(i) = x(idxRange(end));
        bounds{i} = [idxRange(1), idxRange(end)];
    end
end
