function sic = lnPIP25_predict(ip25, sterol, index, unit, bayes, plotOpt)
% lnPIP25_inverse model  inverse PSM for SIC given ln(PIP25).
% Written by matt osman (mo549@cam.ac.uk), Oct 2025
%  adapted from original python code by Chung Yan (Crystal) Fu (cyf25@cam.ac.uk)
%
% Inputs:
% ip25, sterol -- paired non-negative vectors (both N x 1, same units!)
% index -- input either 'dino' or 'bras' for calibration type
% unit -- input either 'toc' or 'sed' for sediment analysis type
% bayes -- (optional) path to posterior matrix file (.npy or .mat with avgdPDF)
%       .. if empty, auto-discovers by indexing what's in current folder (prompts if multiple!)
% plotOpt -- (optional, default=false) plot PDF (input ip25/sterol scalar values) or series (vector)
%
% Output:
%   sic -- (N x 1000) posterior draws of MAM mean sea-ice concentration in [0,1]
%
% Examples .. 
%   ip25   = rand(20,1) * 0.09; % ~ U(0, 0.09) = surrogate ip25 concentration (0–0.09 mg/g toc)
%   sterol = rand(20,1) * 9.0; % ~ U(0, 9.0) = surrogate sterol concentration (0–9 mg/g toc)
%   sic = lnPIP25_predict(ip25, sterol, 'dino', 'toc'); 
%   sic = lnPIP25_predict(ip25(3), sterol(3), 'bras', 'toc', [], true); % for PDF plotting (scalar only) 
%   sic = lnPIP25_predict(ip25, sterol, 'bras', 'toc', [], true); % for time series plotting 
%
% Dependencies: 
%   needs readNPY.m, available in the "npy-matlab" package 
%       available here: https://github.com/kwikteam/npy-matlab
%   if you're plotting, you need cbrewer! 
%       available here: https://uk.mathworks.com/matlabcentral/fileexchange/58350-cbrewer2
%       MBO NOTE: this appears to have been changed to cbrewer2, not cbrewer at
%       some point ... i'm antiquated so i am hard-calling to the 'cbrewer'
%       subdirectory -- you may need to update line96 if you're replacing with cbrewer2
% 
% For details, see:
% Fu, C. Y., Osman, M. B., & Aquino-López, M. A. (2025). Bayesian calibration 
%     for the Arctic sea ice biomarker IP25. Paleoceanography and Paleoclimatology, 
%     40, e2024PA005048. https://doi.org/10.1029/2024PA005048
% -------------------------------------------------------------------------

addpath(genpath('npy-matlab')) % make sure this folder is in your current directory! 

% some preliminary input checks per usual ... 
if nargin < 5, bayes = []; end
if nargin < 6, plotOpt = false; end

% validate some things!
ip25   = ip25(:);
sterol = sterol(:);
if numel(ip25) ~= numel(sterol), error('ip25 and sterol must have the same length.'); end
if any(ip25<0) || any(sterol<0),  error('ip25/sterol must be non-negative.'); end
index = lower(string(index));
unit  = lower(string(unit));
if ~any(index==["dino","bras"]), error('index must be ''dino'' or ''bras''.'); end
if ~any(unit==["toc","sed"]),    error('unit must be ''toc'' or ''sed''.');   end

% added constants (inferred analytical precision offsets)
min_ip25_toc = 0.01049169444307899;
min_ip25_sed = 6.212970072165696e-05;
min_dino_toc = 0.21311343280837622;
min_dino_sed = 0.0008602308083664112;
min_bras_toc = 0.14608475206034371;
min_bras_sed = 0.0005896700301752964;

% pre-treat concentrations w constants (units + index dependent)
switch unit
    case 'toc'
        ip25 = ip25 + min_ip25_toc;
        sterol = sterol + (index=="dino")*min_dino_toc + (index=="bras")*min_bras_toc;
    case 'sed'
        ip25 = ip25 + min_ip25_sed;
        sterol = sterol + (index=="dino")*min_dino_sed + (index=="bras")*min_bras_sed;
end

% compute lnPIP, boom-diggity
lnPIP = log(ip25 ./ (ip25 + sterol));
N = numel(lnPIP);

% load posterior matrix (reminder: rows: SIC grid (sans 0 and 1), cols: lnPIP grid)
avgdPDF = load_MAM_matrix(index, bayes);

% build grids to match the input matrix (easier than loading in)
[nSIC, nL]   = size(avgdPDF);
lnPIP_grid   = linspace(-12, 0, nL)';
sic_grid = linspace(0, 1, nSIC+2)';
sic_grid = sic_grid(2:end-1); % note! make sic_grid length exactly nSIC (interior of [0,1])

% for each lnPIP, interpolate PDF across SIC and sample 1000 draws
ndraw = 1000;
sic = nan(N, ndraw);
for i = 1:N
    v = lnPIP(i);
    % bracket v on lnPIP_grid
    ui = find(lnPIP_grid >= v, 1, 'first'); if isempty(ui), ui = nL; end
    li = max(1, ui-1);
    lo = lnPIP_grid(li); hi = lnPIP_grid(ui);
    w  = (v - lo) / max(hi - lo, eps);

    upCol = avgdPDF(:, ui);
    loCol = avgdPDF(:, li);
    upPDF = upCol / sum(upCol);
    loPDF = loCol / sum(loCol);
    pdf1d = (1-w)*loPDF + w*upPDF;
    pdf1d = pdf1d / sum(pdf1d);

    % inverse-CDF sampling on discrete grid
    cdf = cumsum(pdf1d); cdf = cdf / cdf(end);
    r = rand(1, ndraw);
    idx = arrayfun(@(u) find(cdf >= u, 1, 'first'), r);
    sic(i,:) = sic_grid(idx);
end

% some optional plotting capabilities
if ~plotOpt, return; end

% default hdi's to plot (these are the same as in Fu et al., 2025)
hdiMass  = [0.15, 0.35, 0.55, 0.75, 0.95, 0.9999]; hdiMass  = fliplr(hdiMass);
hdiLabel = ["BaySIC","15% HDI","35% HDI","55% HDI","75% HDI","95% HDI",">95% HDI"]; hdiLabel = fliplr(hdiLabel);

% load in some nice colours
cd cbrewer/
warning('off','all');
if strcmp(index,'dino')
    colors = cbrewer('seq','Reds',length(hdiMass)+1);
else
    colors = cbrewer('seq','Blues',length(hdiMass)+1);
end
colors(colors<0)=0; colors(colors>1)=1;
warning('on','all');
cd ..

if N == 1 % scalar plotting options ... produces smooth PDF over SIC w specified HDI's + MAP location
    
    % regenerate an analytic PDF (matches the one used for sampling)
    v = lnPIP(1);
    ui = find(lnPIP_grid >= v, 1, 'first'); if isempty(ui), ui = nL; end
    li = max(1, ui-1);
    lo = lnPIP_grid(li); hi = lnPIP_grid(ui);
    w  = (v - lo) / max(hi - lo, eps);

    upCol = avgdPDF(:, ui);
    loCol = avgdPDF(:, li);
    upPDF = upCol / sum(upCol);
    loPDF = loCol / sum(loCol);
    pdfVals = (1-w)*loPDF + w*upPDF;
    pdfVals = pdfVals / sum(pdfVals);

    % HDIs on SIC grid
    [loB, hiB, bounds] = computeHDI(pdfVals, sic_grid, hdiMass);

    figure('Color','w','Position',[100 100 400 300]); hold on
    for iB = 1:numel(bounds)
        seg = bounds{iB}(1):bounds{iB}(2);
        fill([sic_grid(seg).' fliplr(sic_grid(seg).')], ...
             [pdfVals(seg).' zeros(1, numel(seg))], ...
             colors(iB,:), 'FaceAlpha',1, 'EdgeColor','none');
    end
    [~, imax] = max(pdfVals);
    plot([sic_grid(imax), sic_grid(imax)], [0, max(pdfVals)], 'color', colors(iB+1,:), 'LineWidth',1.0);
    text(sic_grid(imax)-0.01, max(pdfVals), sprintf('%.2f', sic_grid(imax)), 'HorizontalAlignment','right');
    xlabel('MAM SIC (%)'); ylabel('Probability');
    xlim([0 1]); ylim([0 max(pdfVals)*1.1]); set(gca,'YTick',[]); box on
    legend(gca, arrayfun(@(k) sprintf('%s', hdiLabel(k)), 1:numel(bounds)+1, 'UniformOutput',false), 'Location', 'best', ...
           'Box','off','Orientation','vertical','FontSize',10);

else % for vector-based plotting  ... gives a reconstructed SIC series w specified HDI's + MAP location
    K = numel(hdiMass);
    hmAsc = sort(hdiMass,'ascend');

    mapVals    = nan(N,1);
    loValsAsc  = nan(N,K);
    hiValsAsc  = nan(N,K);

    % compute HDIs for each point from interpolated PDFs
    for i = 1:N
        v = lnPIP(i);
        ui = find(lnPIP_grid >= v, 1, 'first'); if isempty(ui), ui = nL; end
        li = max(1, ui-1);
        lo = lnPIP_grid(li); hi = lnPIP_grid(ui);
        w  = (v - lo) / max(hi - lo, eps);

        upCol = avgdPDF(:, ui);
        loCol = avgdPDF(:, li);
        upPDF = upCol / sum(upCol);
        loPDF = loCol / sum(loCol);
        pdfVals = (1-w)*loPDF + w*upPDF;
        pdfVals = pdfVals / sum(pdfVals);

        [~, imax] = max(pdfVals);
        mapVals(i) = sic_grid(imax);

        [loA, hiA] = computeHDI(pdfVals, sic_grid, hmAsc);
        loValsAsc(i,:) = loA(:).'; hiValsAsc(i,:) = hiA(:).';
    end

    % create figure and axis handle
    figure('Color','w','Position',[100 100 max(500,30*N) 250]); hold on
    xIdx = 1:N;  f = gobjects(K+1,1);

    for ii = 1:K
        j = find(hmAsc == hdiMass(ii), 1, 'first');
        if j == 1
            f(ii) = fill([xIdx, fliplr(xIdx)], ...
                         [loValsAsc(:,1).', fliplr(hiValsAsc(:,1).')], ...
                         colors(ii,:), 'EdgeColor','none', 'FaceAlpha', 1);
        else
            fill([xIdx, fliplr(xIdx)], ...
                 [loValsAsc(:,j-1).', fliplr(loValsAsc(:,j).')], ...
                 colors(ii,:), 'EdgeColor','none', 'FaceAlpha', 1);
            f(ii) = fill([xIdx, fliplr(xIdx)], ...
                         [hiValsAsc(:,j).',  fliplr(hiValsAsc(:,j-1).')], ...
                         colors(ii,:), 'EdgeColor','none', 'FaceAlpha', 1);
        end
    end

    f(K+1) = plot(xIdx, mapVals, '-', 'Color', colors(K+1,:), 'LineWidth', 1.2);
    xlabel('input data index'); ylabel('MAM SIC (%)');
    box on; xlim([1 N]); ylim([0 1]); xticks(xIdx);
    legend(gca, f, hdiLabel, 'Box','off','Orientation','vertical','FontSize',10,'Location','eastoutside');
end

end % end of main function, woot! 

% ------------------------- helpers --------------------------------------
function avgdPDF = load_MAM_matrix(index, matFile)
    index = lower(string(index));
    if ~isempty(matFile)
        avgdPDF = try_load_matrix(matFile); return
    end
    [thisDir, ~, ~] = fileparts(mfilename('fullpath'));
    % look for files with index name and typical prefix
    cands = [dir(fullfile(thisDir, sprintf('server_MAM_%s_*.npy', index))); ...
             dir(fullfile(thisDir, sprintf('server_MAM_%s_*.mat', index)))];
    if isempty(cands)
        % fallback: any file containing index
        cands = [dir(fullfile(thisDir, sprintf('*%s*.npy', index))); ...
                 dir(fullfile(thisDir, sprintf('*%s*.mat', index)))];
        if isempty(cands)
            error('No MAM matrix found for "%s" in %s.', index, thisDir);
        end
    end
    if numel(cands) > 1
        names = {cands.name};
        [sel, ok] = listdlg('PromptString',sprintf('Multiple "%s" matrices found. Choose one:', index), ...
                            'SelectionMode','single','ListString',names,'ListSize',[600,150]);
        if ~ok, error('File selection cancelled.'); end
        chosen = fullfile(thisDir, names{sel});
    else
        chosen = fullfile(thisDir, cands(1).name);
    end
    avgdPDF = try_load_matrix(chosen);
end

function avgdPDF = try_load_matrix(pathstr)
    [~,~,ext] = fileparts(pathstr);
    switch lower(ext)
        case '.mat'
            S = load(pathstr);
            if isfield(S,'avgdPDF'), avgdPDF = S.avgdPDF;
            else
                fn = fieldnames(S);
                got = false;
                for i=1:numel(fn)
                    if isnumeric(S.(fn{i})), avgdPDF = S.(fn{i}); got = true; break; end
                end
                if ~got, error('No numeric variable found in MAT file "%s".', pathstr); end
            end
        case '.npy'
            if exist('readNPY','file') == 2
                avgdPDF = readNPY(pathstr); % make sure the folder npy-matlab is in your current directory! 

            else
                error('Cannot read NPY file "%s". Add readNPY.m to your path or provide a .mat.', pathstr);
            end
        otherwise
            error('Unsupported matrix file: %s', ext);
    end
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
