function [meanSIC, monthsUsed, meta] = calc_meanSIC(sic_climo, varargin)
% calc_meanSIC  Mean SIC over the 3 months ending at the first post-maximum decrease.
% Written by Matt Osman (mo549@cam.ac.uk), Oct 2025
%  adapted from original Python by C.Y. (Crystal) Fu
%
% summary
%   For a 12-month SIC climatology, computes the mean SIC of the three months
%   ending at the first post-maximum decrease in SIC (after rounding to 5%).
%   Works on a single site (vector) or a grid. In grid mode, it finds the
%   nearest grid cell to the imput site location, and (if that cell is constant/invalid),
%   automatically switches to the nearest SIC-variable cell.
%
% inputs
%   sic_climo : EITHER a 12-element vector (fractional SIC in [0,1]; NaN ok)
%               OR a 3-D grid with one dimension of length 12.
%               The function auto-reorients grid inputs to lon x lat x 12.
%
%   varargin  : (grid mode only; order matters)
%               site_lat   – scalar (degrees)
%               site_lon   – scalar (degrees; 0–360deg or −180…180deg both ok)
%               all_lat    – vector Ny x 1 or mesh Ny x Nx (degrees)
%               all_lon    – vector 1 x Nx or mesh Ny x Nx (degrees; 0–360° ok)
%
% outputs
%   meanSIC     : scalar — mean SIC over months [m−2, m−1, m], where m is the
%                 month of first post-maximum decrease (1..12; Jan=1).
%   monthsUsed  : 1 x 3 integer months (1..12), chronological.
%   meta        : struct with diagnostics (grid mode includes coords/indices)
%       .flag_no_clear_decrease : true if multi-max scan fell back to first max
%       .usedNearestVariable    : true if the nearest cell was constant and
%                                 a variable cell was selected instead
%       .iy, .ix                : selected grid indices (lat row, lon col)
%       .lat, .lon              : selected grid coordinates (degrees)
%       .distance               : km to selected grid cell
%       .reoriented             : true if any dimension reordering was applied
%       .movedMonthsTo3rd       : true if month dim moved to 3rd
%       .swappedLonLat          : true if lat x lon swapped to lon x lat
%       .sic                    : 12 x 1 SIC series used for the calculation
%
% some notes on function's behaviour / user rules 
%   * vlues must be proportions in [0,1]; NaNs allowed (ignored in means).
%   * logic for 'first SIC decrease' is performed on SIC rounded to the nearest 0.05.
%   * "First decrease" uses a 24-month (doubled) stack to scan forward from each
%     annual maximum to the first local minimum (last non-increasing step)....
%     among multiple maxima, choose the one with the largest drop.
%   * grid inputs are normalised to lon x lat x 12 (months in 3rd dimension). If
%     provided as,e .g., lat x lon x 12 or 12 x lat x lon or otherwise, they are 
%     swapped automatically.
%   * mixed longitude conventions are supported (0–360deg vs −180…180deg); distances
%     use a wrapped Δlon in [−/pi,/pi], so nearest-cell selection is robust.
%
% usage (varargin examples)
%   % 1) Vector mode (single site series)
%   [meanSIC, monthsUsed] = calc_meanSIC(sic_climo);  % sic_climo is 12 x 1 or 1 x 12 climatology of   
%                           mean monthly SIC values for a site
%
%   % 2) Grid mode (sicGrid == 3D matrix of SIC values, lon x lat x 12), with auto-reorientation if needed
%   % site_lat, site_lon, all_lat, all_lon must be provided in that order..
%   % site_lon can be between -180 to 360; all_lat and all_lon match matrix dimensions of input sicGrid
%   [meanSIC, monthsUsed, meta] = calc_meanSIC(sicGrid, site_lat, site_lon, all_lat, all_lon);
%
% other examples (using provided 3D SIC climatology data in example_sic_data)
%   % Example A — Vector from a grid cell (convert % to fraction)
%   load example_sic_data/curvilinear_sic_climo.mat % provides SIC (lon x lat x 12 or similar), lat, lon
%   sic_vec = squeeze(SIC(114,174,:)) ./ 100; % 12 x 1, convert from %
%   [meanSIC, monthsUsed] = calc_meanSIC(sic_vec);
%
%   % Example B — Full grid with varargin, example must search for nearest SIC grid with monthly variabilty
%   load example_sic_data/tripolar_sic_climo.mat  % SIC, lat (Ny), lon (Nx) or meshes
%   site_lat = lat(74);
%   site_lon = lon(359);  % works whether lon is 0–360 or −180..180
%   [meanSIC, monthsUsed, meta] = calc_meanSIC(SIC./100, site_lat, site_lon, lat, lon);
%   % Inspect chosen cell and series:
%   meta.sic, meta.lat, meta.lon, meta.distance
%
% NOTE SOME IMPORTANT WARNINGS / ERRORS
%   • Errors if all values are NaN or annual mean is constant (i.e., no sea ice or no change across year).
%   • Errors if no dimension equals 12 or grid sizes don’t match coord sizes.
%   • Warns (flag_no_clear_decrease=true) if a unique decrease month isn’t clear.
%
% see also (i.e., this function helps generate the first argument needed for):
%   lnPIP25_forward_model.m (forward PSM )  
%
% Reference (please cite when using this fucntion)
%   Fu, C.Y., Osman, M.B., & Aquino-Lopez, M.A. (2025).
%   Bayesian calibration for the Arctic sea ice biomarker IP25.
%   Paleoceanography and Paleoclimatology, 40, e2024PA005048. https://doi.org/10.1029/2024PA005048
% -------------------------------------------------------------------------

% Ensure SIC is a proportion
if ~isnumeric(sic_climo) || any((~isnan(sic_climo(:))) & (sic_climo(:)<0 | sic_climo(:)>1))
    error('sic_climo must be numeric with all non-NaN values in [0,1].');
end

meta = struct('reoriented',false,'movedMonthsTo3rd',false,'swappedLonLat',false);

% --- detect month dimension (size 12) ---
sz  = size(sic_climo);
nd  = ndims(sic_climo);
d12 = find(sz==12,1);
if isempty(d12), error('Input must include a dimension of length 12.'); end

% Count non-singleton dims
ns = nnz(sz~=1);

% if vector input for sic_climo, make sure oriented as a 12 x 1 input
if ns==1
    % Make sure it's 12 x 1 for vector path
    if d12~=1
        sic_climo = permute(sic_climo,[d12 setdiff(1:nd,d12,'stable')]);
        meta.reoriented = true;
    end
    sic_climo = reshape(sic_climo,12,1);
    nd = ndims(sic_climo);
    
% if 3D matrix input for sic_climo, make sure oriented as a lon x lat x 12 (fixes dimensions to input lon and lat vectors)  
elseif nd>=3
    % We need site_lat, site_lon, all_lat, all_lon to enforce lon/lat orientation
    if numel(varargin) < 4
        error('Grid mode requires: site_lat, site_lon, all_lat, all_lon.');
    end
    site_lat = varargin{1};
    site_lon = varargin{2};
    all_lat  = varargin{3};
    all_lon  = varargin{4};

    % 1) Move months to 3rd dimension (… x … x 12)
    if d12~=3
        sic_climo = permute(sic_climo,[setdiff(1:nd,d12,'stable') d12]);
        meta.reoriented = true;
        meta.movedMonthsTo3rd = true;
    end

    % enforce lon x lat x 12 using coord sizes
    if isvector(all_lon) && isvector(all_lat)
        nx = numel(all_lon);  ny = numel(all_lat);
    else
        % meshgrid style: [lat x lon]
        [ny, nx] = size(all_lon);
    end

    % if data are lat x lon x 12, swap to lon x lat x 12
    if size(sic_climo,1)==ny && size(sic_climo,2)==nx
        sic_climo = permute(sic_climo,[2 1 3]);  % -> lon x lat x 12
        meta.reoriented = true;
        meta.swappedLonLat = true;
    elseif size(sic_climo,1)==nx && size(sic_climo,2)==ny
        % already lon x lat x 12
    else
        error('SIC grid (%d x %d) does not match lon/lat sizes (%d x %d).', size(sic_climo,1),size(sic_climo,2),nx,ny);
    end
    nd = ndims(sic_climo);
else
    error('Expect a 12-element vector or a grid with a 12-month axis.');
end

% ===================== VECTOR MODE =====================
if nd == 2 && numel(sic_climo) == 12
    v = sic_climo(:);
    [meanSIC, monthsUsed, meta_vec] = handle_vector(v);
    meta.flag_no_clear_decrease = meta_vec.flag_no_clear_decrease;
    return

% ===================== GRID MODE (lon x lat x 12) =======================
elseif nd == 3
    % We already validated we have site_lat, site_lon, all_lat, all_lon above
    site_lat = varargin{1};
    site_lon = varargin{2};
    all_lat  = varargin{3};
    all_lon  = varargin{4};

    % Normalize lat/lon inputs to LAT(lat x lon), LON(lat x lon)
    [Ny, Nx] = deal(size(sic_climo,2), size(sic_climo,1)); % lat count, lon count
    [LAT, LON] = normalize_grid_latlon(all_lat, all_lon, Ny, Nx);

    % Find nearest grid cell (Haversine on degrees arrays)
    dist_all = haversine_deg(site_lat, site_lon, LAT(:), LON(:));
    [~, p0] = min(dist_all);
    [iy0, ix0] = ind2sub([Ny, Nx], p0);  % iy: lat row, ix: lon col
    v0 = squeeze(sic_climo(ix0, iy0, :));  % 12x1
    distClosest = dist_all(p0); % grab the closest distance

    usedNearestVariable = false;
    % flag_no_clear_decrease = false;
    iy = iy0; ix = ix0; v = v0;

    if is_constant_series(v0) % find nearest grid cell w variable sea ice..
        isVar = false(Ny*Nx,1);
        for p = 1:Ny*Nx
            [iy_p, ix_p] = ind2sub([Ny, Nx], p);
            vv = squeeze(sic_climo(ix_p, iy_p, :));
            isVar(p) = ~is_constant_series(vv(:));
        end
        if any(isVar)
            dist_var = dist_all; dist_var(~isVar) = Inf;
            [~, p1] = min(dist_var);
            [iy, ix] = ind2sub([Ny, Nx], p1);
            v = squeeze(sic_climo(ix, iy, :));
            usedNearestVariable = true;
            distClosest = dist_var(p1); % grab the distance to closest variable sea ice
        else
            meanSIC = NaN; monthsUsed = [NaN NaN NaN];
            meta = struct('flag_no_clear_decrease',true, ...
                          'usedNearestVariable',false, ...
                          'iy',iy0,'ix',ix0,'lat',LAT(iy0,ix0),'lon',LON(iy0,ix0), ...
                          'site_lat',site_lat,'site_lon',site_lon, ...
                          'reoriented',meta.reoriented, ...
                          'movedMonthsTo3rd',meta.movedMonthsTo3rd, ...
                          'swappedLonLat',meta.swappedLonLat,...
                          'distance',distClosest, ...
                          'sic',v(:));
            warning('SIC is constant at all grid cells; cannot identify first decrease.');
            return
        end
    end

    [meanSIC, monthsUsed, submeta] = handle_vector(v(:));
    flag_no_clear_decrease = submeta.flag_no_clear_decrease;

    meta = struct('flag_no_clear_decrease',flag_no_clear_decrease, ...
                  'usedNearestVariable',usedNearestVariable, ...
                  'iy',iy,'ix',ix,'lat',LAT(iy,ix),'lon',LON(iy,ix), ...
                  'site_lat',site_lat,'site_lon',site_lon, ...
                  'reoriented',meta.reoriented, ...
                  'movedMonthsTo3rd',meta.movedMonthsTo3rd, ...
                  'swappedLonLat',meta.swappedLonLat, ...
                  'distance',distClosest, ...
                  'sic',v(:));

    if flag_no_clear_decrease
        warning(['No clear month of first SIC decrease detected at site; ', ...
                 '3D input preferred for fallback to variable grid.']);
    end
    return
else
    error('Unsupported input shape for sic_climo.');
end
end


%% ----- helper functions -----

function [meanSIC, months, meta] = handle_vector(v)
    meta.flag_no_clear_decrease = false;
    if all(isnan(v)) || numel(v)~=12
        meanSIC = NaN; months = [NaN NaN NaN]; meta.flag_no_clear_decrease = true; return;
    end

    annual_sic = mean(v,'omitnan');
    if isnan(annual_sic) || annual_sic == 0
        error('SIC must not be all NaN or all zeros. Provide a 12 x lat x lon field instead.');
    end

    vq = round(v(:)/0.05)*0.05; % quantize to 5%
    vs = [vq; vq]; % stacked for wrap
    vmax = max(vq,[],'omitnan');
    Imax = find(vq == vmax);

    if numel(Imax) == 12
        error('SIC is constant throughout the year. Provide a 12 x lat x lon field instead.');
    elseif numel(Imax) == 1
        decIdx = Imax(1);
    else
        K = numel(Imax);
        Imin = nan(K,1);
        for k = 1:K
            start = Imax(k);
            m = 1;
            next_sic = vs(start + m);
            if next_sic ~= vs(start)
                while next_sic <= vs(start + m - 1)
                    m = m + 1;
                    next_sic = vs(start + m);
                end
                Imin(k) = start + m - 1;
            end
        end
        valid = ~isnan(Imin);
        if ~any(valid)
            decIdx = Imax(1);
            meta.flag_no_clear_decrease = true;
        else
            diffs = vs(Imin(valid)) - vs(Imax(valid));
            [~, j] = min(diffs);
            decIdx = Imax(valid); decIdx = decIdx(j);
        end
    end

    m3 = decIdx; m2 = wrap12(m3-1); m1 = wrap12(m3-2);
    months = [m1 m2 m3];
    meanSIC = mean(v(months),'omitnan');
end

function tf = is_constant_series(v)
    v=v(:); v=v(~isnan(v));
    tf = numel(v)<2 || (max(v)-min(v))<0.05;
end

function idx = wrap12(k), idx = mod(k-1,12)+1; end

function [LAT,LON] = normalize_grid_latlon(all_lat,all_lon,Ny,Nx)
    % Return LAT/LON as [Ny x Nx] meshes (lat rows, lon cols)
    if isvector(all_lat) && isvector(all_lon)
        [LON,LAT] = meshgrid(all_lon(:).', all_lat(:));
    else
        LAT = all_lat; LON = all_lon;
    end
    if ~isequal(size(LAT),[Ny Nx]) || ~isequal(size(LON),[Ny Nx])
        error('LAT/LON shape mismatch with SIC grid (lat x lon = %d x %d).', Ny, Nx);
    end
end

function d = haversine_deg(lat1,lon1,lat2,lon2)
    % the great-circle distance is invariant to whether lon inputs are in 0–360 or −180…180!
    R = 6371; 
    lat1r = deg2rad(lat1); lon1r = deg2rad(lon1);
    lat2r = deg2rad(lat2); lon2r = deg2rad(lon2);
    dlat = lat2r - lat1r;
    dlon = mod(lon2r - lon1r + pi, 2*pi) - pi;
    a = sin(dlat/2).^2 + cos(lat1r).*cos(lat2r).*sin(dlon/2).^2;
    d = 2*R*asin(min(1,sqrt(a)));
end
