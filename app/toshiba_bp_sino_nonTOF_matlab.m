function bpimg = toshiba_bp_sino_nonTOF(scanner, sino, image_size, voxel_size, number_of_radial_bins)
% single bed only!
warning('FOR SINGLE BED ONLY!');
disp('sino must be in [nrads x nplanes x nangs]. Use Toshiba original raw data format');

%scanner = buildPET('toshiba');

nring = scanner.getNumberOfCrystalRings();

number_of_planes = nring * nring;
nrads = number_of_radial_bins;

xp = scanner.getDefaultSinogramCrystalPairs();
nrads_default = scanner.system_parms.number_of_projections_per_angle;
xp = reshape(xp, 2, nrads_default, size(xp,2) / nrads_default);
nang = size(xp, 3);

if nrads > nrads_default
    error('FOV is too big!');
end
%
xp = xp(:, ceil(nrads_default/2) + (- floor(nrads/2) : floor(nrads/2)), :);
xp = reshape(xp, 2, nrads * nang);

% Toshiba's sinogram plane ordering
ring_pair = 0;
for n = 1 : nring
%   
    if n==1
        ii = [1 : nring; 1 : nring];
        ring_pair = ii;
    else
        i1 = [1 : (nring-n+1); n : (n+(nring-n+1)-1)];
        i2 = flipud(i1);        
        ring_pair = [ring_pair, i1, i2];
    end
end

%
seg=nrads * number_of_planes * 64;
nseg = nang / 64;

i0 = int16([0 0 0 0 0]');
lm = repmat(i0, 1, seg);

bpimg = zeros(image_size);

for m = 1 : nseg
    fprintf('processing seg#%d ...\n', m);

    gid = (0 : (seg-1)) + (m-1)*seg;
    phi = floor(gid / (nrads * number_of_planes));
    r0 = mod(gid, nrads*number_of_planes);
    nr = floor(r0 / nrads);
    rad = mod(r0, nrads);
    nx = rad + phi * nrads;
    
    lm(1,:) = xp(1, nx+1) - 1;
    lm(3,:) = xp(2, nx+1) - 1;
    lm(2,:) = ring_pair(2, nr+1) - 1; % switch !!!
    lm(4,:) = ring_pair(1, nr+1) - 1;
    
    prjs = sino(gid + 1);
    bp = scanner.doListModeBackProjectionNonTOF(prjs, image_size, voxel_size, lm);
    
    % sum up!
    bpimg = bpimg + bp;
    
end
