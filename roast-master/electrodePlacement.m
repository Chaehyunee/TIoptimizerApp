function hdrInfo = electrodePlacement(P1,P2,T2,elecNeeded,options,uniTag)
% hdrInfo = electrodePlacement(P1,P2,T2,elecNeeded,options,uniTag)
%
% Place electrodes on the scalp surface. options.elecPara contains all the options
% info for each electrode.
% 
% (c) Yu (Andy) Huang, Parra Lab at CCNY
% yhuang16@citymail.cuny.edu
% April 2018

[dirname,baseFilename] = fileparts(P1);
if isempty(dirname), dirname = pwd; end
[~,baseFilenameRasRSPD] = fileparts(P2);
if isempty(T2)
    baseFilenameRasRSPD = [baseFilenameRasRSPD '_T1orT2'];
else
    baseFilenameRasRSPD = [baseFilenameRasRSPD '_T1andT2'];
end

elecPara = options.elecPara;

indP = elecPara(1).indP;
indN = elecPara(1).indN;
indC = elecPara(1).indC;

% %% can be any non-ras head (to be consistent with user-provided coordinates)
% Enforcing RAS in the first step starting from ROAST v3.0
% landmarks_original = getLandmarks(P2,T2);
% 
% [perm,iperm,isFlipInner,isFlipOutter] = how2getRAS(landmarks_original);

template = load_untouch_nii([dirname filesep baseFilenameRasRSPD '_masks.nii']);
% Load the scalp mask; template is used for saving the results with the same header info as the input
pixdim = template.hdr.dime.pixdim(2:4);
dim = size(template.img);
v2w = [template.hdr.hist.srow_x;template.hdr.hist.srow_y;template.hdr.hist.srow_z;0 0 0 1];
hdrInfo = struct('pixdim',pixdim,'dim',dim,'v2w',v2w);
% keep the header info for use later
% scalp_original = template.img==5;
% scalp = changeOrientationVolume(scalp_original,perm,isFlipInner);
scalp = template.img==5;

if ~isempty(indP) || ~isempty(indN)
%     landmarks = changeOrientationPointCloud(landmarks_original,perm,isFlipInner,size(scalp));
    landmarks = getLandmarks(P2,T2);
end

if ~isempty(indC)
    fid = fopen([dirname filesep baseFilename '_customLocations']);
    capInfo_C = textscan(fid,'%s %f %f %f');
    fclose(fid);
    elecLoc_C = cell2mat(capInfo_C(2:4));
    % update customized elec coords according to options of isNonRAS, resamp, and zeroPad
    if options.isNonRAS
        [elecLoc_C,perm] = convertToRASpointCloud(P1,elecLoc_C);
    else
        perm = [1 2 3];
    end
    if options.resamp
        data = load_untouch_nii(P1);
        temp = data.hdr.dime.pixdim(2:4);
        temp = temp(perm);
        elecLoc_C = elecLoc_C.*repmat(temp,size(elecLoc_C,1),1);
    end
    if options.zeroPad>0
        elecLoc_C = elecLoc_C + options.zeroPad;
    end
    
    elecLoc_C = elecLoc_C(indC,:);
%     elecLoc_C = changeOrientationPointCloud(elecLoc_C,perm,isFlipInner,size(scalp));
end

scalp_surface = mask2EdgePointCloud(scalp,'erode',ones(3,3,3));

%% fit cap position on the individual's head
if ~isempty(indP)
   switch lower(elecPara(1).capType)
       case {'1020','1010','1005'}
           load('cap1005FullWithExtra.mat','capInfo');
           isBiosemi = 0;
           isEGI = 0;
       case 'biosemi'
           load('capBioSemiFullWithExtra.mat','capInfo');
           isBiosemi = 1;
           isEGI = 0;
       case 'egi'
           load('capEGIfull.mat','capInfo');
           isBiosemi = 0;
           isEGI = 1;
   end
   [electrode_coord_P,center_P]= fitCap2individual(scalp,scalp_surface,landmarks,P2,capInfo,indP,isBiosemi,isEGI);
else
    electrode_coord_P = []; center_P = [];
end

if ~isempty(indN)
    if any(landmarks(5:6,3)<=0)
        error('MRI does not cover the neck, so cannot place electrodes on the neck. Consider using ''zeroPadding'' option to extend the input MRI.');
    else
        [electrode_coord_N,center_N] = placeNeckElec(scalp,scalp_surface,landmarks,indN);
    end
else
    electrode_coord_N = []; center_N = [];
end

if ~isempty(indC)
    [~,indOnScalpSurf] = map2Points(elecLoc_C,scalp_surface,'closest');
    electrode_coord_C = scalp_surface(indOnScalpSurf,:);
else
    electrode_coord_C = [];
end

%% head clean up for placing electrodes
[scalp_clean,scalp_filled] = cleanScalp(scalp,scalp_surface);
temp1 = scalp_filled(:,:,[1 end]);
temp2 = scalp_filled(:,[1 end],:);
temp3 = scalp_filled([1 end],:,:);
if any([temp1(:);temp2(:);temp3(:)])
    warning('Scalp touches image boundary. Electrodes may go out of image boundary. ROAST can continue but results may not be accurate. It is recommended that you expand the input MRI by specifying the ''zeroPadding'' option.');
end    
scalp_clean_surface = mask2EdgePointCloud(scalp_clean,'erode',ones(3,3,3));

disp('calculating gel amount for each electrode...')
if ~isempty(indP)
    %     elec_range_P = zeros(size(electrode_coord_P,1),100);
    [~,indOnScalpSurf] = project2ClosestSurfacePoints(electrode_coord_P,scalp_clean_surface,center_P);
    elec_range_P = indOnScalpSurf(1:100,:);
    %     for i=1:size(elec_range_P,1)
    %         elec_range_P(i,:) = indOnScalpSurf(1:100,i);
    % Get some points on the scalp surface that are close to the exact
    % location of each electrode for the calculation of local normal vector
    % for each electrode in the following step
    %     end
else
    elec_range_P = [];
end

if ~isempty(indN)
    % Get local scalp points for neck electrodes
    %     elec_range_N = zeros(size(electrode_coord_N,1),100);
    [~,indOnScalpSurf] = project2ClosestSurfacePoints(electrode_coord_N,scalp_clean_surface,center_N);
    elec_range_N = indOnScalpSurf(1:100,:);
    %     for i=1:size(elec_range_N,1)
    %         elec_range_N(i,:) = indOnScalpSurf(1:100,i);
    %     end
else
    elec_range_N = [];
end

if ~isempty(indC)
    [~,elec_range_C] = map2Points(electrode_coord_C,scalp_clean_surface,'closer',100);
    % [~,elec_range_C] = map2Points(electrode_coord_C,scalp_surface,'closer',100);
else
    elec_range_C = [];
end

% elecPool = cat(1,elecPool_P,elecPool_N,elecPool_C);
electrode_coord = cat(1,electrode_coord_P,electrode_coord_N,electrode_coord_C);
% electrode_center = cat(1,repmat(center_P,size(electrode_coord_P,1),1),...
%     repmat(center_N,size(electrode_coord_N,1),1),repmat(center_C,size(electrode_coord_C,1),1));
elec_range = cat(1,elec_range_P',elec_range_N',elec_range_C');

%% placing and model the electrodes
resolution = mean(pixdim);
% mean() here to handle anisotropic resolution; ugly. Maybe just
% resample MRI to isotropic in the very beginning?
[elec_C,gel_C] = placeAndModelElectrodes(electrode_coord,elec_range,scalp_clean_surface,scalp_filled,elecNeeded,elecPara,resolution,0,uniTag);

%% generate final results (elec and gel masks, and their coordinate ranges)
disp('constructing electrode and gel volume to be exported...')
% for i = 1:length(elec_C)
%     if ~isempty(elec_C{i})
%         elec_C{i} = changeOrientationPointCloud(elec_C{i},iperm,isFlipOutter,size(scalp_original));
%         gel_C{i} = changeOrientationPointCloud(gel_C{i},iperm,isFlipOutter,size(scalp_original));
%     end
% end

% [volume_elec,volume_elecLabel] = generateElecMask(elec_C,size(scalp_original),elecNeeded,1);
% [volume_gel,volume_gelLabel] = generateElecMask(gel_C,size(scalp_original),elecNeeded,0);
% volume_elec_C = generateElecMask(elec_C,size(scalp_original),elecNeeded,1);
% volume_gel_C = generateElecMask(gel_C,size(scalp_original),elecNeeded,0);
volume_elec_C = generateElecMask(elec_C,size(scalp),elecNeeded,1);
volume_gel_C = generateElecMask(gel_C,size(scalp),elecNeeded,0);

disp('final clean-up...')
volume_elec = volume_elec_C>0;
volume_gel = volume_gel_C>0;
volume_gel = xor(volume_gel,volume_gel & volume_elec); % remove the gel that overlaps with the electrode
for i=1:6
    volume_tissue = template.img==i;
    volume_gel = xor(volume_gel,volume_gel & volume_tissue);
end % remove the gel that goes into other tissue masks

% disp('saving the results...')
template.fileprefix = [dirname filesep baseFilename '_' uniTag '_mask_elec'];
template.hdr.hist.descrip = 'electrode mask';
% template.img = uint8(volume_elec)*255;
template.img = uint8(volume_elec_C.*volume_elec);
save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_mask_elec.nii']);
template.fileprefix = [dirname filesep baseFilename '_' uniTag '_mask_gel'];
template.hdr.hist.descrip = 'gel mask';
% template.img = uint8(volume_gel)*255;
template.img = uint8(volume_gel_C.*volume_gel);
save_untouch_nii(template,[dirname filesep baseFilename '_' uniTag '_mask_gel.nii']);

% save([dirname filesep baseFilename '_' uniTag '_labelVol.mat'],'volume_elecLabel','volume_gelLabel');
[~,baseFilenameRasRSPD] = fileparts(P2);
if ~exist([dirname filesep baseFilenameRasRSPD '_header.mat'],'file')
    save([dirname filesep baseFilenameRasRSPD '_header.mat'],'hdrInfo');
end