function [] = generatePatches(inputImage, labelImage, inputPatchName, labelPatchName, inputPatchDirectory, labelPatchDirectory, patchSize, numberOfPatches)
%GENERATEPATCHES Generates patches of given input-output pair
%   Generates patches for input-output pair at the specified patch size,
%   avoiding the generation of patches with mostly 0 "don't know" pixels

% Create input/labelPatchDirectory for generated patch output. If directory
% already exists, mkdir will produce warning

mkdir(inputPatchDirectory);
mkdir(labelPatchDirectory);

% Compute the "valid" pixels to identify areas of the data where patches
% can be generated by computing the Euclidian distance transform of the
% binary image where known pixels are 0's and unknown pixels are 1's. This
% computes the pixels that are at a distance of at least the image diameter
% away from an unknown pixel. Rotating the image maximizes the number of
% valid patches.

diameter = patchSize/2;

% Pad the input and label image with 0's around the border so generation of
% patches at the border of the image is avoided

inputImage_padded = padarray(inputImage, [1 1], 0);
labelImage_padded = padarray(labelImage, [1 1], 0);

% Compute collection of valid pixels

valid = false(size(label_image)+2);

for degree = 0:5:360
    labelImage_rotated = padarray(imrotate(labelImage_padded, degree, 'crop'), [1 1], 0);
    
    distanceTransform = bwdist(labelImage_rotated == 0, 'chessboard');
    
    valid = valid | imrotate(distanceTransform > diameter, -degree, 'crop');
end

% Define batches and index trackers for generating patches

numberOfPartitions = 25;
batchSize = ceil(numberOfPatches, numberOfPartitions);

maxLabel = double(max(labelImage(:) - 1));

patchesToGenerate = numberOfPatches;
totalPatchIndex = 0;

% PATCH GENERATION

[y,x] = ind2sub(size(valid), find(valid)); % Convert indices of valid pixels to subscripts

if isempty(find(valid, 2))
    warning('Patches with the requested patch size could not be generated',...
        'Request Declined');
end

while patchesToGenerate > 0
    strcat(tostring(patchesToGenerate), " patches left to generate");
   
    if isempty(numberOfPatches)
        numberOfPatches = size(y, 1);
    end
    
    patchesInBatch = min(batchSize, patchesToGenerate);
    patchesToGenerate = patchesToGenerate - patchesInBatch;
    
    validSample = randsample(length(y), patchesInBatch);
    
    for batchPatchIndex = 1:patchesInBatch
       totalPatchIndex = totalPatchIndex + 1;
       
       yPixel = y(validSample(batchPatchIndex));
       xPixel = x(validSample(batchPatchIndex));
       
       try
           rotatedDiameter = ceil(diameter*sqrt(2));
           selectedPatch = labelImage(yPixel - rotatedDiameter:yPixel + rotatedDiameter,...
               xPixel - rotatedDiameter:xPixel + rotatedDiameter);
       catch
           keyboard
       end
       
       % Making 1000 different rotations between 1 and 360 degrees of the
       % patch
       
       randomRotations = (randperm(1000)/1000)*360;
       
       for rotation = randomRotations       
           rotatedPatch = zeros(patchSize, patchSize);
           rotatedPatch = imrotate(selectedPatch, rotation);
           rotatedCenter = round(size(rotatedPatch)/2);  
           rotatedPatch = rotatedPatch(...
               rotatedCenter(1) - patchSize/2:rotatedCenter(1) + patchSize/2-1, ...
               rotatedCenter(2) - patchSize/2:rotatedCenter(2) + patchSize/2-1);
           
           if all(rotatedPatch(:) ~= 0)
               break
           end
       end
       
       selectedPatch = rotatedPatch;
       selectedInput = inputImage(...
           yPixel - rotatedDiameter: yPixel + rotatedDiameter - 1,...
           xPixel - rotatedDiameter: yPixel + rotatedDiameter - 1);
       selectedInput = imrotate(selectedInput, rotation);
       selectedInput = selectedInput(...
           rotatedCenter(1) - patchSize/2:rotatedCenter(1) + patchSize/2-1, ...
           rotatedCenter(2) - patchSize/2:rotatedCenter(2) + patchSize/2-1);
       
       % Save patches
       
       selectedPatch_uint8 = uint8(255*(double(selectedPatch)-1)/maxLabel);
       
       imwrite(selectedPatch_uint8, fullfile(labelPatchDirectory, strcat(labelPatchName, '_', numstr(totalPatchIndex, '%05d'), '.png')));
       imwrite(selectedInput, fullfile(emPatchDirectory, strcat(inputPatchName, '_', numstr(totalPatchIndex, '%05d'), '.png')));
    end    
end

