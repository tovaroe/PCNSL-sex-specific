function [H_count, DAB_count] = CD3_counting(image, showfigs)
%CD3_COUNTING Returns the numbers of CD3-stained and Haematoxylin-stained cells    
%
%   [H_count, DAB_count] = CD3_counting() lets the user choose the path to
%   the image file and shows figures for each segmentation step. It retuns
%   the number of Haematoxylin-stained cells (H_count) and DAB-stained
%   cells (DAB_count).
%
%   [...] = CD3_counting(IM) performs the cell counting on an image IM and
%   does not show figures for each segmentation step.
%
%   [...] = CD3_counting(IM, SHOWFIGS) shows figures for each segmentation
%   step if SHOWFIGS =/= 0.
%
%
% Copyright (c) 2019 Thomas Roetzer, MedUni Vienna
% thomas.roetzer@meduniwien.ac.at
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%

if (nargin==1)
    showfigs=false;
end
if (nargin==0)
    showfigs = true;
    [file, path] = uigetfile('*.tif*');
    image = imread([path filesep file]);
end

%% Segment cells
res = segment_cells(image, ...
    'stain', 'H DAB 2', ...
    'channels', {1 2}, ...
    'phansalkar_k', [0.15 0.15], ...
    'otsuThresholdFactor', [1 1.1], ...
    'globalThreshold', [220 200], ...
    'cellSizeMin', [0 0], ...
    'survivalNeighbors', [2 2]);

%% Cell counting

H = res{1};
DAB = res{2};

% DAB -------------------------------------------------------------------

H_labeled = labelmatrix(bwconncomp(H));
unique_overlap_objects = unique(H_labeled(DAB & H_labeled~=0));
overlap_objects = ismember(H_labeled, unique_overlap_objects);
DAB = DAB | overlap_objects;
DAB = watershed_segmentation(DAB,1.8);

STATS = regionprops(DAB, 'PixelIdxList', 'Area');
STATS = STATS([STATS.Area] >= 6);
DAB_count = length(STATS);


if showfigs
    DAB = zeros(size(DAB));
    DAB(cat(1,STATS.PixelIdxList)) = 1;
    boundaries = bwboundaries(DAB);
    DAB_boundary_image = image;
    for k=1:length(boundaries)
        boundary = boundaries{k};
        for j=1:length(boundary)
            px = boundary(j,:);
            DAB_boundary_image(px(1), px(2), :)=255;
        end
    end
    figure(20), imshow(DAB_boundary_image), title('DAB boundaries');
end

% Haematoxylin ----------------------------------------------------------
H = H | DAB;
H = watershed_segmentation(H);
STATS = regionprops(H, 'PixelIdxList', 'Area');
STATS = STATS([STATS.Area] >= 4);
H_count = length(STATS);


if showfigs
    H = zeros(size(H));
    H(cat(1, STATS.PixelIdxList)) = 1;
    boundaries = bwboundaries(H);
    H_boundary_image = image;
    for k=1:length(boundaries)
        boundary = boundaries{k};
        for j=1:length(boundary)
            px = boundary(j,:);
            H_boundary_image(px(1), px(2), :)=255;
        end
    end
    figure(21), imshow(H_boundary_image), title('Haematoxylin boundaries');
end