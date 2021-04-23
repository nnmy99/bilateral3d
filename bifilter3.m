% BFILTER3 Three dimensional bilateral filter.
%    This function implements 3-D bilateral filter on gray volume.
%    B = bifilter3(A,w,sigma) 
%    A: 3D input volume of type double
%    w: the half-size of the Gaussian filter used in Bilarteral filter.
%    sigma: this is an (1,2) array parameter. Sigma(1) is the spatial
%    standard deviation. Sigma(2) is the intensity standard deviation.

% Developed by Nguyen Ngoc My, Kyoto Institue of Technology, Japan
% Based on 2D Bilateral Filtering of Douglas R. Lanman, Brown University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function B = bifilter3(A,w,sigma)
% Parameters validation.
if isempty(A) || ~exist('A','var')
   error('Input image A is undefined or invalid.');
end
if ~isfloat(A)
   error(['Input image A must be a double precision matrix.']);      
end

% Verify bilateral filter window size.
if ~exist('w','var') || isempty(w) || numel(w) ~= 1 || w < 1
   w = 5; %default w
end
w = ceil(w);

% Verify bilateral filter standard deviations.
if ~exist('sigma','var') || isempty(sigma) || ...
      numel(sigma) ~= 2 || sigma(1) <= 0 || sigma(2) <= 0
   sigma = [3 0.1]; %default sigma
end

% Apply either grayscale or color bilateral filtering.
B = bfltGray(A,w,sigma(1),sigma(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements bilateral filtering for grayscale images.
function B = bfltGray(A,w,sigma_d,sigma_r)
% Pre-compute Gaussian distance weights.
[X,Y,Z] = meshgrid(-w:w,-w:w,-w:w);
G = exp(-(X.^2+Y.^2+Z.^2)/(2*sigma_d^2));
% Create waitbar.
h = waitbar(0,'Applying bilateral filter...');
set(h,'Name','Bilateral Filter Progress');
% Apply bilateral filter.
dim = size(A);
B = zeros(dim);
for i = 1:dim(1)
   for j = 1:dim(2)
        for k = 1:dim(3)
             % Extract local region.
             iMin = max(i-w,1);
             iMax = min(i+w,dim(1));
             jMin = max(j-w,1);
             jMax = min(j+w,dim(2));
             kMin = max(k-w,1);
             kMax = min(k+w,dim(3));
             I = A(iMin:iMax,jMin:jMax,kMin:kMax);
             % Compute Gaussian intensity weights.
             H = exp(-(I-A(i,j,k)).^2/(2*sigma_r^2));
             % Calculate final weight and bilateral filter response.
             F = H.*G((iMin:iMax)-i+w+1,(jMin:jMax)-j+w+1,(kMin:kMax)-k+w+1);
             B(i,j,k) = sum(F(:).*I(:))/sum(F(:));
        end
   end
   waitbar(i/dim(1));
end
% Close waitbar.
close(h);
