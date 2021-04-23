% BFILTER3 Three dimensional bilateral filtering.
%    This function implements 3-D bilateral filtering using
%    the method outlined in:
%
%       C. Tomasi and R. Manduchi. Bilateral Filtering for 
%       Gray and Color Images. In Proceedings of the IEEE 
%       International Conference on Computer Vision, 1998. 
%
%    B = bfilter2(A,W,SIGMA) performs 3-D bilateral filtering
%    for the grayscale image A. A should be a double
%    precision matrix of size NxMx1 with normalized values in
%    the closed interval [0,1]. The half-size of the Gaussian
%    bilateral filter window is defined by W. The standard
%    deviations of the bilateral filter are given by SIGMA,
%    where the spatial-domain standard deviation is given by
%    SIGMA(1) and the intensity-domain standard deviation is
%    given by SIGMA(2).
%
% Douglas R. Lanman, Brown University, September 2006.
% dlanman@brown.edu, http://mesh.brown.edu/dlanman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pre-process input and select appropriate filter.
function B = bifilter3(A,w,sigma)
% Verify that the input image exists and is valid.
if ~exist('A','var') || isempty(A)
   error('Input image A is undefined or invalid.');
end
if ~isfloat(A)
   error(['Input image A must be a double precision matrix.']);      
end
% Verify bilateral filter window size.
if ~exist('w','var') || isempty(w) || ...
      numel(w) ~= 1 || w < 1
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