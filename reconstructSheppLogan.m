figure % create a new figure
colormap bone %change the colormap to bone
image = imread('Modified_Shepp_Logan.png'); % read the image
subplot(1,3,1) % write the next set of commands to the left third of the picture
imagesc(image); % scale and display image
title('Modified Shepp Logan') % name of image
xlabel('X') % X label
ylabel('Y') % y label
freq = 10; % [1/degree]
thetas = 0:1/freq:360-1/freq; % length of theta vector is 360
subplot(1,3,2) % write the next set of commands to the middle third of the picture

nAP  = length(thetas); %number of Angular Projections
nPP = size(image,1); % number of Parallel Projections

sinogram = zeros(nPP,nAP); % create sinogram

% loop over the 3600 of angles
for i = 1:length(thetas)
   tmpImage      = imrotate(image,-thetas(i),'bilinear','crop'); % rotate image 360 deg
   sinogram(:,i) = sum(tmpImage,2); % fill sinogram
   imagesc(sinogram);% scale and display image
   title('Sinogram')
    xlabel('\alpha')
    ylabel('# parallel projection')
   drawnow % scale and display image continously at each iteration
end
subplot(1,3,3) % write the next set of commands to the right third of the picture
numOfParallelProjections = size(sinogram,1); % obtain size of sinogram
numOfAngularProjections                       = length(thetas); % obtain length of thetas
thetas = (pi/180)*thetas; % convert thetas to radians
BPI = zeros(numOfParallelProjections,numOfParallelProjections);% set up the backprojected image
midindex = floor(numOfParallelProjections/2) + 1; % find the middle index of the projections
[xCoords,yCoords] = meshgrid(ceil(-numOfParallelProjections/2):ceil(numOfParallelProjections/2-1));% set up the coords of the image
rampFilter = [floor(numOfParallelProjections/2):-1:0 1:ceil(numOfParallelProjections/2-1)]'; % set up filter


for i = 1:numOfAngularProjections % loop over each projection
    rotCoords = round(midindex + xCoords*sin(thetas(i)) + yCoords*cos(thetas(i)));  % figure out which projections to add to which spots
    indices   = find((rotCoords > 0) & (rotCoords <= numOfParallelProjections)); % check which coords are in bounds
    newCoords = rotCoords(indices); % figure out which coords in bounds to add to which spots
    filteredProfile = real( ifft( ifftshift( rampFilter.*fftshift(fft(sinogram(:,i)) ) ) ) ); % filter
    BPI(indices) = BPI(indices) + filteredProfile(newCoords)./numOfAngularProjections; % summation
    imagesc(BPI) % scale and display image
    drawnow % scale and display image continously at each iteration
    title('Reconstructed Shepp Logan')
    xlabel('X')
    ylabel('Y')
end