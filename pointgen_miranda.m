%%
% save signature as .png file with white text on black background for ease
% of signature outlining 
%%
% Read image and display it.
I = imread('bwsig.png');
imshow(I)
%%
% Convert image to binary black&white 
BW = im2bw(I);
imshow(BW)
%%
% following finds row/col coord of border pixel (AKA where the signature
% starts) this will be beginning of trace function 
dim = size(BW)
col = round(dim(2)/2)-90;
row = min(find(BW(:,col)))
%%
% |bwtraceboundary| will trace all pixels starting at the point found above
% Call |bwtraceboundary| to trace the boundary from starting pt whose
% coords are defined above (these must be input to the function, must also
% input direction of first step (any cardinal direction) 
boundary = bwtraceboundary(BW,[row, col],'E');
%%
% Show the original image first, then display the trace over this image
% (can check any points where image is not tracing properly) 
imshow(I)
hold on;
plot(boundary(:,2),boundary(:,1),'g','LineWidth',3);
%%
% |bwboundaries| function traces all outlines. 
% To ensure that|bwboundaries| only traces the outside use |imfill| to fill the area
% inside each letter. |bwboundaries| returns a cell array, where each cell
% contains the row/column coordinates for an object in the image.
BW_filled = imfill(BW,'holes');
boundaries = bwboundaries(BW_filled);
%%
% Plot the borders on the original grayscale image using
% the coordinates returned by |bwboundaries| .
for k=1:10
   b = boundaries{k};
   plot(b(:,2),b(:,1),'g','LineWidth',3);
end
