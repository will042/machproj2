% Point Generator (William)

% Click points on image, tracing signature in the process. Press enter key
% at the end of the initial. Continue clicking through the path while
% tracing the connection between the 1st and 2nd initial. Press enter to
% record the path. Repeat this process through the entire image. This will
% generate 6 sets of x-y data: one for each initial and one for each
% connecting path.

% Points will be saved in a .mat file. This file can be loaded in the
% camgen.m file to generate the cams.


sig = imread('signature_hector.png'); % Must be 2:1 (x:y) Resolution!

signew = imresize(sig, [600 1200]);
imshow(signew)

[x1,y1]=getpts(); % 1st Initial
[xa,ya]=getpts(); % Connection B/W 1st & 2nd Initial
[x2,y2]=getpts(); % 2nd Initial
[xb,yb]=getpts(); % Connection B/W 2st & 3rd Initial
[x3,y3]=getpts(); % 3rd Initial
[xc,yc]=getpts(); % Connection B/W 3rd & 1st Initial




%Store x-y points in structure

points.x1 = x1;
points.xa = xa;
points.x2 = x2;
points.xb = xb;
points.x3 = x3;
points.xc = xc;

points.y1 = y1;
points.ya = ya;
points.y2 = y2;
points.yb = yb;
points.y3 = y3;
points.yc = yc;

% save('har.mat', 'points')
                