% Cam Profile Generator (Coded by William)
% (SolidWorks Export & Pressure Angle coded by Hector)

% Specify points for signature to load below.
% Define cam base radiuses & z-follower displacement.
% Press F5 to run.


load('har.mat') % Specify point set to load here.

yDPI = 600; % Shoudn't have to change this if image was properly resized in pointgen.m


% Specify Base Radiuses & Max Z Follower Displacement (pen lift height)

Lh = .25;   % Max Z Follower Displacement (inches)
Brx = 2;  % X Cam Base Radius (inches)
Bry = 2;  % Y Cam Base Radius (inches)
Brz = 2;  % Z Cam Base Radius (inches)



%% Do not modify below this line

% Flip Coordinates of Y Values to Normal Coordinate System

points.y1 = yDPI - points.y1;
points.y2 = yDPI - points.y2;
points.y3 = yDPI - points.y3;
points.ya = yDPI - points.ya;
points.yb = yDPI - points.yb;
points.yc = yDPI - points.yc;


% Define span for each initial, and regions were pen is lifted

range.r1 = 1:numel(points.x1);
range.ra = 1:numel(points.xa);
range.r2 = 1:numel(points.x2);
range.rb = 1:numel(points.xb);
range.r3 = 1:numel(points.x3);
range.rc = 1:numel(points.xc);

span.s1 = 1 : length(range.r1);
span.sa = length(range.r1) : length(range.r1)+length(range.ra);
span.s2 = length(range.r1)+length(range.ra) : length(range.r1)+length(range.ra)+length(range.r2);
span.sb = length(range.r1)+length(range.ra)+length(range.r2) : length(range.r1)+length(range.ra)+length(range.r2)+length(range.rb);
span.s3 = length(range.r1)+length(range.ra)+length(range.r2)+length(range.rb) : length(range.r1)+length(range.ra)+length(range.r2)+length(range.rb)+length(range.r3);
span.sc = length(range.r1)+length(range.ra)+length(range.r2)+length(range.rb)+length(range.r3) : length(range.r1)+length(range.ra)+length(range.r2)+length(range.rb)+length(range.r3)+length(range.rc);


% Define a vector of all displacement points for X & Y

x_disp_points = [points.x1; points.xa; points.x2; points.xb; points.x3; points.xc];
y_disp_points = [points.y1; points.ya; points.y2; points.yb; points.y3; points.yc];

range_full = 1:numel(x_disp_points);    % Length of vectors
range_full = range_full';               % Convert to column vector



% Build Smoothing Functions for Data

% Generate Fourier Series fits for each initial and lifted region

fsfit.x1fit = fit( range.r1', points.x1,'fourier6');    % 1st Initial
fsfit.y1fit = fit( range.r1', points.y1,'fourier6');    
fsfit.xafit = fit( range.ra', points.xa,'fourier2');    % Lift
fsfit.yafit = fit( range.ra', points.ya,'fourier2');
fsfit.x2fit = fit( range.r2', points.x2,'fourier8');    % 2nd Initial
fsfit.y2fit = fit( range.r2', points.y2,'fourier8');
fsfit.xbfit = fit( range.rb', points.xb,'fourier2');    % Lift
fsfit.ybfit = fit( range.rb', points.yb,'fourier2');
fsfit.x3fit = fit( range.r3', points.x3,'fourier5');    % 3rd Initial
fsfit.y3fit = fit( range.r3', points.y3,'fourier5');
fsfit.xcfit = fit( range.rc', points.xc,'fourier3');    % Lift
fsfit.ycfit = fit( range.rc', points.yc,'fourier3');


% Combine the FS fits for each section into a single column vector

x_disp_fsfit = [fsfit.x1fit(range.r1)
                fsfit.xafit(range.ra)
                fsfit.x2fit(range.r2)
                fsfit.xbfit(range.rb)
                fsfit.x3fit(range.r3)
                fsfit.xcfit(range.rc)];
     
y_disp_fsfit = [fsfit.y1fit(range.r1)
                fsfit.yafit(range.ra)
                fsfit.y2fit(range.r2)
                fsfit.ybfit(range.rb)
                fsfit.y3fit(range.r3)
                fsfit.ycfit(range.rc)];


% Refit the data set with a smoothing spline, to eliminate discontinuities
% between each fit
            
x_fit = fit(range_full, x_disp_fsfit, 'smoothingspline', 'SmoothingParam',.02);
y_fit = fit(range_full, y_disp_fsfit, 'smoothingspline', 'SmoothingParam',.02);


x_disp_smoothfit = x_fit(range_full)/yDPI;  % Used to convert to inches!
y_disp_smoothfit = y_fit(range_full)/yDPI;
% 
% x_disp_smoothfit(end) = x_disp_smoothfit(1);
% y_disp_smoothfit(end) = y_disp_smoothfit(1);


% Plot Points & Fit Data

figure(1)
    plot(x_disp_smoothfit(span.s1),y_disp_smoothfit(span.s1),...
         x_disp_smoothfit(span.sa),y_disp_smoothfit(span.sa),...
         x_disp_smoothfit(span.s2),y_disp_smoothfit(span.s2),...
         x_disp_smoothfit(span.sb),y_disp_smoothfit(span.sb),...
         x_disp_smoothfit(span.s3),y_disp_smoothfit(span.s3),...
         x_disp_smoothfit(span.sc),y_disp_smoothfit(span.sc),'LineWidth',1); hold
         scatter(x_disp_points/yDPI, y_disp_points/yDPI,25,'.','LineWidth',2); hold
    grid
    axis([0 2 0 1])
    title('Input Points & Generated Parametric Fit')
    xlabel('(inches)')
    ylabel('(inches)')
    legend('1st Initial', 'Lift', '2nd Initial', 'Lift', '3rd Initial', 'Lift')



%% Generate Cam Profiles

theta = linspace(0,2*pi,length(range_full))';    
thetad = theta/pi*180;


% Define Angular Displacement for Lifted Zones
angdisp.sa = theta(span.sa(end))-theta(span.sa(1));
angdisp.sb = theta(span.sb(end))-theta(span.sb(1));
angdisp.sc = theta(span.sc(end))-theta(span.sc(1));


% Harmonic Rise & Fall Function for z-axis cam
lift = @(x,angdisp,starting_point) (1-cos(2*pi*(x-starting_point)/angdisp) ) /2*Lh;

% Define Cam Profile
z_points = zeros(1,length(range_full));
z_points(span.s1)= 0;
z_points(span.sa)= lift(theta(span.sa),angdisp.sa,theta(span.sa(1)));
z_points(span.s2)= 0;
z_points(span.sb)= lift(theta(span.sb),angdisp.sb,theta(span.sb(1)));
z_points(span.s3)= 0;
z_points(span.sc)= lift(theta(span.sc),angdisp.sc,theta(span.sc(1)));

% Points are already smooth--used Harmonic Function. No need to interpolate
% or fit data.

z_disp_smooth = z_points';  



% Generate Final Cam Profiles.

x_cam_profile = x_disp_smoothfit + Brx;
y_cam_profile = y_disp_smoothfit + Bry;
z_cam_profile = z_disp_smooth    + Brz;



% Draw Cams

figure(2)
subplot(1,3,1)
    polarplot(theta(span.s1),x_cam_profile(span.s1),...
              theta(span.sa),x_cam_profile(span.sa),...
              theta(span.s2),x_cam_profile(span.s2),...
              theta(span.sb),x_cam_profile(span.sb),...
              theta(span.s3),x_cam_profile(span.s3),...
              theta(span.sc),x_cam_profile(span.sc),'LineWidth',2)
    title('X Cam Profile')
subplot(1,3,2)
    polarplot(theta(span.s1),y_cam_profile(span.s1),...
              theta(span.sa),y_cam_profile(span.sa),...
              theta(span.s2),y_cam_profile(span.s2),...
              theta(span.sb),y_cam_profile(span.sb),...
              theta(span.s3),y_cam_profile(span.s3),...
              theta(span.sc),y_cam_profile(span.sc),'LineWidth',2)
    title('Y Cam Profile')
subplot(1,3,3)
    polarplot(theta(span.s1),z_cam_profile(span.s1),...
              theta(span.sa),z_cam_profile(span.sa),...
              theta(span.s2),z_cam_profile(span.s2),...
              theta(span.sb),z_cam_profile(span.sb),...
              theta(span.s3),z_cam_profile(span.s3),...
              theta(span.sc),z_cam_profile(span.sc),'LineWidth',2)    
    title('Z Cam Profile')
    legend('1st Initial', 'Lift', '2nd Initial', 'Lift', '3rd Initial', 'Lift','Location','Best')

    
    
    
%% Analyze Cams

% Kinematic Coefficients / s,v,a,j 

 x_s = x_disp_smoothfit;
 y_s = y_disp_smoothfit;
 z_s = z_disp_smooth;

x_k1 = diff(x_cam_profile)./diff(theta);
y_k1 = diff(y_cam_profile)./diff(theta);
z_k1 = diff(z_cam_profile)./diff(theta);

x_k2 = diff(x_k1)./diff(theta(1:length(theta)-1));
y_k2 = diff(y_k1)./diff(theta(1:length(theta)-1));
z_k2 = diff(z_k1)./diff(theta(1:length(theta)-1));

x_j = diff(x_k2)./diff(theta(1:length(theta)-2));
y_j = diff(y_k2)./diff(theta(1:length(theta)-2));
z_j = diff(z_k2)./diff(theta(1:length(theta)-2));

% Radius of Curvature and Pressure Angle

x_rho = ((x_s(1:end-2)+Brx).^2 + x_k1(1:end-1).^2).^(3/2) ./ ...
        ((x_s(1:end-2)+Brx).^2 + 2.*x_k1(1:end-1).^2 - x_k2.*(x_s(1:end-2)+Brx));
x_pa = atan(x_k1./(x_s(1:end-1)+Brx))*180/pi;

y_rho = ((y_s(1:end-2)+Bry).^2 + y_k1(1:end-1).^2).^(3/2) ./ ...
        ((y_s(1:end-2)+Bry).^2 + 2.*y_k1(1:end-1).^2 - y_k2.*(y_s(1:end-2)+Bry));
y_pa = atan(y_k1./(y_s(1:end-1)+Bry))*180/pi;

z_rho = ((z_s(1:end-2)+Brz).^2 + z_k1(1:end-1).^2).^(3/2) ./ ...
        ((z_s(1:end-2)+Brz).^2 + 2.*z_k1(1:end-1).^2 - z_k2.*(z_s(1:end-2)+Brz));
z_pa = atan(z_k1./(z_s(1:end-1)+Brz))*180/pi;



% Plot Data

figure(3)
set(gcf, 'Position',  [100, 00, 1500, 700])
subplot(2,4,1)
    plot(thetad,x_cam_profile-Brx,'LineWidth',1)
    title('X Cam Displacement vs. Theta')
    ylabel('Displacement (inches)')
    xlabel('\theta (deg)')
    axis([0 360 0 2])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,2)
    plot(thetad(1:length(thetad)-1),x_k1,'LineWidth',1)
    title('X Cam: 1st Order KC')
    xlabel('\omega (deg/s)')
    axis([0 360 -4 2])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,3)
    plot(thetad(1:length(thetad)-2),x_k2,'LineWidth',1)
    title('X Cam: 2nd Order KC')
    xlabel('\alpha (deg/s^2)')
    axis([0 360 -20 20])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,4)
    plot(thetad(1:length(thetad)-3),x_j,'LineWidth',1)
    title('X Cam: Jerk')
    xlabel('d\alpha (deg/s^3)')
%     axis([0 360 -20 20])
    xticks([0 60 120 180 240 300 360])
    xlim([0 360])
    grid
subplot(2,4,6)
    plot(thetad(1:end-2),x_rho,'LineWidth',1)
    title('X Cam: Radius of Curvature')
    ylabel('R.O.C. (inches)')
    xlabel('\theta (deg)')
    axis([0 360 -20 20])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,7)
    plot(thetad(1:end-1),x_pa,'LineWidth',1)
    title('X Cam: Pressure Angle')
    ylabel('Pressure Angle (deg)')
    xlabel('\theta (deg)')
    xlim([0 360])
%     axis([0 360 -3 3])
    xticks([0 60 120 180 240 300 360])
    grid
    
figure(4)
set(gcf, 'Position',  [100, 00, 1500, 700])
subplot(2,4,1)
    plot(thetad,y_cam_profile-Bry,'LineWidth',1)
    title('Y Cam Displacement vs. Theta')
    ylabel('Displacement (inches)')
    xlabel('\theta (deg)')
    axis([0 360 0 1])
    grid
subplot(2,4,2)
    plot(thetad(1:length(thetad)-1),y_k1,'LineWidth',1)
    title('Y Cam: 1st Order KC')
    xlabel('\omega (deg/s)')
    axis([0 360 -4 2])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,3)
    plot(thetad(1:length(thetad)-2),y_k2,'LineWidth',1)
    title('Y Cam: 2nd Order KC')
    xlabel('\alpha (deg/s^2)')
    axis([0 360 -20 20])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,4)
    plot(thetad(1:length(thetad)-3),y_j,'LineWidth',1)
    title('Y Cam: Jerk')
    xlabel('d\alpha (deg/s^3)')
%     axis([0 360 -20 20])
    xlim([0 360])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,6)
    plot(thetad(1:end-2),y_rho,'LineWidth',1)
    title('Y Cam: Radius of Curvature')
    ylabel('R.O.C. (inches)')
    xlabel('\theta (deg)')
    axis([0 360 -20 20])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,7)
    plot(thetad(1:end-1),y_pa,'LineWidth',1)
    title('Y Cam: Pressure Angle')
    ylabel('Pressure Angle (deg)')
    xlabel('\theta (deg)')
%     axis([0 360 -3 3])
    xlim([0 360])
    xticks([0 60 120 180 240 300 360])
    grid
    
    
    
figure(5)
set(gcf, 'Position',  [100, 00, 1500, 700])
subplot(2,4,1)
    plot(thetad,z_cam_profile-Brz,'LineWidth',1)
    title('Z Cam Displacement vs. Theta')
    ylabel('Displacement (inches)')
    xlabel('\theta (deg)')
    axis([0 360 -.4 .4])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,2)
    plot(thetad(1:length(thetad)-1),z_k1,'LineWidth',1)
    title('Z Cam: 1st Order KC')
    xlabel('\omega (deg/s)')
    axis([0 360 -3 3])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,3)
    plot(thetad(1:length(thetad)-2),z_k2,'LineWidth',1)
    title('Z Cam: 2nd Order KC')
    xlabel('\alpha (deg/s^2)')
    axis([0 360 -20 20])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,4)
    plot(thetad(1:length(thetad)-3),z_j,'LineWidth',1)
    title('Z Cam: Jerk')
    xlabel('d\alpha (deg/s^3)')
    axis([0 360 -20 20])
    xticks([0 60 120 180 240 300 360])
%     xlim([0 360])
    grid
subplot(2,4,6)
    plot(thetad(1:end-2),z_rho,'LineWidth',1)
    title('Z Cam: Radius of Curvature')
    ylabel('R.O.C. (inches)')
    xlabel('\theta (deg)')
    axis([0 360 -10 10])
    xticks([0 60 120 180 240 300 360])
    grid
subplot(2,4,7)
    plot(thetad(1:end-1),z_pa,'LineWidth',1)
    title('Z Cam: Pressure Angle')
    ylabel('Pressure Angle (deg)')
    xlabel('\theta (deg)')
%     axis([0 360 -3 3])
    xlim([0 360])
    xticks([0 60 120 180 240 300 360])
    grid    
    
    
    
% X-Y (Cartesian) Cam Profile for Solid Model Generation (Hector)
% Note: 3rd columns in the below matricies are due to SolidWorks coordinate
% file loading requirements. (Need X Y Z)

x_cam_cart(:,1) = x_cam_profile.*cos(theta);
x_cam_cart(:,2) = x_cam_profile.*sin(theta);
x_cam_cart(:,3) = .75;

y_cam_cart(:,1) = y_cam_profile.*cos(theta);
y_cam_cart(:,2) = y_cam_profile.*sin(theta);
y_cam_cart(:,3) = 1.5;

z_cam_cart(:,1) = z_cam_profile.*cos(theta);
z_cam_cart(:,2) = z_cam_profile.*sin(theta);
z_cam_cart(:,3) = 0;

% Export Data to .txt files loadable by CAD software   (Hector)

save('x_cam_cart.txt', 'x_cam_cart','-ascii', '-tabs');
save('y_cam_cart.txt', 'y_cam_cart','-ascii', '-tabs');
save('z_cam_cart.txt', 'z_cam_cart','-ascii', '-tabs');