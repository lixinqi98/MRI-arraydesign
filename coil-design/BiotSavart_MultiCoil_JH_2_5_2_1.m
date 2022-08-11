close all
clear all
echo off

% This code generates the B0 maps for a multi-coil array assuming that:
% - each coil is a circular single-loop coil
% - each ring is rotated by ang/2 along z
% - the distance between adjacent coils is the same in each direction
tic
rr = 0.24/2; % radius of each ring (m)
FOV = 0.21; % FOV of the B0 maps to shim (m)
nx0 = 70;
ny0 = nx0;
res = FOV/nx0; % resolution of the B0 maps to shim (m)
nz0 = 40; % number of slices to shim
th = 0.003;
height = 2*rr;
I = 1; % reference current in each coil (A)

range = max(2*rr,FOV);

[x,y,z] = meshgrid(-(range-res)/2:res:(range-res)/2,...
    -(range-res)/2:res:(range-res)/2,-(height-th)/2:th:(height-th)/2);

N7 = 8;
N8 = 16;
N = 32+N7+N8;
Bz = zeros([size(x,1) size(x,2) size(x,3) N]);
coils = zeros([size(x,1) size(x,2) size(x,3)]);
cordf = zeros(N,3,360);

rr1 = 0.20/2;
%% Ellipse

uv = [0 1 0];

N1 = 2;

a = 0.06;
b = 0.04;

x1 = a*cos((0:359)*pi/180);
x2 = a*cos((1:360)*pi/180);

y1 = rr*ones(1,size(x1,2))*0.9;
y2 = y1;

z1 = b*sin((0:359)*pi/180);
z2 = b*sin((1:360)*pi/180);

coordi1 = [x1; y1; z1];
coordi2 = [x2; y2; z2];


for jr = 1:N1
    
    c = [0 0 1]'.*rr*0;
    zoffset0 = repmat(c,1,size(coordi1,2));
    
    if jr==1
        angle00z = -22/180*pi;
    else
        angle00z = 22/180*pi;
    end
    
    
    rot3 = [cos(angle00z) -sin(angle00z) 0; sin(angle00z) cos(angle00z) 0; 0 0 1];
    coordi1_2 = rot3*coordi1-zoffset0.*1; coordi2_2 = rot3*coordi2-zoffset0.*1;
    
    
    
    x1f = reshape(coordi1_2(1,:),[1 size(coordi1_2,2)]);
    y1f = reshape(coordi1_2(2,:),[1 size(coordi1_2,2)]);
    z1f = reshape(coordi1_2(3,:),[1 size(coordi1_2,2)]);
    
    x2f = reshape(coordi2_2(1,:),[1 size(coordi1_2,2)]);
    y2f = reshape(coordi2_2(2,:),[1 size(coordi1_2,2)]);
    z2f = reshape(coordi2_2(3,:),[1 size(coordi1_2,2)]);
    
    dlx = x2f-x1f;
    dly = y2f-y1f;
    dlz = z2f-z1f;
    %
%     for k = 1:360
%         rx = x - (x1f(k)+x2f(k))/2;
%         ry = y - (y1f(k)+y2f(k))/2;
%         rz = z - (z1f(k)+z2f(k))/2;
%         r3 = sqrt(rx.^2 + ry.^2 + rz.^2).^3;
%         Bz(:,:,:,jr) = Bz(:,:,:,jr) + ...
%             1e-7*I*(dlx(k)*ry - dly(k)*rx) ./ r3;
%     end
    
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])
    
    cordf(jr,1,:) = x1f;
    cordf(jr,2,:) = y1f;
    cordf(jr,3,:) = z1f;
    
    cordf1(jr,1,:) = x2f;
    cordf1(jr,2,:) = y2f;
    cordf1(jr,3,:) = z2f;
    
end
%%

% radius of each coil (m)
rc = 0.0475;

x1 = rc*cos((0:359)*pi/180);
x2 = rc*cos((1:360)*pi/180);

z1 = rc*sin((0:359)*pi/180);
z2 = rc*sin((1:360)*pi/180);




%% 1 raws

N2  = 5;

rr01 = (rr+rr1)/2;
% angle1y00 = 30;
% angle1y01 = 150;
angle1y0  = [0 linspace(40,180-40,N2-2) 0]*pi/180;
anglezz = 55;
angle1z0 = [1.15 0.98 0.9 0.98 -1.15].*anglezz*pi/180;
anglexx = -5;
angle1x0 = [0 1 2 1 0].*anglexx*pi/180;

for jr = N1+1:N1+N2
    
    y1 = rr01*ones(1,size(x1,2));
    y2 = y1;
    coordi1 = [x1; y1; z1];
    coordi2 = [x2; y2; z2];
    
    angle1z = angle1z0(jr-N1);
    angle1y = angle1y0(jr-N1);
    angle1x = angle1x0(jr-N1);
    
    %     rot3 = [cos(angle1z) -sin(angle1z) 0; sin(angle1z) cos(angle1z) 0; 0 0 1];
    %     rot2 = [cos(angle1y) 0 sin(angle1y); 0 1 0;  -sin(angle1y) 0 cos(angle1y)];
    %     rot1 = [1 0 0;0 cos(angle1x) -sin(angle1x); 0 sin(angle1x) cos(angle1x)];
    
    c = [0 0 1]'.*0.005;
    zoffset0 = repmat(c,1,size(coordi1,2));
    
    
    if jr== N1+1 | jr== N1+N2
        rot3 = [cos(angle1z) -sin(angle1z) 0; sin(angle1z) cos(angle1z) 0; 0 0 1];
        coordi1_2 = rot3*(coordi1-zoffset0); coordi2_2 = rot3*(coordi2-zoffset0);
    else
        
        rot3 = [cos(angle1z) -sin(angle1z) 0; sin(angle1z) cos(angle1z) 0; 0 0 1];
        rot2 = [cos(angle1y) 0 sin(angle1y); 0 1 0;  -sin(angle1y) 0 cos(angle1y)];
        rot1 = [1 0 0;0 cos(angle1x) -sin(angle1x); 0 sin(angle1x) cos(angle1x)];
        coordi1_2 = rot1*rot2*rot3*coordi1; coordi2_2 = rot1*rot2*rot3*coordi2;
    end
    coordi1_1 =  coordi1_2;
    
    
    %     inv_rot2 = [cos(angle1y) 0 -sin(angle1y); 0 1 0;  sin(angle1y) 0 cos(angle1y)];
    
    
    
    x1f = reshape(coordi1_2(1,:),[1 size(coordi1_1,2)]);
    y1f = reshape(coordi1_2(2,:),[1 size(coordi1_1,2)]);
    z1f = reshape(coordi1_2(3,:),[1 size(coordi1_1,2)]);
    
    x2f = reshape(coordi2_2(1,:),[1 size(coordi1_1,2)]);
    y2f = reshape(coordi2_2(2,:),[1 size(coordi1_1,2)]);
    z2f = reshape(coordi2_2(3,:),[1 size(coordi1_1,2)]);
    
    dlx = x2f-x1f;
    dly = y2f-y1f;
    dlz = z2f-z1f;
    
%     for k = 1:360
%         rx = x - (x1f(k)+x2f(k))/2;
%         ry = y - (y1f(k)+y2f(k))/2;
%         rz = z - (z1f(k)+z2f(k))/2;
%         r3 = sqrt(rx.^2 + ry.^2 + rz.^2).^3;
%         Bz(:,:,:,jr) = Bz(:,:,:,jr) + ...
%             1e-7*I*(dlx(k)*ry - dly(k)*rx) ./ r3;
%     end
    
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])
    
    cordf(jr,1,:) = x1f;
    cordf(jr,2,:) = y1f;
    cordf(jr,3,:) = z1f;
    
    cordf1(jr,1,:) = x2f;
    cordf1(jr,2,:) = y2f;
    cordf1(jr,3,:) = z2f;
    
    jr
    
end

%% 3 rows

N3  = 8;

rr02 = rr1+(rr-rr1)*0.1;
% angle1y00 = 30;
% angle1y01 = 150;
angle2y0  = [0 linspace(0,180,N3-2) 0].*pi/180;
anglezz = 90;
angle2z0 = [1.05 1.07 0.97 1 1 0.97 -1.07 -1.05].*anglezz*pi/180;
anglexx = -5;
angle2x0 = [0 0.5 2 4 4 2 0.5 0].*anglexx*pi/180;

% anglezz = 55;
% angle1z0 = [1 0.95 0.9 0.95 -1].*anglezz*pi/180;
% anglexx = -5;
% angle1x0 = [0 0.8 2 0.8 0].*anglexx*pi/180;


for jr = N1+N2+1:N1+N2+N3
    
    y1 = rr01*ones(1,size(x1,2));
    y2 = y1;
    coordi1 = [x1; y1; z1];
    coordi2 = [x2; y2; z2];
    
    angle2z = angle2z0(jr-N1-N2);
    angle2y = angle2y0(jr-N1-N2);
    angle2x = angle2x0(jr-N1-N2);
    
    rot3 = [cos(angle2z) -sin(angle2z) 0; sin(angle2z) cos(angle2z) 0; 0 0 1];
    rot2 = [cos(angle2y) 0 sin(angle2y); 0 1 0;  -sin(angle2y) 0 cos(angle2y)];
    rot1 = [1 0 0;0 cos(angle2x) -sin(angle2x); 0 sin(angle2x) cos(angle2x)];
    
    c = [0 0 1]'.*0.06;
    zoffset0 = repmat(c,1,size(coordi1,2));
    
    
    if jr== N1+N2+1 | jr== N1+N2+N3
        %     rot3 = [cos(angle2z) -sin(angle2z) 0; sin(angle2z) cos(angle2z) 0; 0 0 1];
        coordi1_2 = rot3*coordi1-zoffset0; coordi2_2 = rot3*coordi2-zoffset0;
    elseif jr== N1+N2+2 | jr== N1+N2+N3-1
        %     rot3 = [cos(angle2z) -sin(angle2z) 0; sin(angle2z) cos(angle2z) 0; 0 0 1];
        coordi1_2 = rot3*coordi1-zoffset0.*0; coordi2_2 = rot3*coordi2-zoffset0.*0;
        
    else
        coordi1_2 = rot1*rot2*rot3*coordi1; coordi2_2 = rot1*rot2*rot3*coordi2;
    end
    coordi1_1 =  coordi1_2;
    
    
    %     inv_rot2 = [cos(angle1y) 0 -sin(angle1y); 0 1 0;  sin(angle1y) 0 cos(angle1y)];
    
    
    
    x1f = reshape(coordi1_2(1,:),[1 size(coordi1_1,2)]);
    y1f = reshape(coordi1_2(2,:),[1 size(coordi1_1,2)]);
    z1f = reshape(coordi1_2(3,:),[1 size(coordi1_1,2)]);
    
    x2f = reshape(coordi2_2(1,:),[1 size(coordi1_1,2)]);
    y2f = reshape(coordi2_2(2,:),[1 size(coordi1_1,2)]);
    z2f = reshape(coordi2_2(3,:),[1 size(coordi1_1,2)]);
    
    dlx = x2f-x1f;
    dly = y2f-y1f;
    dlz = z2f-z1f;
    
%     for k = 1:360
%         rx = x - (x1f(k)+x2f(k))/2;
%         ry = y - (y1f(k)+y2f(k))/2;
%         rz = z - (z1f(k)+z2f(k))/2;
%         r3 = sqrt(rx.^2 + ry.^2 + rz.^2).^3;
%         Bz(:,:,:,jr) = Bz(:,:,:,jr) + ...
%             1e-7*I*(dlx(k)*ry - dly(k)*rx) ./ r3;
%     end
    
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])
    
    cordf(jr,1,:) = x1f;
    cordf(jr,2,:) = y1f;
    cordf(jr,3,:) = z1f;
    
    cordf1(jr,1,:) = x2f;
    cordf1(jr,2,:) = y2f;
    cordf1(jr,3,:) = z2f;
    
    jr
    
end

%% 3 raw (Back 180)

N4  = 5;

rr03 = rr;
% angle1y00 = 30;
% angle1y01 = 150;
angle3y0  = [0 0 0 0 0].*pi/180;
anglezz = 180;
angle3z0 = [1 1 1 1 1].*anglezz*pi/180;
% anglexx = ;
angle3x0 = [0 linspace(-10,85,N4-1)].*pi/180;
% angle3x0 = [0 0 30 60 90].*pi/180;

for jr = N1+N2+N3+1:N1+N2+N3+N4
    
    x1 = rc*cos((0:359)*pi/180);
    x2 = rc*cos((1:360)*pi/180);
    
    z1 = rc*sin((0:359)*pi/180);
    z2 = rc*sin((1:360)*pi/180);
    
    y1 = rr03*ones(1,size(x1,2));
    y2 = y1;
    coordi1 = [x1; y1; z1];
    coordi2 = [x2; y2; z2];
    
    angle3z = angle3z0(jr-N1-N2-N3);
    angle3x = angle3x0(jr-N1-N2-N3);
    angle3y = angle3y0(jr-N1-N2-N3);
    
    rot3 = [cos(angle3z) -sin(angle3z) 0; sin(angle3z) cos(angle3z) 0; 0 0 1];
    rot2 = [cos(angle3y) 0 sin(angle3y); 0 1 0;  -sin(angle3y) 0 cos(angle3y)];
    rot1 = [1 0 0; 0 cos(angle3x) -sin(angle3x); 0 sin(angle3x) cos(angle3x)];
    
    
    cz = [0 0 1]'.*0.07;
    cy = [0 1 0]'.*0.01;
    
    zoffset0 = repmat(cz,1,size(coordi1,2));
    yoffset0 = repmat(cy,1,size(coordi1,2));
    
    
    if jr== N1+N2+N3+1
        coordi1_2 = rot3*coordi1-zoffset0.*1+yoffset0 ; coordi2_2 = rot3*coordi2-zoffset0.*1+yoffset0;
        coordi1_1 = coordi1_2;
        
        %     elseif jr== N1+N2+N3+2
        %         coordi1_2 = rot3*coordi1-zoffset0.*1 ; coordi2_2 = rot3*coordi2-zoffset0.*1 ;
        %          coordi1_1 = coordi1_2;
        
        
    else
        
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot3*coordi1_1; coordi2_2 = rot3*coordi2_1;
        
    end
    
    x1f = reshape(coordi1_2(1,:),[1 size(coordi1_1,2)]);
    y1f = reshape(coordi1_2(2,:),[1 size(coordi1_1,2)]);
    z1f = reshape(coordi1_2(3,:),[1 size(coordi1_1,2)]);
    
    x2f = reshape(coordi2_2(1,:),[1 size(coordi1_1,2)]);
    y2f = reshape(coordi2_2(2,:),[1 size(coordi1_1,2)]);
    z2f = reshape(coordi2_2(3,:),[1 size(coordi1_1,2)]);
    
    dlx = x2f-x1f;
    dly = y2f-y1f;
    dlz = z2f-z1f;
    
%     for k = 1:360
%         rx = x - (x1f(k)+x2f(k))/2;
%         ry = y - (y1f(k)+y2f(k))/2;
%         rz = z - (z1f(k)+z2f(k))/2;
%         r3 = sqrt(rx.^2 + ry.^2 + rz.^2).^3;
%         Bz(:,:,:,jr) = Bz(:,:,:,jr) + ...
%             1e-7*I*(dlx(k)*ry - dly(k)*rx) ./ r3;
%     end
    %
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])
    
    cordf(jr,1,:) = x1f;
    cordf(jr,2,:) = y1f;
    cordf(jr,3,:) = z1f;
    
    cordf1(jr,1,:) = x2f;
    cordf1(jr,2,:) = y2f;
    cordf1(jr,3,:) = z2f;
    
    jr
    
    clear rot1 rot2 rot3
    
end

%% 4 raw

N5 = 8;

rr04 = rr1+(rr-rr1)*1;
% angle1y00 = 30;
% angle1y01 = 150;
% angle4y0  = [0 30 60 -60 -30 -0].*pi/180;
anglezz = 145;
angle4z0 = [0.8 1.01 1.03 1.05 -0.8 -1.01 -1.03 -1.05].*anglezz*pi/180;
angle4x0 = [55 30 5 0 55 30 5 0].*pi/180;

anglexx = 0;
angle4x01 = [anglexx anglexx anglexx anglexx anglexx anglexx anglexx anglexx].*pi/180;



for jr = N1+N2+N3+N4+1:N1+N2+N3+N4+N5
    
    y1 = rr04*ones(1,size(x1,2));
    y2 = y1;
    coordi1 = [x1; y1; z1];
    coordi2 = [x2; y2; z2];
    
    angle4z = angle4z0(jr-N1-N2-N3-N4);
    %     angle4y = angle4y0(jr-N1-N2-N3-N4);
    angle4x = angle4x0(jr-N1-N2-N3-N4);
    angle4x_1 = angle4x01(jr-N1-N2-N3-N4);
    
    rot3 = [cos(angle4z) -sin(angle4z) 0; sin(angle4z) cos(angle4z) 0; 0 0 1];
    %     rot2 = [cos(angle4y) 0 sin(angle4y); 0 1 0;  -sin(angle4y) 0 cos(angle4y)];
    rot1 = [1 0 0;0 cos(angle4x) -sin(angle4x); 0 sin(angle4x) cos(angle4x)];
    rot1_1 = [1 0 0;0 cos(angle4x_1) -sin(angle4x_1); 0 sin(angle4x_1) cos(angle4x_1)];
    
    cz = [0 0 1]'.*0.05;
    cy = [0 1 0]'.*0.01;
    zoffset0 = repmat(cz,1,size(coordi1,2));
    yoffset0 = repmat(cy,1,size(coordi1,2));
    
    
    if jr == N1+N2+N3+N4+4 | jr== N1+N2+N3+N4+N5
        
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot1_1*rot3*coordi1_1-zoffset0 ; coordi2_2 = rot1_1*rot3*coordi2_1-zoffset0 ;
    else
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot1_1*rot3*coordi1_1 ; coordi2_2 = rot1_1*rot3*coordi2_1 ;
    end
    
    
    
    x1f = reshape(coordi1_2(1,:),[1 size(coordi1_1,2)]);
    y1f = reshape(coordi1_2(2,:),[1 size(coordi1_1,2)]);
    z1f = reshape(coordi1_2(3,:),[1 size(coordi1_1,2)]);
    
    x2f = reshape(coordi2_2(1,:),[1 size(coordi1_1,2)]);
    y2f = reshape(coordi2_2(2,:),[1 size(coordi1_1,2)]);
    z2f = reshape(coordi2_2(3,:),[1 size(coordi1_1,2)]);
    
    dlx = x2f-x1f;
    dly = y2f-y1f;
    dlz = z2f-z1f;
    
%     for k = 1:360
%         rx = x - (x1f(k)+x2f(k))/2;
%         ry = y - (y1f(k)+y2f(k))/2;
%         rz = z - (z1f(k)+z2f(k))/2;
%         r3 = sqrt(rx.^2 + ry.^2 + rz.^2).^3;
%         Bz(:,:,:,jr) = Bz(:,:,:,jr) + ...
%             1e-7*I*(dlx(k)*ry - dly(k)*rx) ./ r3;
%     end
    %
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])
    
    cordf(jr,1,:) = x1f;
    cordf(jr,2,:) = y1f;
    cordf(jr,3,:) = z1f;
    
    cordf1(jr,1,:) = x2f;
    cordf1(jr,2,:) = y2f;
    cordf1(jr,3,:) = z2f;
    
    jr
    
end

%% 5 raw

N6 = 4;

rr05 = rr1+(rr-rr1)*0.7;
% angle1y00 = 30;
% angle1y01 = 150;
% angle5y0  = [0 40 70 -70 -40 -0].*pi/180;
anglezz = 135;
angle5z0 = [0.95 0.9 -0.95 -0.9].*anglezz*pi/180;
anglexx = 5;
angle5x0 = [-2 35 -2 35].*pi/180;



for jr = N1+N2+N3+N4+N5+1:N1+N2+N3+N4+N5+N6
    
    y1 = rr05*ones(1,size(x1,2));
    y2 = y1;
    coordi1 = [x1; y1; z1];
    coordi2 = [x2; y2; z2];
    
    angle5z = angle5z0(jr-N1-N2-N3-N4-N5);
    %     angle5y = angle5y0(jr-N1-N2-N3-N4-N5);
    angle5x = angle5x0(jr-N1-N2-N3-N4-N5);
    
    rot3 = [cos(angle5z) -sin(angle5z) 0; sin(angle5z) cos(angle5z) 0; 0 0 1];
    %     rot2 = [cos(angle5y) 0 sin(angle5y); 0 1 0;  -sin(angle5y) 0 cos(angle5y)];
    rot1 = [1 0 0;0 cos(angle5x) -sin(angle5x); 0 sin(angle5x) cos(angle5x)];
    
    
    cz = [0 0 1]'.*0.02;
    cy = [0 1 0]'.*0.01;
    
    zoffset0 = repmat(cz,1,size(coordi1,2));
    yoffset0 = repmat(cy,1,size(coordi1,2));
    
    
    if jr== N1+N2+N3+N4+N5+1 | N1+N2+N3+N4+N5+3
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot3*coordi1_1-zoffset0 ; coordi2_2 = rot3*coordi2_1-zoffset0 ;
    else
        
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot3*coordi1_1; coordi2_2 = rot3*coordi2_1;
    end
    
    x1f = reshape(coordi1_2(1,:),[1 size(coordi1_1,2)]);
    y1f = reshape(coordi1_2(2,:),[1 size(coordi1_1,2)]);
    z1f = reshape(coordi1_2(3,:),[1 size(coordi1_1,2)]);
    
    x2f = reshape(coordi2_2(1,:),[1 size(coordi1_1,2)]);
    y2f = reshape(coordi2_2(2,:),[1 size(coordi1_1,2)]);
    z2f = reshape(coordi2_2(3,:),[1 size(coordi1_1,2)]);
    
    dlx = x2f-x1f;
    dly = y2f-y1f;
    dlz = z2f-z1f;
    
%     for k = 1:360
%         rx = x - (x1f(k)+x2f(k))/2;
%         ry = y - (y1f(k)+y2f(k))/2;
%         rz = z - (z1f(k)+z2f(k))/2;
%         r3 = sqrt(rx.^2 + ry.^2 + rz.^2).^3;
%         Bz(:,:,:,jr) = Bz(:,:,:,jr) + ...
%             1e-7*I*(dlx(k)*ry - dly(k)*rx) ./ r3;
%     end
    
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])
    
    cordf(jr,1,:) = x1f;
    cordf(jr,2,:) = y1f;
    cordf(jr,3,:) = z1f;
    
    cordf1(jr,1,:) = x2f;
    cordf1(jr,2,:) = y2f;
    cordf1(jr,3,:) = z2f;
    
    jr
    
end

%% Modification

%% Mona: change from 2*4 to 3*4

N7 = 12;

% radius of each coil (m)
rc0 = 0.025;


x1 = rc0*cos((0:359)*pi/180);
x2 = rc0*cos((1:360)*pi/180);

y1 = rr*ones(1,size(x1,2));
y2 = y1;

z1 = rc0*sin((0:359)*pi/180);
z2 = rc0*sin((1:360)*pi/180);

coordi1 = [x1; y1; z1];
coordi2 = [x2; y2; z2];

for jr = N1+N2+N3+N4+N5+N6+1:N1+N2+N3+N4+N5+N6+N7
    
    rc  = rc0;
    
    angle01 = atan(rc/rr)*2;
    angle02 = atan(rc/rr)*4;
    
    % for each ring, generate the B0 map for the first coil
    
    zz = 0.04;
    c = [0 0 1]'.*zz;
    zoffset1 = repmat(c,1,size(coordi1,2));
    
    if jr == N1+N2+N3+N4+N5+N6+1
        angle1 = - angle01*0.4;
        angle2 =  angle02*3/6;
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        rot2 = [1 0 0; 0 cos(angle2) -sin(angle2); 0 sin(angle2) cos(angle2)];
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot2*coordi1_1-zoffset1.*0; coordi2_2 = rot2*coordi2_1-zoffset1.*0;
        
    elseif jr == N1+N2+N3+N4+N5+N6+2
        angle1 = angle01*0.4;
        angle2 =  angle02*3/6;
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        rot2 = [1 0 0; 0 cos(angle2) -sin(angle2); 0 sin(angle2) cos(angle2)];
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot2*coordi1_1-zoffset1.*0; coordi2_2 = rot2*coordi2_1-zoffset1.*0;
        
        zz = 0.04;
        c = [0 0 1]'.*zz;
        zoffset1 = repmat(c,1,size(coordi1,2));
        
    elseif jr == N1+N2+N3+N4+N5+N6+3
        angle1 = -angle01*0.4;
        angle2 =  angle02*5/6;
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        rot2 = [1 0 0; 0 cos(angle2) -sin(angle2); 0 sin(angle2) cos(angle2)];
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot2*coordi1_1-zoffset1*0; coordi2_2 = rot2*coordi2_1-zoffset1*0;
        
        
    elseif jr == N1+N2+N3+N4+N5+N6+4
        angle1 = angle01*0.4;
        angle2 = angle02*5/6;
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        rot2 = [1 0 0; 0 cos(angle2) -sin(angle2); 0 sin(angle2) cos(angle2)];
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot2*coordi1_1-zoffset1*0; coordi2_2 = rot2*coordi2_1-zoffset1*0;
        
    elseif jr == N1+N2+N3+N4+N5+N6+5
        angle1 = -angle01*0.4;
        angle2 =  angle02*7/6;
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        rot2 = [1 0 0; 0 cos(angle2) -sin(angle2); 0 sin(angle2) cos(angle2)];
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot2*coordi1_1-zoffset1*0; coordi2_2 = rot2*coordi2_1-zoffset1*0;
        
        
    elseif jr == N1+N2+N3+N4+N5+N6+6
        angle1 = angle01*0.4;
        angle2 = angle02*7/6;
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        rot2 = [1 0 0; 0 cos(angle2) -sin(angle2); 0 sin(angle2) cos(angle2)];
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot2*coordi1_1-zoffset1*0; coordi2_2 = rot2*coordi2_1-zoffset1*0;
        
    elseif jr == N1+N2+N3+N4+N5+N6+7
        angle1 = -angle01*0.4;
        angle2 =  angle02*9/6;
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        rot2 = [1 0 0; 0 cos(angle2) -sin(angle2); 0 sin(angle2) cos(angle2)];
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot2*coordi1_1-zoffset1*0; coordi2_2 = rot2*coordi2_1-zoffset1*0;
        
        
    elseif jr == N1+N2+N3+N4+N5+N6+8
        angle1 = angle01*0.4;
        angle2 = angle02*9/6;
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        rot2 = [1 0 0; 0 cos(angle2) -sin(angle2); 0 sin(angle2) cos(angle2)];
        
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = rot2*coordi1_1-zoffset1*0; coordi2_2 = rot2*coordi2_1-zoffset1*0;
        
    end
    
    
    x1f = reshape(coordi1_2(1,:),[1 size(coordi1_1,2)]);
    y1f = reshape(coordi1_2(2,:),[1 size(coordi1_1,2)]);
    z1f = reshape(coordi1_2(3,:),[1 size(coordi1_1,2)]);
    
    x2f = reshape(coordi2_2(1,:),[1 size(coordi1_1,2)]);
    y2f = reshape(coordi2_2(2,:),[1 size(coordi1_1,2)]);
    z2f = reshape(coordi2_2(3,:),[1 size(coordi1_1,2)]);
    
    dlx = x2f-x1f;
    dly = y2f-y1f;
    dlz = z2f-z1f;
    %
%     for k = 1:360
%         rx = x - (x1f(k)+x2f(k))/2;
%         ry = y - (y1f(k)+y2f(k))/2;
%         rz = z - (z1f(k)+z2f(k))/2;
%         r3 = sqrt(rx.^2 + ry.^2 + rz.^2).^3;
%         Bz(:,:,:,jr) = Bz(:,:,:,jr) + ...
%             1e-7*I*(dlx(k)*ry - dly(k)*rx) ./ r3;
%     end
%     
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])
    
    cordf(jr,1,:) = x1f;
    cordf(jr,2,:) = y1f;
    cordf(jr,3,:) = z1f;
    
    cordf1(jr,1,:) = x2f;
    cordf1(jr,2,:) = y2f;
    cordf1(jr,3,:) = z2f;
    
    jr
    
end


%% Side 4cm, 16

N8 = 16;

rcc = 0.02;
x1 = rcc*cos((0:359)*pi/180);
x2 = rcc*cos((1:360)*pi/180);

rrs = 0.11;

y1 = rrs*ones(1,size(x1,2));
y2 = y1;

z1 = rcc*sin((0:359)*pi/180);
z2 = rcc*sin((1:360)*pi/180);

coordi1 = [x1; y1; z1];
coordi2 = [x2; y2; z2];

for jr = N1+N2+N3+N4+N5+N6+N7+1:N1+N2+N3+N4+N5+N6+N7+N8
    
    angle0001 = 82;
    angle0002 = 98;
    
    angle01 = [angle0001 angle0001 angle0001 angle0001 angle0002 angle0002 angle0002 angle0002...
        -angle0001 -angle0001 -angle0001 -angle0001 -angle0002 -angle0002 -angle0002 -angle0002]*pi/180;
    angle02 = 20/180*pi;
    
    c = [0 0 1]'.*0.05;
    zoffset0 = repmat(c,1,size(coordi1,2));
    
    nnn = N1+N2+N3+N4+N5+N6+N7;
    rot2 = [1 0 0 ; 0 cos(angle02) -sin(angle02); 0 sin(angle02) cos(angle02)];
    
    %% Right
    
    
    if jr ==N1+N2+N3+N4+N5+N6+N7+1
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1-zoffset0.*0.8; coordi2_2 = coordi2_1-zoffset0.*0.8;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+2
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1-zoffset0.*0.2; coordi2_2 = coordi2_1-zoffset0.*0.2;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+3
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1+zoffset0.*0.4; coordi2_2 = coordi2_1+zoffset0.*0.4;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+4
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*rot2*coordi1; coordi2_1 = rot1*rot2*coordi2;
        coordi1_2 = coordi1_1+zoffset0.*0.2; coordi2_2 = coordi2_1+zoffset0.*0.2;
        
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+5
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1-zoffset0.*0.8; coordi2_2 = coordi2_1-zoffset0.*0.8;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+6
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1-zoffset0.*0.2; coordi2_2 = coordi2_1-zoffset0.*0.2;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+7
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1+zoffset0.*0.4; coordi2_2 = coordi2_1+zoffset0.*0.4;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+8
        angl   e1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*rot2*coordi1; coordi2_1 = rot1*rot2*coordi2;
        coordi1_2 = coordi1_1+zoffset0.*0.2; coordi2_2 = coordi2_1+zoffset0.*0.2;
        
        
        
        %% Left
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+9
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1-zoffset0.*0.8; coordi2_2 = coordi2_1-zoffset0.*0.8;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+10
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1-zoffset0.*0.2; coordi2_2 = coordi2_1-zoffset0.*0.2;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+11
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1+zoffset0.*0.4; coordi2_2 = coordi2_1+zoffset0.*0.4;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+12
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*rot2*coordi1; coordi2_1 = rot1*rot2*coordi2;
        coordi1_2 = coordi1_1+zoffset0.*0.2; coordi2_2 = coordi2_1+zoffset0.*0.2;
        
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+13
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1-zoffset0.*0.8; coordi2_2 = coordi2_1-zoffset0.*0.8;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+14
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1-zoffset0.*0.2; coordi2_2 = coordi2_1-zoffset0.*0.2;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+15
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*coordi1; coordi2_1 = rot1*coordi2;
        coordi1_2 = coordi1_1+zoffset0.*0.4; coordi2_2 = coordi2_1+zoffset0.*0.4;
        
    elseif jr ==N1+N2+N3+N4+N5+N6+N7+16
        angle1 = angle01(jr-nnn);
        rot1 = [cos(angle1) -sin(angle1) 0; sin(angle1) cos(angle1) 0; 0 0 1];
        coordi1_1 = rot1*rot2*coordi1; coordi2_1 = rot1*rot2*coordi2;
        coordi1_2 = coordi1_1+zoffset0.*0.2; coordi2_2 = coordi2_1+zoffset0.*0.2;
    end
    
    
    
    x1f = reshape(coordi1_2(1,:),[1 size(coordi1_1,2)]);
    y1f = reshape(coordi1_2(2,:),[1 size(coordi1_1,2)]);
    z1f = reshape(coordi1_2(3,:),[1 size(coordi1_1,2)]);
    
    x2f = reshape(coordi2_2(1,:),[1 size(coordi1_1,2)]);
    y2f = reshape(coordi2_2(2,:),[1 size(coordi1_1,2)]);
    z2f = reshape(coordi2_2(3,:),[1 size(coordi1_1,2)]);
    
    dlx = x2f-x1f;
    dly = y2f-y1f;
    dlz = z2f-z1f;
    
%     for k = 1:360
%         rx = x - (x1f(k)+x2f(k))/2;
%         ry = y - (y1f(k)+y2f(k))/2;
%         rz = z - (z1f(k)+z2f(k))/2;
%         r3 = sqrt(rx.^2 + ry.^2 + rz.^2).^3;
%         Bz(:,:,:,jr) = Bz(:,:,:,jr) + ...
%             1e-7*I*(dlx(k)*ry - dly(k)*rx) ./ r3;
%     end
    %
    %   figure(100), plot3(x1f,y1f,z1f,'.'), axis([-0.3 0.3 -0.3 0.3 -0.3 0.3]), hold on
    %   plot3(x1,y1,z1,'r.')
    %   set(gcf,'Color',[1 1 1])
    
    cordf(jr,1,:) = x1f;
    cordf(jr,2,:) = y1f;
    cordf(jr,3,:) = z1f;
    
    cordf1(jr,1,:) = x2f;
    cordf1(jr,2,:) = y2f;
    cordf1(jr,3,:) = z2f;
    
    jr
    
end



%% Up

[nx ny nz nc]= size(Bz);
nxi = round(nx/2)-round(nx0/2); nxf = round(nx/2)+round(nx0/2)-1;
nyi = round(ny/2)-round(ny0/2); nyf = round(ny/2)+round(ny0/2)-1;
nzi = round(nz/2)-round(nz0/2); nzf = round(nz/2)+round(nz0/2)-1;

rhh = -0.03;
offsetz = round(rhh/th)-1;

rhy = 0.015;
offsety =  round(rhy/res)-1;

% check with figure (13) and figure (14)
Bz_f = Bz(nxi-offsety:nxf-offsety,nyi:nyf,nzi-offsetz:nzf-offsetz,1:N);
%%
figure(100)

R = rr;
[X,Y,Z] = cylinder(R,100);
Z = (Z - 1)*height/2+0.01;
test1 = surf(X,Y,Z,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
set(test1,'FaceColor',[0 0 0],'FaceAlpha',0.3)
hold on

% rh = 0.06;

R1=rr;
[X1,Y1,Z1] = sphere;
Z1(Z1<-0.001) = NaN;
X1 = X1*R1; Y1 = Y1*R1; Z1 = Z1*R1;
% Z1(Z1> rh) = NaN;
test11 = surf(X1,Y1,Z1,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
set(test11 ,'FaceColor',[0 0 0],'FaceAlpha',0.3)
hold on

% load('B0_afterHOS_m.mat')
% B0_0 = logical(abs((B0)));
% 
% for j=1:size(B0_0,3)
%     B0_1(:,:,j) = rot90(B0_0(:,:,j),2);
% end
% [x00 y00 z00] = size(B0_1);
% [x11,y11,z11] = ind2sub(size(B0_1),find(B0_1 ~= 0));
% 
% x11f = (x11-round(x00/2))*res;
% y11f = (y11-round(y00/2))*res-rhy;
% z11f = (z11-round(z00/2))*th-rhh;
% 
% test111 = plot3(x11f,y11f,z11f,'b.','MarkerSize',6)
% % set(test111 ,'FaceColor',[0 0 0],'FaceAlpha',0.4)
% hold on
list = [5,8,9,14,15];
N = 51;
cordf(list,:,:) = [];
for k = 1:N
    
%     if k==5 | k==8 | k==9 | k==14 | k==15
%         continue;
%         
%     elseif k < N1+1;
%         color1 = [1 1 0];
%     elseif k < N1+N2+1
%         color1 = [1 0.5 0];
%     elseif k < N1+N2+N3+1
%         color1 = [0.9 0.3 0];
%     elseif k < N1+N2+N3+N4+1
%         color1 = [0.3 1 0.1];
%     elseif k < N1+N2+N3+N4+N5+1
%         color1 = [0.4 1 0.2];
%     elseif k < N1+N2+N3+N4+N5+N6+1
%         color1 = [0.2 0.9 0.1];
%     elseif k < N1+N2+N3+N4+N5+N6+N7+1
%         color1 = [1 0 1];
%     elseif k < N1+N2+N3+N4+N5+N6+N7+N8+1
%         color1 = [0 1 1];
%     end
    
    xx1f = reshape(cordf(k,1,:),[1 360]);
    yy1f = reshape(cordf(k,2,:),[1 360]);
    zz1f = reshape(cordf(k,3,:),[1 360]);
    
%     plot3(xx1f,yy1f,zz1f,'.-','Color',color1), axis([-range*1.1/2 range*1.1/2 -range*1.1/2 range*1.1/2 -range*1.1/2 range*1.1/2]), hold on
%     
        plot3(xx1f,yy1f,zz1f,'.-'), axis([-range*1.1/2 range*1.1/2 -range*1.1/2 range*1.1/2 -range*1.1/2 range*1.1/2]), hold on
    
end
hold on 

cordf_center = mean(cordf,3);
for p =1:size(cordf_center,1)
    xc = cordf_center(p,1);
    yc = cordf_center(p,2);
    zc = cordf_center(p,3);
    
    text(xc,yc,zc,num2str(p),'FontSize',12,'Color',[1 1 1])
end

grid on
set(gcf,'Color',[1 1 1])
toc
%%

toc

% save('Bz.mat','Bz')
% % save('Bz_f.mat','Bz_f')
% save('cordf.mat','cordf')

%% (2) 6cm frontal part, 8 of them
%
% x11 = rr/res.*cos((0:359)*pi/180)+round(nx/2);
% y11 = rr/res.*sin((0:359)*pi/180)+round(nx/2);
% 
% % save(strcat('Bz_rc',num2str(rc*1000),'_rr',num2str(rr*1000)),'Bz');
% m = 40; % upper & lower limits for the B0 maps (Hz)
% figure(11);
% for jr = 1:N
%     subplot(4,8,jr)
%     imshow(Bz(:,:,round(0.5*(height/th)),jr)*42.57e6,[-m m]),  colormap jet
%     hold on
%     plot(x11,y11,':','Color',[0.5 0 0.9],'LineWidth',1.2)
% %     set(f ,'FaceColor',[0 0 0],'FaceAlpha',0.4,'edgealpha',0.4)
% %         mesh(x11,y11)
% end
% set(gcf,'color',[1 1 1])

