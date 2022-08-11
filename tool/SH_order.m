function [B_collect, S_collect, title_collect] = SH_order(order,x,y,z)
% ======================================================================= %
% SH_order.m
% ======================================================================= %
% This function is going to construct the Bz (coils) field for shimming 
% based on the selected order. Then ouput the B field and the sensitivity
% of all coils.
%
% Input: 
% =======
%   order      : The order of SH shimming, it can be 2 (9 coils), 3 (16 coils), 
%                4 (25 coils), 5 (36 coils), 6 (49 coils).              
%
% Output:
% =======
%   B_collect  : Bz (magnetic field of each coils) in a cell format with
%                size (# of coils,1)
%   s_collect  : sensitivity of all coils in a vector format with size 
%                (# of coils,1)
%
% Reference:
% ==========
% Paper: Modeling Real Shim Fields for Very High Degree (and Order) B0 Shimming 
% of the Human Brain at 9.4 T. (Paul Chang,1,2y* Sahar Nassirpour,1,2y and 
% Anke Henning)
% ======================================================================= %

% ------------------------- 0th order Bz in SH -------------------------- %
Bo_Z0 = 1;     % shimming term/symmetry: Z0
s0 = 6058;     % sensitivity of Bo (Hz/cm^n/A)

% ------------------------- 1st order Bz in SH -------------------------- %
B1_Z = z;  s1_Z = 1;   % shimming term (cm) & sensitivity (Hz/cm^n/A)/symmetry: Z
B1_X = x;  s1_X = 1;   % shimming term (cm) & sensitivity (Hz/cm^n/A)/symmetry: X
B1_Y = y;  s1_Y = 1;   % shimming term (cm) & sensitivity (Hz/cm^n/A)/symmetry: Y

% ------------------------- 2nd order Bz in SH -------------------------- %       
B2_Z2 = z.^2 - (1/2).*(x.^2 + y.^2);    s2_Z2 = 6.942;   % (cm); (Hz/cm^n/A)
B2_ZX = z .* x;                         s2_ZX = 24.15;   % (cm); (Hz/cm^n/A)
B2_ZY = z .* y;                         s2_ZY = 24.15;   % (cm); (Hz/cm^n/A)
B2_C2 = x.^2 - y.^2;                    s2_C2 = 3.64;    % (cm); (Hz/cm^n/A)
B2_S2 = 2 * x .* y;                     s2_S2 = 3.64;    % (cm); (Hz/cm^n/A)

% ------------------------- 3rd order Bz in SH -------------------------- %
B3_Z3 = z.^3 - 3/2.*z.*(x.^2 + y.^2);      s3_Z3 = 0.4923;  % (cm); (Hz/cm^n/A)
B3_Z2X = z.^2.*x - 1/4*x.*(x.^2 + y.^2);   s3_Z2X = 1.0;    % (cm); (Hz/cm^n/A)
B3_Z2Y = z.^2.*y - 1/4*y.*(x.^2 + y.^2);   s3_Z2Y = 1.0;    % (cm); (Hz/cm^n/A)
B3_ZC2 = z .* (x.^2 - y.^2);               s3_ZC2 = 1.77;   % (cm); (Hz/cm^n/A)
B3_ZS2 = z .* (2*x.*y);                    s3_ZS2 = 1.77;   % (cm); (Hz/cm^n/A)
B3_C3 = x.^3 - 3*x.*y.^2;                  s3_C3 = 0.188;   % (cm); (Hz/cm^n/A)
B3_S3 = 3*x.^2.*y - y.^3;                  s3_S3 = 0.188;   % (cm); (Hz/cm^n/A)

% ------------------------- 4th order Bz in SH -------------------------- %
B4_Z4 = z.^4 - 3*z.^2.*(x.^2 + y.^2) + 3/8.*(x.^2 + y.^2).^2; s4_Z4 = 0.04206; % (cm); (Hz/cm^n/A)
B4_Z3X = z.^3.*x - 3/4*z.*x.*(x.^2 + y.^2);                   s4_Z3X = 0.123;  % (cm); (Hz/cm^n/A)
B4_Z3Y = z.^3.*y - 3/4*z.*y.*(x.^2 + y.^2);                   s4_Z3Y = 0.123;  % (cm); (Hz/cm^n/A)
%B4_Z2C2 = z.*(x.^2 - y.^2).*(z.^2 - 1/6*(x.^2 + y.^2));       s4_Z2C2 = 0.093; % (cm); (Hz/cm^n/A)
%B4_Z2S2 = 2*z.* (x.*y).*(z.^2 - 1/6*(x.^2 + y.^2));           s4_Z2S2 = 0.093; % (cm); (Hz/cm^n/A)
B4_Z2C2 = (x.^2 - y.^2).*(z.^2 - 1/6*(x.^2 + y.^2));          s4_Z2C2 = 0.093; % (cm); (Hz/cm^n/A)
B4_Z2S2 = 2*(x.*y).*(z.^2 - 1/6*(x.^2 + y.^2));               s4_Z2S2 = 0.093; % (cm); (Hz/cm^n/A)
%B4_ZC3 = x.* (x.^2-3*y.^2).*(z.^2 - 1/8*(x.^2 + y.^2));       s4_ZC3 = 0.121;  % (cm); (Hz/cm^n/A)
%B4_ZS3 = x.* (x.^2-3*y.^2).*(z.^2 - 1/8*(x.^2 + y.^2));       s4_ZS3 = 0.121;  % (cm); (Hz/cm^n/A)
B4_ZC3 = z .* (x.^3 - 3*x.*y.^2);                             s4_ZC3 = 0.121;  % (cm); (Hz/cm^n/A)
B4_ZS3 = z .* (3*x.^2.*y - y.^3);                             s4_ZS3 = 0.121;  % (cm); (Hz/cm^n/A)
B4_C4 = x.^4 - 6*x.^2.*y.^2 + y.^4;                           s4_C4 = 0.0187;  % (cm); (Hz/cm^n/A)
B4_S4 = 4 * (x.^3.*y - y.^3.*x);                              s4_S4 = 0.0187;  % (cm); (Hz/cm^n/A)

% ------------------------- 5th order Bz in SH -------------------------- %
B5_ZC4 = z.* (x.^4 - 6*x.^2.*y.^2 + y.^4);      s5_ZC4 = 5.71e-3;  % (cm); (Hz/cm^n/A)
B5_ZS4 = 4*z.* (x.^3.*y - y.^3.*x);             s5_ZS4 = 5.71e-3;  % (cm); (Hz/cm^n/A)
B5_C5 = x.* (x.^4 - 10*x.^2.*y.^2 + 5*y.^4);    s5_C5 = 9.9e-4;    % (cm); (Hz/cm^n/A)
B5_S5 = y.* (y.^4 - 10*x.^2.*y.^2 + 5*x.^4);    s5_S5 = 9.9e-4;    % (cm); (Hz/cm^n/A)
% extra term not in the paper
B5_Z5 = z.^5 - 5*z.^3.*(x.^2 + y.^2) + 15/8*z.*(x.^2 + y.^2).^2;           s5_Z5 = 5.71e-3;  % (cm); (Hz/cm^n/A)
B5_Z4X = z.^4.*x - 3/2*z.^2.*x.*(x.^2 + y.^2) + 1/8*x.*(x.^2 + y.^2).^2;   s5_Z4X = 5.71e-3;  % (cm); (Hz/cm^n/A)
B5_Z4Y = z.^4.*y - 3/2*z.^2.*y.*(x.^2 + y.^2) + 1/8*y.*(x.^2 + y.^2).^2;   s5_Z4Y = 5.71e-3;  % (cm); (Hz/cm^n/A)
B5_Z3C2 = z.^3.*(x.^2 - y.^2) - 1/2*z.*(x.^2 - y.^2).*(x.^2 + y.^2);       s5_Z3C2 = 5.71e-3;  % (cm); (Hz/cm^n/A)
B5_Z3S2 = z.^3 .* (2*x.*y) - 1/2*z.*(2*x.*y).*(x.^2 + y.^2);               s5_Z3S2 = 5.71e-3;  % (cm); (Hz/cm^n/A)
B5_Z2C3 = z.^2 .* (x.^3 - 3*x.*y.^2) - 1/8*(x.^3 - 3*x.*y.^2).*(x.^2 + y.^2);    s5_Z2C3 = 5.71e-3;  % (cm); (Hz/cm^n/A)
B5_Z2S3 = z.^2 .* (3*x.^2.*y - y.^3) - 1/8*(3*x.^2.*y - y.^3).*(x.^2 + y.^2);    s5_Z2S3 = 5.71e-3;  % (cm); (Hz/cm^n/A)

% ------------------------- 6th order Bz in SH -------------------------- %
B6_ZC5 = z.* x.* (x.^4 - 10*x.^2.*y.^2 + 5*y.^4);    s6_ZC5 = 3.21e-4; % (cm); (Hz/cm^n/A)
B6_ZS5 = z.* y.* (y.^4 - 10*x.^2.*y.^2 + 5*x.^4);    s6_ZS5 = 3.21e-4; % (cm); (Hz/cm^n/A)
% remaining terms
B6_Z6 = z.^6 - 15/2*z.^4.*(x.^2+y.^2) + 45/8*z.^2.*(x.^2+y.^2).^2 - 5/16*(x.^2+y.^2).^3; s6_Z6 = 3.21e-4; % (Hz/cm^n/A)
B6_Z5X = z.*x.*(z.^4 - 5/2*z.^2.*(x.^2+y.^2) + 5/8*(x.^2+y.^2).^2);  s6_Z5X = 3.21e-4; % (Hz/cm^n/A)
B6_Z5Y = z.*y.*(z.^4 - 5/2*z.^2.*(x.^2+y.^2) + 5/8*(x.^2+y.^2).^2);  s6_Z5Y = 3.21e-4; % (Hz/cm^n/A)
B6_Z4C2 = (x.^2-y.^2).*(z.^4-z.^2.*(x.^2+y.^2)+1/16*(x.^2+y.^2).^2); s6_Z4C2 = 3.21e-4; % (Hz/cm^n/A)
B6_Z4S2 = (2*x.*y).*(z.^4-z.^2.*(x.^2+y.^2)+1/16*(x.^2+y.^2).^2);    s6_Z4S2 = 3.21e-4; % (Hz/cm^n/A)
B6_Z3C3 = x.*z.*(x.^2-3*y.^2).*(z.^2-3/8*(x.^2+y.^2));               s6_Z3C3 = 3.21e-4; % (Hz/cm^n/A)
B6_Z3S3 = y.*z.*(3*x.^2-y.^2).*(z.^2-3/8*(x.^2+y.^2));               s6_Z3S3 = 3.21e-4; % (Hz/cm^n/A)
B6_Z2C4 = z.^2.*(x.^2-y.^2).^2-x.^2.*y.^2.*(4*z.^2-1/2*(x.^2+y.^2))-1/10*(x.^6+y.^6);  s6_Z2C4 = 3.21e-4; % (Hz/cm^n/A)
B6_Z2S4 = 4*x.*y.*(x.^2-y.^2).*(z.^2-1/10*(x.^2+y.^2));              s6_Z2S4 = 3.21e-4; % (Hz/cm^n/A)
B6_C6 = x.^6-15*x.^2.*y.^2.*(x.^2-y.^2)-y.^6;                        s6_C6 = 3.21e-4; % (Hz/cm^n/A)
B6_S6 = (2.*x.*y).*(3*x.^2-y.^2).*(x.^2-3*y.^2);                     s6_S6 = 3.21e-4; % (Hz/cm^n/A)



% ------------ collect B field and sensitivity of each coil ------------- %
switch order
    case 2
        B_collect = {Bo_Z0; B1_Z; B1_X; B1_Y; B2_Z2; B2_ZX; B2_ZY; B2_C2; B2_S2};
        S_collect = [s0; s1_Z; s1_X; s1_Y; s2_Z2; s2_ZX; s2_ZY; s2_C2; s2_S2];
        title_collect = {'Z0';'Z';'X';'Y';'Z2';'ZX';'ZY';'C2';'S2'};
    case 3
        B_collect = {Bo_Z0; B1_Z; B1_X; B1_Y; B2_Z2; B2_ZX; B2_ZY; B2_C2; B2_S2;...
                      B3_Z3; B3_Z2X; B3_Z2Y; B3_ZC2; B3_ZS2; B3_C3; B3_S3}; 
        S_collect = [s0; s1_Z; s1_X; s1_Y; s2_Z2; s2_ZX; s2_ZY; s2_C2; s2_S2;...
                      s3_Z3; s3_Z2X; s3_Z2Y; s3_ZC2; s3_ZS2; s3_C3; s3_S3];
        title_collect = {'Z0';'Z';'X';'Y';'Z2';'ZX';'ZY';'C2';'S2';...
                      'Z3';'Z2X';'Z2Y';'ZC2';'ZS2';'C3';'S3'};
    case 4
        B_collect = {Bo_Z0; B1_Z; B1_X; B1_Y; B2_Z2; B2_ZX; B2_ZY; B2_C2; B2_S2;...
                     B3_Z3; B3_Z2X; B3_Z2Y; B3_ZC2; B3_ZS2; B3_C3; B3_S3; B4_Z4; B4_Z3X;...
                     B4_Z3Y; B4_Z2C2; B4_Z2S2; B4_ZC3; B4_ZS3; B4_C4; B4_S4}; 
        S_collect = [s0; s1_Z; s1_X; s1_Y; s2_Z2; s2_ZX; s2_ZY; s2_C2; s2_S2;...
                     s3_Z3; s3_Z2X; s3_Z2Y; s3_ZC2; s3_ZS2; s3_C3; s3_S3; s4_Z4; s4_Z3X;...
                     s4_Z3Y; s4_Z2C2; s4_Z2S2; s4_ZC3; s4_ZS3; s4_C4; s4_S4];
        title_collect = {'Z0';'Z';'X';'Y';'Z2';'ZX';'ZY';'C2';'S2';...
                      'Z3';'Z2X';'Z2Y';'ZC2';'ZS2';'C3';'S3';'Z4';'Z3X';...
                      'Z3Y';'Z2C2';'Z2S2';'ZC3';'ZS3';'C4';'S4'};
    case 5
        B_collect = {Bo_Z0; B1_Z; B1_X; B1_Y; B2_Z2; B2_ZX; B2_ZY; B2_C2; B2_S2;...
                     B3_Z3; B3_Z2X; B3_Z2Y; B3_ZC2; B3_ZS2; B3_C3; B3_S3; B4_Z4; B4_Z3X;...
                     B4_Z3Y; B4_Z2C2; B4_Z2S2; B4_ZC3; B4_ZS3; B4_C4; B4_S4;...
                     B5_ZC4; B5_ZS4; B5_C5; B5_S5; B5_Z5; B5_Z4X; B5_Z4Y; B5_Z3C2;...
                     B5_Z3S2; B5_Z2C3; B5_Z2S3}; 
                   
        S_collect = [s0; s1_Z; s1_X; s1_Y; s2_Z2; s2_ZX; s2_ZY; s2_C2; s2_S2;...
                     s3_Z3; s3_Z2X; s3_Z2Y; s3_ZC2; s3_ZS2; s3_C3; s3_S3; s4_Z4; s4_Z3X;...
                     s4_Z3Y; s4_Z2C2; s4_Z2S2; s4_ZC3; s4_ZS3; s4_C4; s4_S4;...
                     s5_ZC4; s5_ZS4; s5_C5; s5_S5; s5_Z5; s5_Z4X; s5_Z4Y; s5_Z3C2;...
                     s5_Z3S2; s5_Z2C3; s5_Z2S3];
                 
        title_collect = {'Z0';'Z';'X';'Y';'Z2';'ZX';'ZY';'C2';'S2';...
                      'Z3';'Z2X';'Z2Y';'ZC2';'ZS2';'C3';'S3';'Z4';'Z3X';...
                      'Z3Y';'Z2C2';'Z2S2';'ZC3';'ZS3';'C4';'S4';...
                      'ZC4';'ZS4';'C5';'S5';'Z5';'Z4X';'Z4Y';'Z3C2';'Z3S2';...
                      'Z2C3';'Z2S3'};                  
    case 6
        B_collect = {Bo_Z0; B1_Z; B1_X; B1_Y; B2_Z2; B2_ZX; B2_ZY; B2_C2; B2_S2;...
                     B3_Z3; B3_Z2X; B3_Z2Y; B3_ZC2; B3_ZS2; B3_C3; B3_S3; B4_Z4; B4_Z3X;...
                     B4_Z3Y; B4_Z2C2; B4_Z2S2; B4_ZC3; B4_ZS3; B4_C4; B4_S4;...
                     B5_ZC4; B5_ZS4; B5_C5; B5_S5;B5_Z5; B5_Z4X; B5_Z4Y; B5_Z3C2;...
                     B5_Z3S2; B5_Z2C3; B5_Z2S3; B6_ZC5; B6_ZS5; B6_Z6; B6_Z5X; ...
                     B6_Z5Y; B6_Z4C2; B6_Z4S2; B6_Z3C3; B6_Z3S3; B6_Z2C4; B6_Z2S4;...
                     B6_C6; B6_S6}; 
                    
        S_collect = [s0; s1_Z; s1_X; s1_Y; s2_Z2; s2_ZX; s2_ZY; s2_C2; s2_S2;...
                     s3_Z3; s3_Z2X; s3_Z2Y; s3_ZC2; s3_ZS2; s3_C3; s3_S3; s4_Z4; s4_Z3X;...
                     s4_Z3Y; s4_Z2C2; s4_Z2S2; s4_ZC3; s4_ZS3; s4_C4; s4_S4;...
                     s5_ZC4; s5_ZS4; s5_C5; s5_S5; s5_Z5; s5_Z4X; s5_Z4Y; s5_Z3C2;...
                     s5_Z3S2; s5_Z2C3; s5_Z2S3; s6_ZC5; s6_ZS5; s6_Z6; s6_Z5X; ...
                     s6_Z5Y; s6_Z4C2; s6_Z4S2; s6_Z3C3; s6_Z3S3; s6_Z2C4; s6_Z2S4; ...
                     s6_C6; s6_S6];
                    
        title_collect = {'Z0';'Z';'X';'Y';'Z2';'ZX';'ZY';'C2';'S2';...
                      'Z3';'Z2X';'Z2Y';'ZC2';'ZS2';'C3';'S3';'Z4';'Z3X';...
                      'Z3Y';'Z2C2';'Z2S2';'ZC3';'ZS3';'C4';'S4';...
                      'ZC4';'ZS4';'C5';'S5';'Z5';'Z4X';'Z4Y';'Z3C2';'Z3S2';...
                      'Z2C3';'Z2S3';'ZC5';'ZS5';'Z6';'Z5X';'Z5Y';'Z4C2';...
                      'Z4S2';'Z3C3';'Z3S3';'Z2C4';'Z2S4';'C6';'S6'};
                     
end


