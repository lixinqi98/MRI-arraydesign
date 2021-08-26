function [Bz] = BiotSavart100Points(coilcordf,N_coils,x,y,z,current)
% generate the Bz for each coil
% ensure the 100 points in sequence
% Created. 24/09/2019 Mona
% 
% Params: 
%       @coilcordf: the position of coils
%       @N_coils  : number of coils
%       @x        : 
%       @y        :
%       @z        :
%       @current  : 

    Bzc = double(zeros([size(x,1) size(x,2) size(x,3)]));
    [dim, nums, n] = size(coilcordf);
    tic
    h=waitbar(0,'please wait');
    for i = 1:N_coils
        % attention the coordinate
        x1f = coilcordf(1,:,i);
        y1f = coilcordf(3,:,i);
        z1f = coilcordf(2,:,i);
        I = current;


        for k = 1:nums
            if k == nums
               m = 1;
            else
                m = k+1;
            end
            dlx = x1f(m)-x1f(k);
            dly = y1f(m)-y1f(k);

            rx = x - (x1f(k)+x1f(m))/2;
            ry = y - (y1f(k)+y1f(m))/2;
            rz = z - (z1f(k)+z1f(m))/2;
            rsq = sqrt(rx.^2 + ry.^2 + rz.^2);
            %for faster implementation
            r3=rsq.*rsq.*rsq;
            %sumB(k,:) = 1e-7*I*(dlx(k)*ry - dly(k)*rx)./r3;
            Bzc = Bzc+ 1e-7*I*(dlx*ry - dly*rx)./r3;
        end

        Bz(:,:,:,i)=Bzc;
        Bzc = double(zeros([size(x,1) size(x,2) size(x,3)]));
        str=['processing...',num2str(i/N_coils*100),'%'];
        waitbar(i/N_coils,h,str)
    end
    delete(h)
    toc
end
