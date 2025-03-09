function [u v] = warp_flow_Nitish(im1,im2,pyr_lev,warp_num,tempmat,beta,lmbda,wndcff)
%THIS FUNCTION TAKES INPUT IMAGES AS IM1 AND IM2, THE NUMBER OF
%PYRAMID LEVELS PYR_LEV, THE NUMBER OF WARPING TIMES IN EACH
%PYRAMID LEVEL GIVEN BY WARP_NUM.
if size(im1,3)>1
    im1 = rgb2gray(im1);
    im2 = rgb2gray(im2);
end
if ~isfloat(im1)
    im1 = double(im1);
    im2 = double(im2);
end
[mNew,nNew] = size(im1);%[m n]=[384 384]
for lev = pyr_lev:-1:1%pyr_lev=5
    tmp_im1 = imresize(im1,[ceil(0.5^(lev-1)*mNew) ceil(0.5^(lev-1)*nNew)]);
    tmp_im2 = imresize(im2,[ceil(0.5^(lev-1)*mNew) ceil(0.5^(lev-1)*nNew)]);
%     lev 
%     size(tmp_im2)
%     size(tmp_im1)
    [new_m ,new_n] = size(tmp_im1);
    if lev == pyr_lev
        u = zeros(new_m,new_n);
        v = zeros(new_m,new_n);
    end
    %lmbda = 100;%imresize(lmbda,[new_m new_n]);
    u = imresize(u,[new_m new_n]);
    v = imresize(v,[new_m new_n]);
    delta_u = zeros(new_m,new_n);
    delta_v = delta_u;
    tmp_im1 = imwarp(tmp_im1,-u,-v,'true');
    for warping_number = 1:1:warp_num
        disp(size(u))
        [fx,fy,ft] = computeDerivatives(tmp_im1,tmp_im2);
        %[row ,col] = size(tmp_im1);
        fxN = padarray(fx,[1 1],'both');
        fyN = padarray(fy,[1 1],'both');
        ftN = padarray(ft,[1 1],'both');
        IN =  padarray(tmp_im1,[1 1],'both');
        
        fxN2 = fxN.*fxN;
        fyN2 = fyN.*fyN;
        fxyN2 = fxN.*fyN;
        fxtN2 = fxN.*ftN;
        fytN2 = fyN.*ftN;
        
        [m,n] = size(fxN);
        M = m-2;
        N = n-2;
        
        fxN2E = zeros(M,N);
        fyN2E = zeros(M,N);
        fxyN2E = zeros(M,N);
        fxtN2E = zeros(M,N);
        fytN2E = zeros(M,N);
        for i = 2:1:m-1
            for j = 2:1:n-1
                sumfx = 0;
                sumfy = 0;
                sumfxy = 0;
                sumfxt = 0;
                sumfyt = 0;
                DSum = 0;
                IVal = IN(i,j);
                for r = -1:1:1
                    for t = -1:1:1
                        sumfx = sumfx + exp(-abs(IN(i+r,j+t) - IVal))*fxN2(i+r,j+t);
                        sumfy = sumfy + exp(-abs(IN(i+r,j+t) - IVal))*fyN2(i+r,j+t);
                        sumfxy = sumfxy + exp(-abs(IN(i+r,j+t) - IVal))*fxyN2(i+r,j+t);
                        sumfxt = sumfxt + exp(-abs(IN(i+r,j+t) - IVal))*fxtN2(i+r,j+t);
                        sumfyt = sumfyt + exp(-abs(IN(i+r,j+t) - IVal))*fytN2(i+r,j+t);
                        DSum = DSum + exp(-abs(IN(i+r,j+t) - IVal));
                    end
                end
                fxN2E(i-1,j-1) = sumfx/DSum;
                %max(max(fxN2E))
                fyN2E(i-1,j-1) = sumfy/DSum;
                %max(max(fyN2E))
                fxyN2E(i-1,j-1) = sumfxy/DSum;
                %max(max(fxyN2E))
                fxtN2E(i-1,j-1) = sumfxt/DSum;
                %max(max(fxtN2E))
                fytN2E(i-1,j-1) = sumfyt/DSum;
                %max(max(fytN2E))
            end
        end
        for ite = 1:1:1500
            uAvg=conv2(delta_u,tempmat,'same');
            vAvg=conv2(delta_v,tempmat,'same');
            D1=fxN2E+beta^2+lmbda*wndcff;
            D2=fyN2E+beta^2+lmbda*wndcff;
            delta_u=(lmbda.*uAvg-fxyN2E.*delta_v-fxtN2E)./D1;
            delta_v=(lmbda.*vAvg-fxyN2E.*delta_u-fytN2E)./D2;

        end
        delta_u = medfilt2(delta_u,[5 5]);
        delta_v = medfilt2(delta_v,[5 5]);
        u = u + delta_u;
        v = v + delta_v;
        u = medfilt2(u,[5 5]);
        v = medfilt2(v,[5 5]);
        tmp_im1 = imwarp(tmp_im1,-delta_u,-delta_v,'true');
        delta_u = 0*delta_u;
        delta_v = 0*delta_v;
        
    end
end
%figure, imshow(imwarp(im1,u,v,'true')-im2,[])
end