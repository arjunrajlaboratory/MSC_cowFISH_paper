 %% Crop and merge images for Figure 1 showing multiple FISH in same cell.

 % NOTE: Images originally from chondrocyte dedifferation core, donor 4,
 % passage 3, experimentid G30, imaged 2013/10/31, folder /2013-10-31/D4
 % (see metadata spreadsheet for additional details)
 
 
img_tmr = readmm('./raw_images/tmr002.tif');
img_alexa = readmm('./raw_images/alexa002.tif');
img_cy = readmm('./raw_images/cy002.tif');
img_nir = readmm('./raw_images/nir002.tif');
img_dapi = readmm('./raw_images/dapi002.tif');



R = [504 51 425 350];

imgmax_alexa = max(img_alexa.imagedata(:,:,:),[],3);  %MAX MERGE
imgmax_cy = max(img_cy.imagedata(:,:,:),[],3);  %MAX MERGE
imgmax_tmr = max(img_tmr.imagedata(:,:,:),[],3);  %MAX MERGE
imgmax_dapi = max(img_dapi.imagedata(:,:,:),[],3);  %MAX MERGE
imgmax_nir = max(img_nir.imagedata(:,:,:),[],3);  %MAX MERGE

cy1 = imcrop(imgmax_cy,R);
al1 = imcrop(imgmax_alexa,R);
tm1 = imcrop(imgmax_tmr,R);
nir1 = imcrop(imgmax_nir,R);
dp1 = imcrop(imgmax_dapi,R);



% mn = min([tm1(:);cy1(:)]);
% mx = max([tm1(:);cy1(:)]);  
al1 = scale(al1);
dp1 = scale(dp1);
nir1 = scale(nir1)*.7;
tm1 = scale(tm1)*2;%,'intensityScale',[mn mx]);  % Pick same scales for both TMR and CY for a better apples to apples comparison
cy1 = scale(cy1)*2;

RGB = cat(3,nir1 + dp1,nir1,nir1 + dp1);
imwrite(im2uint16(RGB),'Fig1_GAPDH.tif')

RGB = cat(3,al1 + dp1,al1,al1 + dp1);
imwrite(im2uint16(RGB),'Fig1_LPL.tif')

RGB = cat(3,tm1 + dp1,tm1,tm1 + dp1);
imwrite(im2uint16(RGB),'Fig1_OPN.tif')

RGB = cat(3,cy1 + dp1,cy1,cy1 + dp1);
imwrite(im2uint16(RGB),'Fig1_ACAN.tif')

