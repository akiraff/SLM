%% Offline trap uniformity tweak
%% Take a picture with a certain exposure time
exposure_time = 600; % in us
image2D = takePictureThorCam(cam, exposure_time);
img_data = transpose(double(image2D));
%% Get rid of zero order crap such that we can use our old script for peak finding and
% gaussian fitting.
minrow_cp1=200;
maxrow_cp1=900;
mincol_cp1=200;
maxcol_cp1=1000;
%%
data_crop0=(img_data(minrow_cp1:maxrow_cp1,mincol_cp1:maxcol_cp1));
%%
img_data_noZero = data_crop0;
figure(91);clf;
imagesc(img_data_noZero)
axis('equal');
%% check pattern spacing
xy0 = [178, 111];
xy=[284, 116];
%((547-469)^2+(275-195)^2).^0.5*3.45
(sum((xy-xy0).^2))^0.5*3.45
%% record where zero order is.
xo = 864;
yo = 260;
%% Now fit for the array intensity
level=3;
num_sites=36;
gain =0.01;
spacing = 379/3.45;
bGaussian = 1;
bmodguess = 0;
bmanual = 0; 
manualguess = 0;
trapArray=trapArraysCam_GS2D(img_data_noZero, level, num_sites,gain,spacing,bGaussian, bmodguess, bmanual, manualguess);
%%
save("img_data_noZero.mat","img_data_noZero");
%% now save trap array
save("trapArray.mat","trapArray");
%% Now start to sort the trap for uniformity tweak
cropArea=trapArray.cropArea;
minrow_cp2 = cropArea(1);
mincol_cp2 = cropArea(3);
%% More general method to find the end of the row. Note: row distribution is regular
trap=trapArray.Trap;
trap_sort = sortrows(trap,3,'ascend');
xloc=trap_sort(:,2)+mincol_cp1+mincol_cp2;
yloc=trap_sort(:,3)+minrow_cp1+minrow_cp2;
rloc=((xloc-xo).^2+(yloc-yo).^2).^0.5;
trap_sort_rloc=[trap_sort, rloc];
sz = size(trap_sort_rloc);
trap_sort_final=zeros(sz);
row = sz(1);
delta=379/3.45/2;
row_val = trap_sort_rloc(:,3);
row_end = [];
for j = 1:row
    row_cur = row_val(j);
    if j < row
        row_next = row_val(j+1);
        if row_next>row_cur+delta
            row_end = [row_end, j];
        end
    end
end
% populate first row
trap_row_1 = trap_sort_rloc(1:row_end(1),:);
trap_row_1_sort = sortrows(trap_row_1,4,'ascend');
trap_sort_final(1:row_end(1),:) = trap_row_1_sort;
for j = 1:length(row_end)
    row_startnow = row_end(j)+1;
    if j<length(row_end)
        row_endnow = row_end(j+1);
        trap_row=trap_sort_rloc(row_startnow: row_endnow,:);
        trap_row_sort = sortrows(trap_row, 4, 'ascend');
        trap_sort_final(row_startnow: row_endnow,:) =  trap_row_sort;
    elseif j==length(row_end)
        trap_row=trap_sort_rloc(row_startnow:end,:);
        trap_row_sort = sortrows(trap_row, 4, 'ascend');
        trap_sort_final(row_startnow: end,:) =  trap_row_sort;
    end   
end
%% normalize intensity for the sorted trap
trapInt = trap_sort_final(:,1);
trapInt_norm = trapInt/sum(trapInt);
%% Intensity histogram
h=figure(100);clf
histogram(trapInt_norm*36)
axis([0.1 2 0 15]);
nacstools.display.makePretty('figure',h, 'width', 10, 'height', 8, 'textFontSize', 10, 'axisFontSize', 10);
%%
csvwrite('trapInt_norm_R0_0217.csv', trapInt_norm)