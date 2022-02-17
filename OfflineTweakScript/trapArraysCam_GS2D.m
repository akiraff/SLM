classdef trapArraysCam_GS2D
    properties
        NumSites; % Number of sites in each array
        AOD_Amp0; % Initial Amplitude for each AOD tones, frequency ascending
        AOD_Amp; % Amplitude for each AOD tones after optimization
        AOD_Phase0; % Initial Phase for each AOD tones
        AOD_Phase; % Phase for each AOD tones after optimization
        Data; % 2 D matrix with value on each pixel
        Trap; % [peak,x,y] include information of peak values and location, y pos descend
        TrapArea; % This area value corresponds to the integrated area after 2D gaussian fitting
        TrapArea_Error; % This error is defined to be abs(max-min)/trapavg
        Trapsigx;
        Trapsigy;
        cropArea;
        TotalPower; % sum of the values on all trap pixels
        NonUniformity; % (sum((peakvalue-peakvalue_avg).^2)/num_sites).^0.5
        ErrorPercentNow; % Error for current array
        TrapTweak_Ind; % The trap index being tweaked
        Gain; % Gain for single step tweaking for amplitude adjustment
        bGaussian; % 1: do 2D gaussian. 0: use local maximum
    end
    
    methods
        function self = trapArraysCam_GS2D(imgData,level,num_sites,gain,spacing,bGaussian, bmodguess, bmanual, manualguess )
            % imgData: 2D matrix taken by andor
            % level: Define threshold for peak finder extrema and exrema2
            % num_sites: number of sites for optimization
            % gain: proportional gain for RF amp adjustment
           Data0 = imgData;
           maxdata=max(max(Data0));
           threshold = maxdata/level;
           [row,col]=find(Data0 >threshold);
           minrow=max(min(row)-2,1);
           maxrow=max(row)+2;
           mincol=max(min(col)-2,1);
           maxcol=max(col)+2;
          %disp([minrow,maxrow,mincol,maxcol]);
          self.cropArea = [minrow,maxrow,mincol,maxcol];
          data_crop=(Data0(minrow:maxrow,mincol:maxcol));
          self.Data = data_crop;
          self.bGaussian = bGaussian;
          figure(101);clf;
          data_surf=surf(data_crop);
          IMGX=data_surf.XData;
          IMGY=data_surf.YData;
          IMGZ=data_surf.ZData;
         [xMesh,yMesh]=meshgrid(IMGX,IMGY);
         [zmax,imax,~,~] = extrema2(IMGZ);
         hold on  
         plot3(xMesh(imax),yMesh(imax),zmax,'bo')
         %disp(imax)
        for i = 1:length(zmax)
           text(xMesh(imax(i)),yMesh(imax(i)),zmax(i),['  ' num2str(zmax(i))])
        end
        hold off 
        peak=[zmax,xMesh(imax),yMesh(imax)];
        peak_sort=sortrows(peak,'descend');
        self.NumSites = num_sites;
        % start with more sites and remove closely spaced sites
        sites_candi= min(floor(4.0*num_sites),length(zmax));
        trapInfo_mod = peak_sort(1:sites_candi,:);
        % sorting spaceing 1
        repeat = 1;
        sortingSpacing = 1;
        trapInfo_mod1 = peakFinderHelper(trapInfo_mod, repeat, spacing, sortingSpacing);
        % sorting spaceing 2
        sortingSpacing = 2;
        trapInfo_mod2 = peakFinderHelper(trapInfo_mod1, repeat, spacing, sortingSpacing);
        % sorting spacing 3
        sortingSpacing = 3;
       % disp(spacing);
        trapInfo_mod3 = peakFinderHelper(trapInfo_mod2, repeat, spacing, sortingSpacing);
        % sorting spacing 4
        sortingSpacing = 4;
        trapInfo_mod4 = peakFinderHelper(trapInfo_mod3, repeat, spacing, sortingSpacing);
        if bmodguess
           trapInfo_sort = sortrows(trapInfo_mod4,'descend');
        else
           trapInfo_sort = sortrows(trapInfo_mod,'descend');
        end
        if bmanual
           trapInfo =  manualguess;
        else
           trapInfo = trapInfo_sort(1:num_sites,:);
        end
        
        trapInfo_pos = sortrows(trapInfo,3,'descend');% y position down to up, frequency down to up
        
        
        peak_sites = trapInfo_pos(:,1);
        xpeak_sites = trapInfo_pos(:,2);
        ypeak_sites = trapInfo_pos(:,3);
        figure(102);clf;
       % disp(isites)
        plotc(xpeak_sites,ypeak_sites,peak_sites,'d');
        
        if bGaussian
            sz=size(data_crop);
            x=1:sz(2);
            y=1:sz(1);
            [X,Y] = meshgrid(x,y);
            xy=[X(:),Y(:)];
            % 2D gaussian Fit
            coef_guess = initial_guess(trapInfo_pos, num_sites,0);
            [coef, ~] = lsqcurvefit(@gauss2d, coef_guess, xy, data_crop(:), [], [], optimset('Display','off'));
            pks=(coef(1:6:end-1))';
            xpos=(coef(2:6:end-1))';
            ypos=(coef(3:6:end-1))';
            sigx=(coef(4:6:end-1))';
            sigy=(coef(5:6:end-1))';
            trapArea = 2*pks*pi.*sigx.*sigy;
            trapArea_avg = sum(trapArea)/length(trapArea);
            trapArea_err = (max(trapArea)-min(trapArea))/trapArea_avg;
            self.TrapArea = trapArea;
            self.TrapArea_Error = trapArea_err;
            self.Trapsigx = sigx;
            self.Trapsigy = sigy;
            trapInfo_pos=[pks,xpos,ypos];
            peak_sites = pks;
            z=gauss2d(coef, xy);
            gaussFit= reshape(z,size(data_crop));
            % check fitting
            figure(103);clf;
            subplot(1,2,1);
            imagesc(gaussFit);
            title('Gauss Fit');
            subplot(1,2,2);
            imagesc(data_crop);
            title('Raw image');
            data_crop_integral=trapz(data_crop,2);
            gaussFit_integral=trapz(gaussFit,2);
            figure(104);clf;
            plot(y,data_crop_integral,'ro',y,gaussFit_integral,'r');
            legend('data','Fit');
        end
        self.Trap = trapInfo_pos;
        self.TotalPower = sum(peak_sites);
        peak_sites_avg = sum(peak_sites)/num_sites;
        self.NonUniformity = (sum ((peak_sites-peak_sites_avg).^2)./num_sites).^0.5;
        % calculate current error percent for trap uniformity
           int_now = trapInfo_pos(:,1);
           [int_min, ind_min] = min(int_now);
            [int_max, ind_max] = max(int_now);
            int_mean = mean(int_now);
            if abs(int_min-int_mean) > abs(int_max-int_mean)
                int_goal = int_max;
                ind_tweak = ind_min;
            else
                int_goal = int_min;
                ind_tweak = ind_max;
            end
            self.ErrorPercentNow = -(int_goal - int_now(ind_tweak)) / int_goal; 
            self.TrapTweak_Ind = ind_tweak;
            self.Gain = gain;
        end
     
        
        function self=init_AODAmp(self,Amp)
            self.AOD_Amp0 = Amp;
        end
        
        function self=init_AODPhase(self,Phase)
            self.AOD_Phase0 = Phase;
        end
        
        function self = updateArrayAmp(self)
              amp_perc = self.AOD_Amp0/sum(self.AOD_Amp0);
              ap_tmp = amp_perc;   
              err_perc = self.ErrorPercentNow;
              delta = -err_perc * self.Gain;
              ind_tweak = self.TrapTweak_Ind;
              ap_tmp(ind_tweak) = amp_perc(ind_tweak) + delta;
              ap_tmp([1:(ind_tweak-1),(ind_tweak+1):end]) ...
                    = ap_tmp([1:(ind_tweak-1),(ind_tweak+1):end]) - delta / (self.NumSites-1);
              amp_perc = ap_tmp;
              self.AOD_Amp = sum(self.AOD_Amp0)*amp_perc;
              self.AOD_Phase = self.AOD_Phase0;
        end
        
        function self = updateArrayPhase(self, Phase)
              self.AOD_Phase = Phase;
        end   
    end
end

function coef0 = initial_guess(trapdata, n, offset)
% trapdata = trapArray.Trap
xw = 5*ones([n,1]);
pks = trapdata(:,1);
x_locs = trapdata(:,2);
y_locs = trapdata(:,3);
coef0 = [pks, x_locs, y_locs, xw, xw, zeros(size(pks))];
coef0 = reshape(coef0', [1, numel(coef0)]);
coef0(end+1)=offset;
end
function z = gauss2d0(coef, xy)
A = coef(1);
x0 = coef(2);
y0 = coef(3);
sigx = coef(4);
sigy = coef(5);
theta = coef(6);

a = 1/2 * ((cos(theta)/sigx)^2 + (sin(theta)/sigy)^2);
b = 1/4 * (-sin(2*theta)/sigx^2 + sin(2*theta)/sigy^2);
c = 1/2 * ((sin(theta)/sigx)^2 + (cos(theta)/sigy)^2);

xdiff = xy(:, 1) - x0;
ydiff = xy(:, 2) - y0;

z = A*exp(-(a*xdiff.^2+2*b*xdiff.*ydiff + c*ydiff.^2));
end

function z = gauss2d(coef, xy)
z = 0;
l = length(coef);
if mod(l, 6) == 0 || mod(l, 6) == 1
    num = l / 6;
    for i = 1:num
        ind = ((i-1)*6+1):(i*6);
        z = z + gauss2d0(coef(ind), xy);
    end
    if mod(l, 6) == 1
        z = z + coef(end); % offset term
    end
end
%disp(size(z))
end

function z = gauss2d_pal(coef, xy)
l = length(coef);
if mod(l, 6) == 0
    coef_reshape = reshape(coef',[6, l / 6]);
    coef_reshape =  coef_reshape';
    num = l / 6;
end
if mod(l, 6) == 1
    coef_mod = coef(1:end-1);
    coef_reshape = reshape(coef_mod',[6,(l-1) / 6]);
    coef_reshape = coef_reshape';
    num = (l-1)/6;
end
    
    sz= size(xy);
    z_indiv = zeros(num,max(sz));
parfor ind = 1:num
       coef_indiv = coef_reshape(ind,:);
       z_indiv(ind,:) = gauss2d0(coef_indiv, xy);
       disp(ind);
end
z=sum(z_indiv);
if mod(l, 6) == 1
    z = z + coef(end); % offset term
end

z=z';
end