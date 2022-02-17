function trapGood = WrongSitesRemover(trap, wrong_peak)
% This function is used to remove the wrong sites 
% trap is [x,y,z] for peak location
% wrong_peak: peak value of the wrong sites
num_wrong = length(wrong_peak);
disp(num_wrong);
for j = 1:num_wrong
    [row,~]=find(trap == wrong_peak(j)); 
    trap(row,:)=[];
    disp(row);
end
     % check after removing
     peak_sites = trap(:,1);
     xpeak_sites = trap(:,2);
     ypeak_sites = trap(:,3);
     figure(1000);clf;
      plotc(xpeak_sites,ypeak_sites,peak_sites, 'ro');
     % output
     trapGood = trap;
end

