function trapGood = RightSitesAdder(trap, right_trap)
% This function is to add right sites 
  [row, ~]= size(right_trap);
  disp(row);
  for i = 1:row
      trap(end+1,:) = right_trap(i,:);
      disp(trap);
  end
      % check after adding
     peak_sites = trap(:,1);
     xpeak_sites = trap(:,2);
     ypeak_sites = trap(:,3);
     figure(2000);clf;
      plotc(xpeak_sites,ypeak_sites,peak_sites, 'ro');
      trapGood = trap;
end

