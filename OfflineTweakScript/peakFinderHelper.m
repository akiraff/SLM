function out = peakFinderHelper(trapInfo_mod, repeat, spacing, sortingSpacing)
% Goal of this function is to get rid of closely spaced peak points
    for re = 1:repeat
        trap_nonzero = nonzeros(trapInfo_mod);
        peak_sort_mod = sortrows(trapInfo_mod,'descend');
        trapInfo = peak_sort_mod;
        % sort by y, remove closely spaced point
        trapInfo_pos = sortrows(trapInfo,3,'descend');
        trapInfo_mod =  trapInfo_pos;
        xpeak_sites = trapInfo_mod(:,2);
        ypeak_sites = trapInfo_mod(:,3);
        peak_sites = trapInfo_mod(:,1);
        for j = 1:length(trap_nonzero)/3-sortingSpacing
            x_cur =  xpeak_sites(j);
            y_cur =  ypeak_sites(j);
            peak_cur = peak_sites(j);
            x_next = xpeak_sites(j+sortingSpacing);
            y_next = ypeak_sites(j+sortingSpacing);
            peak_next =  peak_sites(j+sortingSpacing);
            d = ((x_cur-x_next)^2 + (y_cur-y_next)^2)^0.5;
            if d < spacing/4
                if peak_cur < peak_next
                    trapInfo_mod(j,:) = 0;
                else
                    trapInfo_mod(j+sortingSpacing,:) = 0;
                end
                disp(d);
                disp("find closely spaced point");
            end
        end
        trap_nonzero = nonzeros(trapInfo_mod);
        trapInfo_mod = reshape(trap_nonzero, [length(trap_nonzero)/3,3]);
        disp(size(trapInfo_mod));
        peak_sort_mod = sortrows(trapInfo_mod,'descend');
        trapInfo = peak_sort_mod;
        % sort by x, remove closely spaced point
        trapInfo_pos = sortrows(trapInfo,2,'descend');
        trapInfo_mod =  trapInfo_pos;
        xpeak_sites = trapInfo_mod(:,2);
        ypeak_sites = trapInfo_mod(:,3);
        for j = 1:length(trap_nonzero)/3-sortingSpacing
            x_cur =  xpeak_sites(j);
            y_cur =  ypeak_sites(j);
            x_next = xpeak_sites(j+sortingSpacing);
            y_next = ypeak_sites(j+sortingSpacing);
            d = ((x_cur-x_next)^2 + (y_cur-y_next)^2)^0.5;
           % disp(d);
            if d < spacing/4
                trapInfo_mod(j,:) = 0;
                 disp("Find closely spaced spots!");
                 disp(d);
            end
        end
        
        trap_nonzero = nonzeros(trapInfo_mod);
        trapInfo_mod = reshape(trap_nonzero, [length(trap_nonzero)/3,3]);
        disp(size(trapInfo_mod));
    end
    out = trapInfo_mod;
end

