function res = peakSAmplitudes(s,threshold,smooth)
    if nargin<3, smooth = 0; end;
    smooths = s;
    if smooth>0
        smooths = gaussianBlur(s,smooth);
    end
    ils = findIslands(smooths,threshold);
    res = [];
    if isempty(ils), return; end
    res = ils(:,1); %+ 1;
end

