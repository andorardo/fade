function islands = findIslands(s, threshold)
    threshold = min(threshold, max(s)/2);
    mask = (s - threshold)>0;
    d = (diff(mask));
    boundaries = find(abs(d)>0);
    islands = [];
    if isempty(boundaries), return; end
    if d(boundaries(1))<0, boundaries = boundaries(2:end); end
    if d(boundaries(end))>0, boundaries = boundaries(1:end-1); end
    islands = reshape(boundaries(:),2,[])'; 
    islands(:,1) = islands(:,1) + 1;
end
