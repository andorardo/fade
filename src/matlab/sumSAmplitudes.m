function res = sumSAmplitudes(s,threshold,smooth)
    if nargin<3, smooth = 0; end;
    smooths = s;
    if smooth>0
        smooths = gaussianBlur(s,smooth);
    end
    slen = length(smooths);
    ils = findIslands(smooths,threshold);
    res = [];
    if isempty(ils), return; end
    spans = ils;
    for i=1:size(spans,1)
        x1 = spans(i,1);
        x2 = spans(i,2);
        while x1>1 && smooths(x1-1)<smooths(x1) && smooths(x1-1)>0
            x1 = x1-1;
        end
        while x2<slen && smooths(x2+1)<smooths(x2) && smooths(x2+1)>0
            x2 = x2+1;
        end
        spans(i,:) = [x1,x2];
    end
    collision = diff(reshape(spans',[],1))==0;
    collision = find(collision(2:2:end));
    smooths(spans(collision,2)) = 0.5*smooths(spans(collision,2));
    res = arrayfun( @(s1,s2) sum(smooths(s1:s2)), spans(:,1), spans(:,2) );
end
