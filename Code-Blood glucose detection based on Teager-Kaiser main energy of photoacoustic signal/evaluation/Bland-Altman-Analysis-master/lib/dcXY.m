function str = dcXY(data,x,y)
str = {};
if ~isempty(data.xName), str = [str, [data.xName ': ' num2str(x)]]; end
if ~isempty(data.yName), str = [str, [data.yName ': ' num2str(y)]]; end
end