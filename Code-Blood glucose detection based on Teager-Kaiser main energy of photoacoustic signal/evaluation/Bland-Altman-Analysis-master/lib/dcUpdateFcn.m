function str = dcUpdateFcn(~,e)
x = e.Position(1);
y = e.Position(2);
data = e.Target.UserData;
str = data.dcFun(data,x,y);
if ~isempty(data.prefix), str = [data.prefix, str]; end
if ~isempty(data.suffix), str = [str, data.suffix]; end
end