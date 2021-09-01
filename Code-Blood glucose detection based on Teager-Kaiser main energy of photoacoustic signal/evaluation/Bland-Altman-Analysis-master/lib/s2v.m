function s2v(s)
%#ok<*NODEF>
fn = fieldnames(s);
for f = fn(:).'
    assignin('caller',f{:},s.(f{:}))
end
end