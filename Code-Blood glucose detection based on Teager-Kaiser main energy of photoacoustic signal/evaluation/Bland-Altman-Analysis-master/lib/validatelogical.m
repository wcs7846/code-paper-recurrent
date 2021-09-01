function tf = validatelogical(L)
%VALIDATELOGICAL True for arrays that can be converted to logical.
%   VALIDATELOGICAL(L) returns logical 1 (true) if L is a logical or
%   numeric array and logical 0 (false) otherwise.
%
%   See also LOGICAL.

try tf = islogical(logical(L)); catch, tf = false; end
end