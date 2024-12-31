function [v,s] = getOr(s, field, default)
%getOr Returns the structure field or a default if either don't exist
%   v = getOr(s, field, [default]) returns the 'field' of the structure 's'
%   or 'default' if the structure is empty of the field does not exist. If
%   default is not specified it defaults to []. 'field' can also be a cell
%   array of fields, in which case it will check for all of them and return
%   the value of the first field that exists, if any (otherwise the default
%   value).
% this was copied from kilosort 2.0 matlab package

if nargin < 3
  default = [];
end

fieldExists = isfield(s, field);
if any(fieldExists) && ~isempty(s.(field))
  if iscellstr(field)
    v = s.(field{find(fieldExists, 1)});
  else
    v = s.(field);
  end
else
  v = default;
end


if nargout > 1
  s.(field) = v;
end

end

