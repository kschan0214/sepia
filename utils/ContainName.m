% return boolean value to check if the input name contains certain string
function bool = ContainName(name,string)
    bool= ~isempty(strfind(lower(name),string));
end