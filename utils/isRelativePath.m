function isRelative = isRelativePath(pathStr)
    % Check if the path is absolute
    if ispc  % Windows
        isRelative = isempty(regexp(pathStr, '^[A-Za-z]:\\', 'once'));
    else     % Unix/Mac
        isRelative = pathStr(1) ~= '/';
    end
end
