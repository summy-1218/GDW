% ===== 子函数 =====
function val = getField(s, name, defaultVal)
if isfield(s, name)
    val = s.(name);
else
    val = defaultVal;
end
end