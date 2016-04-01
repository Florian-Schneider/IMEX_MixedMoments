function b = isIsotropic(f)
s = functions(f);

if strcmp(s.function(3),'x')
    b = 1;
else
    b = 0;
end

end