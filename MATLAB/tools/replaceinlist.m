function sin=replaceinlist(sin,varargin)
sstr=varargin(1:2:end);
rstr=varargin(2:2:end);
for l=1:length(sin)
    for k=1:length(sstr)
        if strcmp(sin{l},sstr{k})
            sin{l}=rstr{k};
        end
    end
end
end