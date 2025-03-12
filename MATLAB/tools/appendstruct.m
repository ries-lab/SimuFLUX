function sout=appendstruct(sin, sadd)
if isempty(sin)
    sout=sadd;
    return
end
sout=sin;
fn=fieldnames(sin);
ln=length(sadd.(fn{1}));
for k=1:length(fn)
    if isfield(sadd,fn{k})
     sout.(fn{k})(end+1:end+ln)=sadd.(fn{k});
    end
end
end