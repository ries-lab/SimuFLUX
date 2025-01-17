function sout=sumstruct(sin, sadd)
if isempty(sin)
    sout=sadd;
    return
end
if isempty(sadd)
    sout=sin;
    return
end
fn=fieldnames(sin);
for k=1:length(fn)
    sout.(fn{k})=sin.(fn{k})+sadd.(fn{k});
end
end