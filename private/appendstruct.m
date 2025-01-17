function sout=appendstruct(sin, sadd)
sout=sin;
fn=fieldnames(sin);
ln=length(sadd.(fn{1}));
for k=1:length(fn)
    sout.(fn{k})(end+1:end+ln)=sadd.(fn{k});
end
end