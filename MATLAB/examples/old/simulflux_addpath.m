function simulflux_addpath
dirs={'PSF_simulation/library', 'tools','settings','examples'};

fh=fileparts(mfilename('fullpath')); %add all folders to serach path
if contains(fh,"examples")
    fh=fileparts(fh);
end
for k=1:length(dirs)
    addpath([fh,filesep,dirs{k}])
end

end