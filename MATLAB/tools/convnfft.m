function out=convnfft(a,b)
afft=fftn(a);
bfft=fftn(b,size(a));
cfft=afft.*bfft;
out=ifftn(cfft);
for k=1:length(size(a)) %dims
    out=ifftshift(out,k);
end

