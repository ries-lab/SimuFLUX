function out=convnfft(a,b)
afft=fftn(a);
bfft=fftn(b,size(a));
cfft=afft.*bfft;
out=ifftn(cfft);