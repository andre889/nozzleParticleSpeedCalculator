%Adapted from Austin's Aerosol Lib 7/27/20
function Cc = getCc(dp, mfp)
		
%from Kim(2005) PSL with uncertainty if interested. http://dx.doi.org/10.6028/jres.110.005.
a = 1.165;
b = 0.486;
c = -0.997;
		kn = 2.0 * mfp ./ dp;
		Cc = 1.0 + kn .* (a + b .* exp(c ./ kn));	
end
