%adapted from Austin's Aerosol Lib 7/27/20

function mfp = getMFP(temp, P)
	
		%based on Willeke's (1976) Temperature Dependence of Particle Slip in a Gaseous Medium. Relations valid from 200k to 1000k
		%https://doi.org/10.1016/0021-8502(76)90024-0

		if(or(temp > 1000,temp < 200))
			print("CRITICAL ERROR: TEMP OUT OF BOUNDS FOR THIS MFP CALCULATION Willeke(1976)/n");
        end
		mfp0 = 67.3e-9;
		T0 = 296.15;
		p0 = 101325.0;

		mfp =  mfp0.* (temp ./ T0).* (p0 ./ P).* (1.0 + 110.4 ./ T0) ./ (1.0 + 110.4 ./ temp);
        
        end