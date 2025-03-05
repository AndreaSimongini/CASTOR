The structure of this directory has to be: 

- sn_name (i.e. SN1987A) 

	-- light_curves 
	
		--- sn_name.dat (i.e. SN1987A.dat) 
		
			----- time [MJD], mag [AB], emag [AB], filter 
	
	-- spectra (if available) 
	
		--- sn_name_epoch.dat (i.e. SN1987A_46850.68.dat) 
		
			----- wave [AA], flux [erg/s/cm2]
			
			
			
			
NB: you MUST respect this scheme and nomenclature, otherwise CASTOR will not read properly your data. 
