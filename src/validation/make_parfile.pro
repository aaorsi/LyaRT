; Routine to create Parameters file,
; where most of the parameters are given in the input (for using galform results)

pro make_parfile,	GEOMNAME = geomname,TEMP = temp,MEAN_NH = mean_nh, V = v,$
					NPHOTONS=nphotons,TAU0=tau0,XCRIT=xcrit,IncDust = incdust,$
					NCells = ncells,FTAG=ftag,PARDIR=pardir,OUTDIR=outdir,$
					Set_tolerance=set_tolerance,nout_max = nabs_max, np_max = np_max,np_min=np_min,$
					redshift = redshift,IncShield=incshield,IncUV=IncUV,Z = Z,$
					Rout = rout, Rinn= rinn,b = b

	;T must be given in log units.

	OutMode = 'Short'
	
	NCellX 	 = '1'
	NCellY   = '1'
	NCellZ   = '1'

	T4		 = Temp/10000.

	If not keyword_set(Set_Tolerance) then  $
		Set_tolerance = 'no'
	If not keyword_set(nout_max) then $
		nabs_max = 100
	If not keyword_set(np_max) then $
		np_max = NPhotons
	If not keyword_set(np_min) then $
		np_min = 100


	sigma_const = 5.8689e-14 
	ID = FTag

	spawn,'mkdir -p '+ParDir
	spawn,'mkdir -p '+OutDir

	OutName  = OutDir + FTag

	gtype = GeomName eq 'HomSlab' ? 0 : 1

	if gtype eq 0 then begin

		zSize = tau0 * sigma_const^(-1)*sqrt(T4)*mean_NH^(-1)*2
		ySize = zSize
		xSize = zSize
		RSphere = zSize/2.
	endif


	if GeomName eq 'HomSphere' then begin

		RSphere = tau0 * sigma_const^(-1)*sqrt(T4)*mean_NH^(-1)
		zSize = RSphere*2.
		xSize = RSphere*2.
		ySize = RSphere*2.
	endif

	if GeomName eq 'Shell_VConst' or GeomName eq 'HomShell' then begin
		
		xSize = 2*ROut
		ySize = 2*ROut
		zSize = 2*ROut
		Rsphere = rout
	endif

	if Geomname eq 'ThinShell' then begin
		
		xSize	= 2 * ROut
		ySize	= 2 * Rout
		zSize	= 2 * Rout	
		RSphere = rout
	endif



	ParName = Ftag+'.par'
				 	
	openw,1,ParDir+ParName

		printf,1,'GeomName		= '+GeomName
		printf,1,'Set_Tolerance = '+Set_Tolerance
		printf,1,'nout_max	= ',nabs_max
		printf,1,'np_max	= ',long(np_max)
		printf,1,'np_min	= ',long(np_min)

		if keyword_set(redshift) then $
		printf,1,'Redshift	= ',redshift
		if keyword_set(IncUV) then $
		printf,1,'IncUV		= ',IncUV 
		if keyword_set(IncShield) then $
		printf,1,'IncShield	= ',fix(IncShield)
		if keyword_set(IncDust) then $
		printf,1,'IncDust  = ',IncDust
		printf,1,'NPhotons  = ',NPhotons
		printf,1,'NCellX	= ',NCellX
		printf,1,'NCellY	= ',NCellY
		printf,1,'NCellZ	= ',NCellZ
		printf,1,'NCells	= ',NCells
		printf,1,'XPeriodic = '+'1'
		printf,1,'YPeriodic = '+'1'
		printf,1,'ZPeriodic = '+'0'
		printf,1,'xSize		= ',double(xSize),format='((A),(G))'
		printf,1,'ySize		= ',double(ySize),format='((A),(G))'
		printf,1,'zSize		= ',double(zSize),format='((A),(G))'
		if keyword_set(Rinn) then $
		printf,1,'R_inner	= ',double(Rinn),format='((A),(G))'

		printf,1,'RSphere	= ',double(RSphere),format='((A),(G))'
		
		printf,1,'mean_nH	= ',mean_nh
		printf,1,'Temp		= ',alog10(T4)
		printf,1,'vmax		= ',V
		if keyword_set(Z) then $
		printf,1,'Z			= ',Z
	
		if keyword_set(b) then $
		printf,1,'b			= ',b

		printf,1,'xcritval  = ',xcrit
		printf,1,'OutShort	= '+OutName
		printf,1,'OutMode	= '+OutMode
		
		close,1

end
	


	
