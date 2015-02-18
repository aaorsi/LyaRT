; to produce a set of parameter files to be used
; by LyaRT to produce results used later to validate
; the code.


pro make_validation_pars

;	RootDir = '/home/jmejia/LyaRT/'
	RootDir ='/home/aaorsi/work/LyaRT/' 
	PDir = RooTDir + '/data/Params/validation/'
;	ODir = RooTDir + '/out/validation/'
	ODir = '/home/aaorsi/work/LyaRT/out/validation/'

	Test_HomSlab	= 1				; Homogeneous Slab
	Test_HomSphere	= 1				; Homogeneous Sphere
	Test_NScat		= 1				; Number of scatterings from an homogeneous slab
	Test_fesc_slab	= 1				; f_esc from a dusty homogeneous slab
	Test_ExpSphere	= 1				; Expanding sphere
	Test_ThinShell  = 1				; An homogeneous thin shell



	if Test_HomSlab then begin
		print,'Testing HomSlab'	
		GeomName 	 = 'HomSlab'
		log_tau0Arr	 = [5.0,6.0,7.0]
		Temp		 = 10.		; [K]
		xcritArr	 = [0.,3.]
		nx			 = n_elements(xcritArr)
		
		for ix = 0,nx-1 do begin
			xcrit = xcritArr[ix]
			ntau		 = n_elements(log_Tau0Arr)		
			Name		 = GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		
			NPhotons	 = 10000
			NParFiles	 = 100
			mean_nh		 = 100.
			spawn,'mkdir -p '+PDir + Name
			ParList		 = PDir + Name+'/file_list'
			openw,2,ParList
		
			for i = 0,ntau-1 do begin			
				ltau = log_Tau0Arr[i]
				t0 = 10^ltau
				for j = 0,nParFiles-1 do begin
		
					FTag = 'log_tau'+strn(ltau,len=3)+'.'+strn(j)
					print,FTag
					make_parfile,GeomName=GeomName,Temp = Temp,Mean_Nh = mean_nh,$
					V = 0.,NPhotons = NPhotons,Tau0 = t0, xcrit = xcrit, NCells = 1,$
					FTag = FTag,ParDir = PDIr + Name+'/',OutDir=ODir + Name+'/',$
					IncDust = 0
					
					printf,2,FTag+'.par'
	
				endfor
			endfor

			close,2
		endfor		

	endif

	if Test_HomSphere then begin
		print,'Testing HomSphere'
		GeomName 	 = 'HomSphere'
		log_tau0Arr	 = [5.,6.,7.]
		Temp		 = 10.		; [K]
		xcrit		 = 3.


		ntau		 = n_elements(log_Tau0Arr)		
		Name		 = GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		
		NPhotons	 = 10000
		NParFiles	 = 100
		mean_nh		 = 100.

		spawn,'mkdir -p '+PDir + Name
		ParList		 = PDir + Name+'/file_list'
		openw,2,ParList
		
		for i = 0,ntau-1 do begin			
			ltau = log_Tau0Arr[i]
			t0 = 10^ltau
			for j = 0,nParFiles-1 do begin
		
				FTag = 'log_tau'+strn(ltau,len=3)+'.'+strn(j)
				print,FTag
				make_parfile,GeomName=GeomName,Temp = Temp,Mean_Nh = mean_nh,$
				V = 0.,NPhotons = NPhotons,Tau0 = t0, xcrit = xcrit, NCells = 1,$
				FTag = FTag,ParDir = PDIr + Name+'/',OutDir=ODir + Name+'/',$
				IncDust = 0
				
				printf,2,FTag+'.par'

			endfor
		endfor

		close,2

	endif


	if Test_NScat then begin
		
		print,'Testing NScat'
		GeomName 	 = 'HomSlab'
		log_tau0Arr	 = [4.0,5.0,6.0,7.0,8.0,9.0]
		Temp		 = 10.		; [K]
		xcrit		 = 0.


		ntau		 = n_elements(log_Tau0Arr)		
		Name		 = 'NScat'+GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		
		NPhotons	 = 50
		NParFiles	 = [40,40,40,40,50,100]
		mean_nh		 = 100.
		spawn,'mkdir -p '+PDir + Name
		ParList		 = PDir + Name+'/file_list'
		openw,2,ParList
		
		for i = 0,ntau-1 do begin			
			ltau = log_Tau0Arr[i]
			t0 = 10^ltau
			for j = 0,nParFiles[i]-1 do begin
		
				FTag = 'log_tau'+strn(ltau,len=3)+'.'+strn(j)
				print,FTag
				make_parfile,GeomName=GeomName,Temp = Temp,Mean_Nh = mean_nh,$
				V = 0.,NPhotons = NPhotons,Tau0 = t0, xcrit = xcrit, NCells = 1,$
				FTag = FTag,ParDir = PDIr + Name+'/',OutDir=ODir + Name+'/',$
				IncDust = 0
				
				printf,2,FTag+'.par'

			endfor
		endfor

		close,2

	endif

	if Test_fesc_slab then begin
		print,'Testing fesc from HomSlab'	

		Albedo		 = 0.39
		z_star		 = 0.02
		ext_zstar	 = 8.87e-20
		sigma_nh	 = 5.8689e-14

		GeomName 	 = 'HomSlab'
		log_tau0Arr	 = 6.
		Temp		 = 10.		; [K]
		xcrit		 = 3.
		T4			 = Temp * 1e-4
		a_par		 = 4.693e-4 * sqrt(1./T4)

		minx		 = -2.0
		maxx		 = 1.0
		binx		 = 0.5
		nx			 = (maxx - minx)/binx + 1.0
				
		xaxis 		 = findgen(nx)*binx + minx

		ntau		 = n_elements(log_Tau0Arr)		
		
		Name		 = 'fesc_'+GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		
		NPhotons	 = 1000
		NParFiles	 = 10
		mean_nh		 = 100.
		spawn,'mkdir -p '+PDir + Name
		ParList		 = PDir + Name+'/file_list'
		openw,2,ParList
		
		for i = 0,ntau-1 do begin			
			ltau = log_Tau0Arr[i]
			t0 = 10^ltau
			for iz = 0,nx-1 do begin
				Z = (10d^(xaxis[iz]) * sigma_nH * sqrt(1./T4))/((1.-Albedo) * ext_zstar * t0 * (a_par*t0)^(1./3)) 
				print,'Z = ',Z
				for j = 0,nParFiles-1 do begin
		
					FTag = 'xvar_'+strn(xaxis[iz],len=4)+'log_tau'+strn(ltau,len=3)+'.'+strn(j)
					print,FTag
					make_parfile,GeomName=GeomName,Temp = Temp,Mean_Nh = mean_nh,$
					V = 0.,NPhotons = NPhotons,Tau0 = t0, xcrit = xcrit, NCells = 1,$
					FTag = FTag,ParDir = PDIr + Name+'/',OutDir=ODir + Name+'/',$
					Z = Z, IncDust = 'Yes'
					printf,2,FTag+'.par'
				endfor
			endfor
		endfor

		close,2

	endif

;	#########

	if Test_ExpSphere then begin
		print,'Testing Expanding Sphere'
		GeomName 	 = 'HomSphere'
		log_tau0Arr	 = 7.06
		Temp		 = 10000.		; [K]
		xcrit		 = 3.
		vmaxarr		 = [0.,20.,200.,2000]
		
		nv			 = n_elements(vmaxarr)	
		Name		 = 'Expanding_'+GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		
		NPhotons	 = 2000
		NParFiles	 = 50 
		mean_nh		 = 100.

		NCells		 = 1000

		spawn,'mkdir -p '+PDir + Name
		ParList		 = PDir + Name+'/file_list'
		openw,2,ParList
		ltau = log_Tau0Arr
		t0 = 10^ltau
		for i = 0,nv-1 do begin			
			vmax	= vmaxarr[i]			
			for j = 0,nParFiles-1 do begin
		
				FTag = 'vmax_'+strn(vmax,len=4)+'log_tau'+strn(ltau,len=3)+'.'+strn(j)
				print,FTag
				make_parfile,GeomName=GeomName,Temp = Temp,Mean_Nh = mean_nh,$
				V = vmax,NPhotons = NPhotons,Tau0 = t0, xcrit = xcrit, NCells = NCells,$
				FTag = FTag,ParDir = PDIr + Name+'/',OutDir=ODir + Name+'/',$
				IncDust = 0
				
				printf,2,FTag+'.par'

			endfor
		endfor

		close,2

	endif

	if Test_ThinShell then begin
		print,'Testing thin shell'
		GeomName 	 = 'ThinShell'
		b			 = 40.0		; kms-1
		column_d	 = 2.0e20
		tau0		 = 3.8e6
		xcrit		 = 0.
		vmax		 = 300.
		RMax		 = 1.17e19
		f			 = 0.9			; RMin/RMax
		RMin		 = RMax * f			
		
		T4			 = (b/12.85)^2
;		T4		     = 1.0
		Temp		 = T4 * 1.e4

		logTemp		 = round(alog10(T4)*100.)/100.
		
		print,'log Temp/1e4K = ',logTemp
	
		Name		 = +GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		sigma_const = 5.8689e-14 
		
		tau2		= sigma_const * (12.85/b) * column_d

		mean_nh		= (1./(RMax*(1 - f))) * column_d

		NPhotons	 = 2000
		NParFiles	 = 30

		NCells		 = 1000
	
		spawn,'mkdir -p '+PDir + Name
		ParList		 = PDir + Name+'/file_list'
		openw,2,ParList
		print,ParList
		ltau = alog10(tau0)
		for j = 0,nParFiles-1 do begin
		
				FTag = 'vmax_'+strn(vmax,len=4)+'log_tau'+strn(ltau,len=3)+'.'+strn(j)
				print,FTag
				make_parfile,GeomName=GeomName,Temp = Temp,Mean_Nh = mean_nh,$
				V = vmax,NPhotons = NPhotons,Tau0 = t0, xcrit = xcrit, NCells = NCells,$
				FTag = FTag,ParDir = PDIr + Name+'/',OutDir=ODir + Name+'/',$
				IncDust = 0,Rout = RMax,Rinn = RMin,b = b
				
				printf,2,FTag+'.par'

		endfor

		close,2

	endif









end		
			

				
