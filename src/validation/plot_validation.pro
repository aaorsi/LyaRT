; This script reads and plots data to be tested 
; against analytical formulae or other tests.

Function J_Neufeld,x,a_tau0
	
	J = (sqrt(6/!pi)/24.) * (x^2/a_tau0) * $
	(1./(cosh(sqrt((!pi^3)/54.)*(abs(x^3))/(a_tau0)))) 
	
	return,J
end	

Function J_Dijkstra, x,a_tau0
	J = double((sqrt(!pi/24.) * (x^2/a_tau0)) * $
	(1./(1 + cosh(sqrt((2*!pi^3)/27.)*(abs(x^3))/(a_tau0)))) )
	
	return,J
end	

Function NScat_Harrington,tau0
	return,1.612*tau0
end

Function fesc_neufeld,xvar
	zeta = 0.525
	zeta_ = sqrt(3)/(zeta * !PI^(5./12))
	return,1./cosh(zeta_ * sqrt(xvar))
end


pro plot_validation

	NMAX			= 1e6
	PlotPS			= 1
		
	Plot_HomSlab	= 1
	Plot_HomSphere	= 1 
	Plot_NScat		= 0
	Plot_fesc		= 0
	Plot_ExpSphere	= 0
	Plot_ThinShell  = 0

	RootDir ='/home/aaorsi/LyaRT/' 
	PDir = RooTDir + '/data/Params/validation/'
	ODir = RooTDir + '/out/validation/'
	PlotDir = 'plots/'
	
	wc	= 0l

	if Plot_HomSlab then begin
		PlotFile	 = 'J_HomSlab'
		GeomName 	 = 'HomSlab'
		log_tau0Arr	 = [5.,6.,7.]
		Temp		 = 10.		; [K]
		xcritarr	 = [0.,3.]
		nx			 = n_elements(xcritarr)

		ntau		 = n_elements(log_Tau0Arr)		
		
		NPhotons	 = 10000
		NParFiles	 = 100
		mean_nh		 = 100.


		xr = [-120,120]
		yr = [0,0.0050]


		if PlotPs then begin
			set_plot,'ps'
			!P.Font = 0
			device,/color,/landscape,/palatino,filename	=PlotDir + PlotFile+'.ps'
		endif else begin
			window,wc,xsiz=850,ysiz=650,retain=2
			wc++
		endelse

		defplotcolors
		col = [!orange,!blue,!red,!dgreen,!purple]
		del = '\Delta'

		plot,indgen(2),/nodata,xr=xr,yr=yr,/xs,/ys,$
		xtitle=textoidl('x = (\nu - \nu_0)/'+del+'\nu_D'),$
		ytitle = textoidl('J(\tau_0,x)'),$
		charsiz=2,position=aspect(.8),xthic=4,ythic=4
	
		for ix = 0,nx-1 do begin
			xcrit 		 = xcritarr[ix]
			Name		 = GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
			OutputDir	 = ODir + Name+'/'

			for i = 0,ntau-1 do begin
				ltau = log_Tau0Arr[i]
				t0 = 10^ltau

				xarr 		= fltarr(NMAX)
				xarr_abs 	= fltarr(NMAX)
				nscat	 	= fltarr(NMAX)
				nout		= 0l
				np			= 0l

				for j = 0,nParFiles-1 do begin
					FTag = 'log_tau'+strn(ltau,len=3)+'.'+strn(j)
					OutFile = OutputDir + FTag
				
					openr,1,OutFile,err=err

					if err ne 0 then begin
						print,'File '+OutFile+' not found. Skipping it...'
						close,1
						continue
					endif
					print,OutFile
					str=''
					readf,1,str

					while ~EOF(1) do begin
						readf,1,str
						res = strsplit(str,/ext)
						nel = n_elements(res)	

						if res[0] eq '#' then continue

						if res[nel-1] eq 'Photon_escaped' then begin
							xarr[nout]  = res[1]		
							nscat[nout] = res[0]
							nout++	
						endif
						np++
					endwhile
					close,1
				endfor

				print,'Nphotons computed =',np

				if np lt 10 then begin
					print,'Only ',np,' photons found. Forget it...'
					continue
				endif

				xarr = xarr[0:nout-1]
				nscat = nscat[0:nout-1]

				if np lt 100 then databin = 5.0
				if np gt 100 and np lt 1000 then databin = 2.5
				if np gt 1000 and np lt 5e3 then databin = 2.0
				if np gt 5e3 and np lt 1e4 then databin = 1.5
				if np gt 1e4 and np lt 5e4 then databin = 1.0
				if np ge 5e4 then databin = .5
				
				T4 = Temp*1e-4
				a = 4.693e-4 * T4^(-1./2)
				a_tau0 = a*t0

				xnmax = 120
				xnmin = -120
				xnbin = 0.1
				Nxneu = (xnmax - xnmin)/xnbin + 1

				xneu = findgen(Nxneu)*xnbin + xnmin
				
				jneu = J_Neufeld(xneu,a_tau0)
				Norm_Neu = int_tabulated(xneu,jneu,/double)

				jneu = jneu/norm_Neu*(1./(4*!PI))
				mspec = max(jneu)
				dr = 1

				plothist,xarr,xdata,ydata,bin=databin,/noplot
				ydata = ((ydata+0.)/(np+0.)/databin)*(1./(4*!PI))

				oplot,xneu,jneu,thic=3
				oplot,xdata,ydata,thic=4 - ix*2,color=col[ix],ps=10

			endfor
		endfor
	
		legend,textoidl(['x_{crit} = '+strn(fix(xcritarr[0])),$
		 				 'x_{crit} = '+strn(fix(xcritarr[1]))]),$
		/right,/top,textcol=col[0:nx-1],charsiz=2,box=0

;		xyouts,10,0.0090,textoidl('\tau_0 = 10^4'),charsiz=1.5
		xyouts,15,0.0045,textoidl('\tau_0 = 10^5'),charsiz=1.5
		xyouts,25,0.0025,textoidl('\tau_0 = 10^6'),charsiz=1.5
		xyouts,60,0.0015,textoidl('\tau_0 = 10^7'),charsiz=1.5


		if PlotPs then device,/close
	endif

;########


	if Plot_HomSphere then begin
		PlotFile	 = 'J_HomSphere'
		GeomName 	 = 'HomSphere'
		log_tau0Arr	 = [5.,6.,7.]
		Temp		 = 10.		; [K]
		xcritarr	 = [3.]
		nx			 = n_elements(xcritarr)

		ntau		 = n_elements(log_Tau0Arr)		
		
		NPhotons	 = 10000
		NParFiles	 = 100
		mean_nh		 = 100.


		xr = [-120,120]
		yr = [0,0.0060]


		if PlotPs then begin
			set_plot,'ps'
			!P.Font = 0
			device,/color,/landscape,/palatino,filename	=PlotDir + PlotFile+'.ps'
		endif else begin
			window,wc,xsiz=850,ysiz=650,retain=2
			wc++
		endelse

		defplotcolors
		col = [!orange,!orange,!orange]
		del = '\Delta'

		plot,indgen(2),/nodata,xr=xr,yr=yr,/xs,/ys,$
		xtitle=textoidl('x = (\nu - \nu_0)/'+del+'\nu_D'),$
		ytitle = textoidl('J(\tau_0,x)'),$
		charsiz=2,position=aspect(.8),xthic=4,ythic=4
	
		for ix = 0,nx-1 do begin
			xcrit 		 = xcritarr[ix]
			Name		 = GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
			OutputDir	 = ODir + Name+'/'

			for i = 0,ntau-1 do begin
				ltau = log_Tau0Arr[i]
				t0 = 10^ltau

				xarr 		= fltarr(NMAX)
				xarr_abs 	= fltarr(NMAX)
				nscat	 	= fltarr(NMAX)
				nout		= 0l
				np			= 0l

				for j = 0,nParFiles-1 do begin
					FTag = 'log_tau'+strn(ltau,len=3)+'.'+strn(j)
					OutFile = OutputDir + FTag
				
					openr,1,OutFile,err=err

					if err ne 0 then begin
						print,'File '+OutFile+' not found. Skipping it...'
						close,1
						continue
					endif
					print,OutFile
					str=''
					readf,1,str

					while ~EOF(1) do begin
						readf,1,str
						res = strsplit(str,/ext)
						nel = n_elements(res)	

						if res[0] eq '#' then continue

						if res[nel-1] eq 'Photon_escaped' then begin
							xarr[nout]  = res[1]		
							nscat[nout] = res[0]
							nout++	
						endif
						np++
					endwhile
					close,1
				endfor

				print,'Nphotons computed =',np

				if np lt 10 then begin
					print,'Only ',np,' photons found. Forget it...'
					continue
				endif

				xarr = xarr[0:nout-1]
				nscat = nscat[0:nout-1]

				if np lt 100 then databin = 5.0
				if np gt 100 and np lt 1000 then databin = 2.5
				if np gt 1000 and np lt 5e3 then databin = 2.0
				if np gt 5e3 and np lt 1e4 then databin = 1.5
				if np gt 1e4 and np lt 5e4 then databin = 1.0
				if np ge 5e4 then databin = .5
				
				T4 = Temp*1e-4
				a = 4.693e-4 * T4^(-1./2)
				a_tau0 = a*t0

				xnmax = 120
				xnmin = -120
				xnbin = 0.1
				Nxneu = (xnmax - xnmin)/xnbin + 1

				xneu = findgen(Nxneu)*xnbin + xnmin
				
				jneu = J_Dijkstra(xneu,a_tau0)
				Norm_Neu = int_tabulated(xneu,jneu,/double)

				jneu = jneu/norm_Neu*(1./(4*!PI))
				mspec = max(jneu)
				dr = 1

				plothist,xarr,xdata,ydata,bin=databin,/noplot
				ydata = ((ydata+0.)/(np+0.)/databin)*(1./(4*!PI))

				oplot,xneu,jneu,thic=3
				oplot,xdata,ydata,thic=4,color=col[i],ps=10

			endfor
		endfor

		xyouts,10,0.0055,textoidl('\tau_0 = 10^5'),charsiz=1.5
		xyouts,25,0.0030,textoidl('\tau_0 = 10^6'),charsiz=1.5
		xyouts,50,0.0015,textoidl('\tau_0 = 10^7'),charsiz=1.5

		if PlotPs then device,/close
	endif


;#######

	if Plot_NScat then begin
		PlotFile	 = 'NScat'
		GeomName 	 = 'HomSlab'
		log_tau0Arr	 = [4.,5.,6.,7.,8.,9.]
		Temp		 = 10.		; [K]
		xcrit		 = 0.
		nx			 = n_elements(xcritarr)

		ntau		 = n_elements(log_Tau0Arr)		
		
		NPhotons	 = 10000
		NParFiles	 = [10,20,30,40,50,100]
		mean_nh		 = 100.


		tau0arr = 10^(log_tau0Arr)


		if PlotPs then begin
			set_plot,'ps'
			!P.Font = 0
			device,/color,/landscape,/palatino,filename	=PlotDir + PlotFile+'.ps'
		endif else begin
			window,wc,xsiz=650,ysiz=650,retain=2
			wc++
		endelse

		defplotcolors
		col = [!orange,!blue,!red,!dgreen,!purple]
		del = '\Delta'

		xr = [3,8.5]
		yr = [3,8.5]

		Name		 = 'NScat'+GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		OutputDir	 = ODir + Name+'/'

		mean_nscat = fltarr(ntau)

		for i = 0,ntau-1 do begin
			ltau = log_Tau0Arr[i]
			t0 = 10^ltau

			xarr 		= fltarr(NMAX)
			xarr_abs 	= fltarr(NMAX)
			nscat	 	= fltarr(NMAX)
			nout		= 0l
			np			= 0l

			for j = 0,nParFiles[i]-1 do begin
				FTag = 'log_tau'+strn(ltau,len=3)+'.'+strn(j)
				OutFile = OutputDir + FTag
			
				openr,1,OutFile,err=err

				if err ne 0 then begin
					print,'File '+OutFile+' not found. Skipping it...'
					close,1
					continue
				endif
				print,OutFile
				str=''
				readf,1,str

				while ~EOF(1) do begin
					readf,1,str
					res = strsplit(str,/ext)
					nel = n_elements(res)	

					if res[0] eq '#' then continue

					if res[nel-1] eq 'Photon_escaped' then begin
						xarr[nout]  = res[1]		
						nscat[nout] = res[0]
						nout++	
					endif
					np++
				endwhile
				close,1
			endfor

			print,'Nphotons computed =',np

			if np lt 10 then begin
				print,'Only ',np,' photons found. Forget it...'
					continue
			endif

			xarr = xarr[0:nout-1]
			nscat = nscat[0:nout-1]

			T4 = Temp*1e-4
			a = 4.693e-4 * T4^(-1./2)
			a_tau0 = a*t0

			mean_nscat[i] = mean(nscat)
		endfor

		plot,indgen(2),/nodata,xr=xr,yr=yr,/xs,/ys,$
		position=aspect(1.0),xtitle=textoidl('log(\tau_0)'),$
		ytitle=textoidl('log(<N_{scat}>)'),charsiz=2,xthic=4,ythic=4

		plotsym,0,/fill
		taup = findgen(20)*0.5
		oplot,taup,alog10(nscat_harrington(10^taup)),lines=2,thic=3
		oplot,alog10(tau0arr),alog10(mean_nscat),ps=8,symsiz=2

		legend,'Harrington (1973)',/left,box=0,lines=2,psp=1.5,$
		charsiz=2,thic=3

		if PlotPs then device,/close
	endif


; #################################


	if Plot_fesc then begin
		PlotFile	 = 'fesc_HomSlab'
		GeomName 	 = 'HomSlab'
		log_tau0Arr	 = [6.]
		Temp		 = 10.		; [K]
		xcritarr	 = 3.

		ntau		 = n_elements(log_Tau0Arr)		
		
		NPhotons	 = 10000
		NParFiles	 = 10
		mean_nh		 = 100.

		minx		 = -2.0
		maxx		 = 1.0
		binx		 = 0.5
		nx			 = (maxx - minx)/binx + 1.0

		xaxis 		 = findgen(nx)*binx + minx

		xr = [-2.0	,	1.1]
		yr = [-3.	,	0.5]


		if PlotPs then begin
			set_plot,'ps'
			!P.Font = 0
			device,/color,/landscape,/palatino,filename	=PlotDir + PlotFile+'.ps'
		endif else begin
			window,wc,xsiz=650,ysiz=650,retain=2
			wc++
		endelse

		defplotcolors
		col = [!orange,!blue,!red,!dgreen,!purple]
		del = '\Delta'

		plot,indgen(2),/nodata,xr=xr,yr=yr,/xs,/ys,$
		xtitle=textoidl('log[(a\tau_0)^{1/3}\tau_a]'),$
		ytitle = textoidl('log(f_{esc})'),$
		charsiz=2,position=aspect(1.0)
	
		xcrit 		 = xcritarr
		Name		 = 'fesc_'+GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		OutputDir	 = ODir + Name+'/'

		ltau = log_Tau0Arr

		fesc = fltarr(nx)
		
		for iz = 0,nx-1 do begin
			t0 = 10^ltau

			xarr 		= fltarr(NMAX)
			xarr_abs 	= fltarr(NMAX)
			nscat	 	= fltarr(NMAX)
			nout		= 0l
			nabs		= 0l
			np			= 0l

			for j = 0,nParFiles-1 do begin
				FTag = 'xvar_'+strn(xaxis[iz],len=4)+'log_tau'+strn(ltau,len=3)+'.'+strn(j)
				OutFile = OutputDir + FTag
				
				openr,1,OutFile,err=err

				if err ne 0 then begin
					print,'File '+OutFile+' not found. Skipping it...'
					close,1
					continue
				endif
				print,OutFile
				str=''
				readf,1,str

				while ~EOF(1) do begin
					readf,1,str
					res = strsplit(str,/ext)
					nel = n_elements(res)	

					if res[0] eq '#' then continue

					if res[nel-1] eq 'Photon_escaped' then begin
						xarr[nout]  = res[1]		
						nscat[nout] = res[0]
						nout++	
					endif
					if res[nel-1] eq 'Absorbed_by_dust' then begin
						xarr_abs[nabs] = res[1]
						nabs++
					endif	
					np++
				endwhile
				close,1
			endfor

			print,'Nphotons computed =',np

			if np lt 10 then begin
				print,'Only ',np,' photons found. Forget it...'
				continue
			endif

			fesc[iz] = (nout + 0.)/(np + 0.)
		endfor
		plotsym,0,thic=5
		xax = findgen(50)*0.1 - 2.0
		fesc_neufeld = fesc_neufeld(10^xax)

		oplot,xax,alog10(fesc_neufeld),color = !orange		

		oplot,xaxis,alog10(fesc),ps=8,symsiz=2
		plotsym,0,/fill
		oplot,xaxis,alog10(fesc),ps=8,symsiz=1.6,color=!yellow
		
		print,fesc

		legend,'Neufeld (1990)',box=0,/right,lines=0,color=!orange,$
		psp=1.5,charsiz=2

		if PlotPs then device,/close
	endif



;#######

	if Plot_ExpSphere then begin
		PlotFile	 = 'J_ExpSphere'
		GeomName 	 = 'HomSphere'
		log_tau0Arr	 = 7.06
		Temp		 = 10000.		; [K]
		xcrit		 = 3.
		vmaxarr		 = [0.,20.,200.,2000.]

		nv			 = n_elements(vmaxarr)
	
		NParFiles	 = 30
		mean_nh		 = 100.

		xr = [-120,50]
		yr = [0,0.0050]


		if PlotPs then begin
			set_plot,'ps'
			!P.Font = 0
			device,/color,/landscape,/palatino,filename	=PlotDir + PlotFile+'.ps'
		endif else begin
			window,wc,xsiz=850,ysiz=650,retain=2
			wc++
		endelse

		defplotcolors
		col = [!orange,!blue,!red,!dgreen,!purple]
		del = '\Delta'

		plot,indgen(2),/nodata,xr=xr,yr=yr,/xs,/ys,$
		xtitle=textoidl('x = (\nu - \nu_0)/'+del+'\nu_D'),$
		ytitle = textoidl('J(\tau_0,x)'),$
		charsiz=2,position=aspect(.8),xthic=3,ythic=3
	
		Name		 = 'Expanding_'+GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		OutputDir	 = ODir + Name+'/'
		ltau = log_Tau0Arr
		t0 = 10^ltau

		for i = 0,nv-1 do begin

			xarr 		= fltarr(NMAX)
			xarr_abs 	= fltarr(NMAX)
			nscat	 	= fltarr(NMAX)
			nout		= 0l
			np			= 0l
			
			vmax		= vmaxarr[i]		
	
			for j = 0,nParFiles-1 do begin
				FTag = 'vmax_'+strn(vmax,len=4)+'log_tau'+strn(ltau,len=3)+'.'+strn(j)
				OutFile = OutputDir + FTag
				
				openr,1,OutFile,err=err

				if err ne 0 then begin
					print,'File '+OutFile+' not found. Skipping it...'
					close,1
					continue
				endif
				print,OutFile
				str=''
				readf,1,str

				while ~EOF(1) do begin
					readf,1,str
					res = strsplit(str,/ext)
					nel = n_elements(res)	

					if res[0] eq '#' then continue

					if res[nel-1] eq 'Photon_escaped' then begin
						xarr[nout]  = res[1]		
						nscat[nout] = res[0]
						nout++	
					endif
					np++
				endwhile
				close,1
			endfor

			print,'Nphotons computed =',np

			if np lt 10 then begin
				print,'Only ',np,' photons found. Forget it...'
				continue
			endif

			xarr = xarr[0:nout-1]
			nscat = nscat[0:nout-1]

			if np lt 100 then databin = 5.0
			if np gt 100 and np lt 1000 then databin = 2.5
			if np gt 1000 and np lt 5e3 then databin = 2.0
			if np gt 5e3 and np lt 1e4 then databin = 1.5
			if np gt 1e4 and np lt 5e4 then databin = 1.0
			if np ge 5e4 then databin = .5
				
			T4 = Temp*1e-4
			a = 4.693e-4 * T4^(-1./2)
			a_tau0 = a*t0



			if vmax eq 0. then begin
				xnmax = 120
				xnmin = -120
				xnbin = 0.1
				Nxneu = (xnmax - xnmin)/xnbin + 1

				xneu = findgen(Nxneu)*xnbin + xnmin
				
				jneu = J_Dijkstra(xneu,a_tau0)
				Norm_Neu = int_tabulated(xneu,jneu,/double)

				jneu = jneu/norm_Neu*(1./(4*!PI))
				mspec = max(jneu)
				dr = 1
				oplot,xneu,jneu,thic=3
			endif


			plothist,xarr,xdata,ydata,bin=databin,/noplot
			ydata = ((ydata+0.)/(np+0.)/databin)*(1./(4*!PI))

			oplot,xdata,ydata,thic=4,color=col[i],ps=10

		endfor

		print,'Reading L09 data'
		datadir = '/home/aaorsi/LyaRT/data/peter/Sphere/'
		dataFiles = datadir + $
		['SpDLAT4v0_x.dat','SpDLAT4v20_x.dat','SpDLAT4v200_x.dat','SpDLAT4v2000_x.dat']
		ndata = n_elements(dataFiles)

		
		col2 = [!tan,!lblue,!pink,!lgreen,!dpink]
		bin = 0.5
		for i = 0,ndata-1 do begin
			rdfloat,datafiles[i],x_l09,/sile
			np = n_elements(x_l09)
			plothist,x_l09,xdata,ydata,bin=bin,/noplot
			ydata = (ydata + 0.)/(np * bin) * (1./(4*!PI))
			oplot,xdata,ydata,color=!dgray,thic=5
;			oplot,xdata,ydata,color=col2[i],thic=3,lines=2
		endfor

		vstr = strarr(nv)
		for i = 0,nv-1 do vstr[i] = strn(fix(vmaxarr[i]))

		legend,textoidl('v_{max}[km/s] = '+vstr),box=0,/left,textcolor = col[0:nv-1],$
		charsiz=1.5 

		if PlotPs then device,/close
	endif

	if Plot_ThinShell then begin
		PlotFile	 = 'J_ThinShell'
		GeomName 	 = 'ThinShell'
		ltau		 = 6.5
		Temp		 = 96897.		; [K]
;		Temp		 = 10000.		; [K]
		xcrit		 = 0.0
		vmax		 = 300.

		NParFiles	 = 30

		xr = [-50,20]
		yr = [0,0.07]


		if PlotPs then begin
			set_plot,'ps'
			!P.Font = 0
			device,/color,/landscape,/palatino,filename	=PlotDir + PlotFile+'.ps'
		endif else begin
			window,wc,xsiz=850,ysiz=650,retain=2
			wc++
		endelse

		defplotcolors
		col = [!orange,!blue,!red,!dgreen,!purple,!dpink]
		del = '\Delta'

	
		Name		 = GeomName + '_T'+strn(Temp,len=5)+'_xcrit'+strn(xcrit,len=3)
		OutputDir	 = ODir + Name+'/'


		xarr 		= fltarr(NMAX)
		xarr_abs 	= fltarr(NMAX)
		nscat	 	= fltarr(NMAX)
		flag_zero	= fltarr(nmax)
		nout		= 0l
		np			= 0l
		
		for j = 0,nParFiles-1 do begin



			FTag = 'vmax_'+strn(vmax,len=4)+'log_tau'+strn(ltau,len=3)+'.'+strn(j)
			OutFile = OutputDir + FTag
				
			openr,1,OutFile,err=err

			if err ne 0 then begin
				print,'File '+OutFile+' not found. Skipping it...'
				close,1
				continue
			endif
			print,OutFile
			str=''
			readf,1,str

			while ~EOF(1) do begin
				readf,1,str
				res = strsplit(str,/ext)
				nel = n_elements(res)	

				if res[0] eq '#' then continue
				if res[nel-1] eq 'Photon_escaped' then begin
				xarr[nout]  = res[1]		
				nscat[nout] = res[0]
				flag_zero[nout] = res[nel-2]
				nout++	
				endif
				np++
			endwhile
			close,1

;			xarr = xarr[0:nout-1]
;			nscat = nscat[0:nout-1]

			if np lt 100 then databin = 5.0
			if np gt 100 and np lt 1000 then databin = 2.5
			if np gt 1000 and np lt 5e3 then databin = 2.0
			if np gt 5e3 and np lt 1e4 then databin = 1.5
			if np gt 1e4 and np lt 5e4 then databin = 1.0
			if np ge 5e4 then databin = .5
	

		endfor
			databin = 1.0
	

			plot,indgen(2),/nodata,xr=xr,yr=yr,/xs,/ys,$
			xtitle=textoidl('x = (\nu - \nu_0)/'+del+'\nu_D'),$
			ytitle = textoidl('J(\tau_0,x)'),$
			charsiz=2,position=aspect(.6)
			plothist,xarr[0:nout-1],xdata,ydata,bin=databin,/noplot
			ydata = ((ydata+0.)/(np+0.)/databin)
	
			oplot,xdata,ydata,thic=4,color=col[0],ps=10

			i0 = where(flag_zero[0:nout-1] eq 1,n0)
			i1 = where(flag_zero[0:nout-1] eq 2,n1)
			i2 = where(flag_zero[0:nout-1] eq 3,n2)
			i3 = where(flag_zero[0:nout-1] ge 4,n3)
			
			if n0 gt 3 then begin

				plothist,xarr[i0],xdata,ydata,bin=databin,/noplot
				ydata = ((ydata+0.)/(np+0.)/databin)
				oplot,xdata,ydata,thic=1,color=col[1]
			endif

			if n1 gt 3 then begin

				plothist,xarr[i1],xdata,ydata,bin=databin,/noplot
				ydata = ((ydata+0.)/(np+0.)/databin)
				oplot,xdata,ydata,thic=1,color=col[2]
			endif

			if n2 gt 3 then begin

				plothist,xarr[i2],xdata,ydata,bin=databin,/noplot
				ydata = ((ydata+0.)/(np+0.)/databin)
				oplot,xdata,ydata,thic=1,color=col[3]
			endif

			if n3 gt 3 then begin

				plothist,xarr[i3],xdata,ydata,bin=databin,/noplot
				ydata = ((ydata+0.)/(np+0.)/databin)
				oplot,xdata,ydata,thic=1,color=col[4]

			endif

		legend,['Thin Shell'],box=0,/right,/top
		legend,textoidl(['Full spectrum',' 0 backscattering','1 backscattering','2 backscatterings','3 backscatterings or more']),$
		lines=0,thic=3,color=col[0:4],box=0,/left
		
	
		print,'Nphotons computed =',np


		if PlotPs then device,/close
	endif


		
end	
