function FM89, x, p
;*****************************************************************************
	bump	= p[0]
	gamma	= p[1]
	c1	= p[2]
	c2	= p[3]
	c3	= p[4]
	c4	= p[5]
	c5	= p[6]
	fit	= fltarr(1, n_elements(x))

	for i = 0, n_elements(x) - 1 do begin
		xlin = c1 + c2 * x[i]
		xbum = x[i]^2 / ((x[i]^2 - bump^2)^2 + (gamma * x[i])^2)
		fit[i] = xlin + c3 * xbum
		if x[i] gt 5.9 then fit[i]=fit[i]+$
			c4*(x[i]-c5)^2
	endfor

	return, fit

;***************************************************************************
end
;***************************************************************************
;MAIN Program
pro fitsed,star_name,nstar=nstar,min_temp=min_temp,max_temp=max_temp,$
	min_mh=min_mh,max_mh=max_mh,min_lg=min_lg,max_lg=max_lg,$
	nhpa=nhpa,ines=ines,nhiue=nhiue,dtemp=dtemp,quiet=quiet,bin=bin,$
	doublewt=doublewt
;***************************************************************************
;Default Values
	fixrv=0
	dobin=0
	if n_elements(nstar) eq 0 then nstar=1
	if n_elements(min_temp) eq 0 then min_temp=[11000.,11000.]
	if n_elements(max_temp) eq 0 then max_temp=[35000.,35000.]
	if n_elements(min_mh) eq 0 then min_mh=[0,0]
	if n_elements(max_mh) eq 0 then max_mh=[0,0]
	if n_elements(min_lg) eq 0 then min_lg=[3.5,3.5]
	if n_elements(max_lg) eq 0 then max_lg=[5.,5.]
	if n_elements(min_ebv) eq 0 then min_ebv=-0.05
	if n_elements(max_ebv) eq 0 then max_ebv=2.
	if n_elements(nhpa) eq 0 then nhpa=0
	if n_elements(ines) eq 0 then ines=0
	if n_elements(nhiue) eq 0 then nhiue=0
	if n_elements(dtemp) eq 0 then dtemp=[1000,1000]
	if n_elements(quiet) eq 0 then quiet=1
	if n_elements(bin) eq 1 then dobin=1
	if n_elements(doublewt) eq 0 then doublewt=0
	if n_elements(min_temp) eq 1 then min_temp=[min_temp,13000]
	if n_elements(max_temp) eq 1 then max_temp=[max_temp,13000]
	if n_elements(min_lg) eq 1 then min_lg=[min_lg,0.]
	if n_elements(max_lg) eq 1 then max_lg=[max_lg,0.]
	if n_elements(min_mh) eq 1 then min_mh=[min_mh,0.]
	if n_elements(max_mh) eq 1 then max_mh=[max_mh,0.]
	tr1=1.
	tr2=1.
	chisq=1e10
	itr=0
	conv=4.8481368110954E-9
	fv_vega=3.67E-9
	fb_vega=6.4E-9
	m_vega=0.
;***************************************************************************
;Get star data
	if nhpa eq 0 and nhiue eq 0 then begin
		if ines eq 1 then restore,'Star_Sav/INES/'+star_name+'_ines.sav'
		if ines eq 0 then restore,'Star_Sav/NEWSIPS/'+star_name+'_newsips.sav'
		wave=iue.wave
		star=iue.flux
		estar=iue.sigma
		;plot,wave,(star),xrange=[1000,3200],psy=3
	endif
	if nhpa eq 1 then begin
		restore,'Star_Sav/NH/'+star_name+'.sav'
		wave=nh.wave
		star=nh.flux
		estar=nh.sigma
		star=star*6.625*2.99792*1e-9/wave^2
		estar=estar*6.625*2.99792*1e-9/wave^2
		;plot,wave,star,xrange=[900,2000],psym=3
	endif
	if nhiue eq 1 then begin	
		restore,'Star_Sav/INES/'+star_name+'_ines.sav'
		iuewave=iue.wave
		iueflux=iue.flux
		iuesigma=iue.sigma
		restore,'Star_Sav/NH/'+star_name+'.sav'
		nhwave=nh.wave
		nhstar=nh.flux
		nhestar=nh.sigma
		i1=min(where(nhwave ge 900))
		i2=max(where(nhwave le 1800))
		nhwave=nhwave[i1:i2]
		nhstar=nhstar[i1:i2]
		nhestar=nhestar[i1:i2]
		nhstar=nhstar*6.625*2.99792/1e9/nhwave
		nhestar=nhestar*6.625*2.99792/1e9/nhwave
		wave=[iuewave,nhwave]
		star=[iueflux,nhstar]
		estar=[iuesigma,nhestar]
		s=sort(wave)
		wave=wave[s]
		star=star[s]
		estar=estar[s]
		;plot,wave,star,xrange=[800,3200],psym=3
	endif
	fv=fv_vega*10^(-0.4*(V-m_vega))

	if dobin eq 1 then begin
		nbin=(max(wave)-min(wave))/bin
		sbin=bin/2
		w1=fltarr(1,nbin)
		f1=fltarr(1,nbin)
		e1=fltarr(1,nbin)
		for i=0,nbin-1 do begin
			if i eq 0 then x=min(wave)+sbin else x=x+bin
			fb=star(where(wave ge x-sbin and wave lt x+sbin))
			eb=estar(where(wave ge x-sbin and wave lt x+sbin))
			wb=wave(where(wave ge x-sbin and wave lt x+sbin))
			f1[i]=total(fb/eb^2)
			w1[i]=total(wb/eb^2)
			e1[i]=total(1/eb^2)
		endfor
		f1=f1/e1
		w1=w1/e1
		e1=sqrt(1/e1)
		s=sort(w1)
		w1=w1[s]
		f1=f1[s]
		e1=e1[s]
		wave=fltarr(1,nbin)
		star=fltarr(1,nbin)
		estar=fltarr(1,nbin)
		wave[0:nbin-1]=w1[0:nbin-1]
		star[0:nbin-1]=f1[0:nbin-1]
		estar[0:nbin-1]=e1[0:nbin-1]
	endif
;***************************************************************************
;FIT RANGE
;Donot fit data in the range 1100A-1300A
;Select only good data, here where ever flux or sigma eq +/-NaN removed
	
	ws=wave(where(FINITE(star*estar) ne 0 and estar gt 0))
	fs=star(where(FINITE(star*estar) ne 0 and estar gt 0))
	es=estar(where(FINITE(star*estar) ne 0 and estar gt 0))
	ws=ws(where(fs*es ne 0))
	fs=fs(where(fs*es ne 0))
	es=es(where(fs*es ne 0))
	
	i1 = min(where(ws gt 900))
	i2 = max(where(ws lt 3200))
	es[0:i1] = 1e6
	es[i2:*] = 1e6
	i1 = min(where(ws gt 1150))
	i2 = max(where(ws lt 1250))
	es[i1:i2] = 1e6
	
	good_data=where(es lt 1)
	w=ws[good_data]
	f=fs[good_data]
	s=es[good_data]
	ndeg=total(good_data)
	deg=n_elements(good_data)
	
	if nstar eq 1 then begin
		min_temp[1]=max_temp[1]
		min_mh[1]=max_mh[1]
		min_lg[1]=max_lg[1]
	endif
	
	wt=(1/s^2)/(total(1/s^2)/n_elements(s))
	avgvar=1/(total(1/s^2)/n_elements(s))
	if doublewt eq 1 then begin
		i1 = min(where(ws gt 1600))
		wt[i1:*] = wt[i1:*]*2
	endif
	wset,0
	plot,w,f
	errplot,w,f-s,f+s
;***************************************************************************
	t1=double(max_temp[0])
	while t1 ge double(min_temp[0]) do begin
		m1=min_mh[0]
		while m1 le max_mh[0] do begin
			for l1=min_lg[0],max_lg[0],0.5 do begin	
				ckdata=readck(temp=long(t1),mh=long(m1),lg=long(l1))
				ckwave=ckdata.wave
				ckflux1=ckdata.flux	
				t2=double(max_temp[1])
				while t2 ge double(min_temp[1]) do begin
					m2=min_mh[1]
					while m2 le max_mh[1] do begin
						for l2=min_lg[1],max_lg[1],0.5 do begin
							ckdata=readck(temp=long(t2),mh=long(m2),lg=long(l2))
							ckflux2=ckdata.flux
							ckflux=ckflux1+(nstar-1)*ckflux2
							ckv=ckflux[min(where(ckwave ge 5500))]
							ckb=ckflux[min(where(ckwave ge 4450))]
							V0=2.5*alog10(fv_vega/ckv)
							B0=2.5*alog10(fb_vega/ckb)
							B_V0=B0-V0
							i1 = min(where(ckwave gt 800))
							i2 = max(where(ckwave lt 3300))
							ckwave = ckwave[i1:i2]
							ck = ckflux[i1:i2]
							quadterp,ckwave,ck,w,ck
							ebv=B_V-B_V0
							kext=(2.5*alog10(fv/f)-2.5*alog10(ckv/ck))/ebv
							n=(max(1e4/w)-min(1e4/w))/0.05
							abin=fltarr(1,n)
							kbin=fltarr(1,n)
							kwt=fltarr(1,n)
							for ibin=0,n-1 do begin
								if ibin eq 0 then xbin=min(1e4/w)+0.025 else xbin=xbin+0.05
								abin[ibin]=1e4/mean(w(where(1e4/w ge xbin-0.025 and 1e4/w lt xbin+0.025)))
								fbin=mean(f(where(1e4/w ge xbin-0.025 and 1e4/w lt xbin+0.025)))
								ckbin=mean(ck(where(1e4/w ge xbin-0.025 and 1e4/w lt xbin+0.025)))
								kbin[ibin]=(2.5*alog10(fv/fbin)-2.5*alog10(ckv/ckbin))/ebv
								if abin[ibin] lt 5.9 then kwt[ibin]=2 else kwt[ibin]=1
							endfor
							sbin=sort(abin)
							abin=abin[sbin]
							kbin=kbin[sbin]
							kwt=kwt[sbin]
							par	= replicate({value: 0.0d, fixed: 0, limited: [0,0], $
								limits: [0.D,0.D], step: 0.0d}, 7)
							par[0:6].value=[4.595,1.051,-0.384,0.739,3.961,0.265,5.9]
							kpar = mpfitfun('FM89', abin, kbin, WEIGHTS = kwt, perror = dkpar,$
								yfit = yfit, PARINFO = par, quiet = 1, dof = DOF)
							rv=4.717/(kpar[3]+0.824)
							kext=FM89(1e4/w,kpar)
							tr=sqrt(3.6/ckv*10^(-9.-0.4*(-ebv*rv+V)))
							ck=ck*tr^2*10^(-0.4*ebv*(kext+rv))
							varnew=total(wt*(f-ck)^2)/(n_elements(w))
							chisq1=varnew/avgvar
							if chisq1 le chisq then begin
								sed_sav=[t1,t2,m1,m2,l1,l2,tr/conv,ebv,rv]
								print,sed_sav,chisq1
								kpar_sav=kpar
								dkpar_sav=dkpar
								chisq=chisq1
								oplot,w,ck,col=cgcolor('red')
							endif
						endfor
						dm2=5
						if m2 eq 0 then dm2=2
						if m2 eq 2 then dm2=3
						m2=m2+dm2
					endwhile
					t2=long(t2)
					if t2 le 13000 then dtemp[1]=250
					t2=t2-dtemp[1]
				endwhile
			endfor
			dm1=5
			if m1 eq 0 then dm1=2
			if m1 eq 2 then dm1=3
			m1=m1+dm1
		endwhile
		t1=long(t1)
		if t1 le 13000. then dtemp[0]=250
		dtemp[1]=1000.
		t1=t1-dtemp[0]
	endwhile
;***************************************************************************
;Print Output
	print,'BESTFIT'
	print,'SED VALUES:',sed_sav
	print,'EXTINCTION:',kpar_sav
	print,'ERROR	 :',dkpar_sav
	print,'CHISQ	:',chisq
;***************************************************************************
;Plotting
;SED Plotting
	if nhpa eq 1 then xrange=[600,1800]
	if nhiue eq 1 then xrange=[900,3200]
	if nhpa eq 0 and nhiue eq 0 then xrange=[1100,3200]
	angstrom = '!6!sA!r!u!9 %!6!n'
	yunit	 = 'Flux ( erg\ cm$\up2$ \ s\ '+angstrom+' )'
	window,0,xs=1024,ys=512

	ckdata=readck(temp=sed_sav[0],mh=sed_sav[2],lg=sed_sav[4])
	ckwave=ckdata.wave
	ckflux1=ckdata.flux	
	ckdata=readck(temp=sed_sav[1],mh=sed_sav[3],lg=sed_sav[5])
	ckflux2=ckdata.flux
	ckflux=ckflux1+(nstar-1)*ckflux2
	ckv=ckflux[min(where(ckwave ge 5500))]
	ckb=ckflux[min(where(ckwave ge 4450))]
	i1 = min(where(ckwave gt 800))
	i2 = max(where(ckwave lt 3300))
	ckwave = ckwave[i1:i2]
	ckflux = ckflux[i1:i2]
	quadterp,ckwave,ckflux,wave,cksedflux
	ckfit=cksedflux*(sed_sav[6]*conv)^2*10^(-0.4*sed_sav[7]*(FM89(1e4/wave,kpar_sav)+sed_sav[8]))
	
	cgplot,wave,star,xtitle='Wavelength '+angstrom,ytitle=yunit,$
		xrange=xrange,title=star_name,yrange=[0,max(star+estar)]
	errplot,wave,star-estar,star+estar,color=cgcolor('black')
	oplot,wave,ckfit,col=255,thick=2
	image=cgSnapshot(filename='Plots/SED/'+star_name+'.png',$
 		quality=100,nodialog=1)
		

	V0=2.5*alog10(fv_vega/ckv)
	B0=2.5*alog10(fb_vega/ckb)
	B_V0=B0-V0
	ebv=B_V-B_V0
	quadterp,ckwave,ckflux,w,ckflux
	kext=(2.5*alog10(fv/f)-2.5*alog10(ckv/ckflux))/ebv
;Extinction curve plotting
	l	= cgGreek('lambda')
	xunit	= l+'$\up-1$'+' ( $\mu$ m $\up-1$ )'
	yunit	= 'k ( '+l+' - V )'
	window,1,xs=512,ys=512
	ymax=max(kext(where(1e4/w ge 2 and 1e4/w le 10)))
	cgplot,1e4/w,kext,xtitle=xunit,ytitle=yunit,aspect=1,title=star_name,$
		yrange=[2,ymax],xrange=[2,10]
	oplot,1e4/w,FM89(1e4/w,kpar_sav),col=255,thick=2
	image=cgSnapshot(filename='Plots/Extinction/'+star_name+'.png',$
 		quality=100,nodialog=1)
;********************************************************************************
;Save Data
	save,ckfit,nstar,sed_sav,kpar_sav,dkpar,filename='Star_Sav/ckfit/'+$
		star_name+'_ckfit.sav'			
end