function readck,temp=temp,mh=mh,lg=lg

	if mh eq -25 then str='ckm25'
	if mh eq -20 then str='ckm20'
	if mh eq -15 then str='ckm15'
	if mh eq -10 then str='ckm10'
	if mh eq -5 then str='ckm05'
	if mh eq 0 then str='ckp00'
	if mh eq 2 then str='ckp02'
	if mh eq 5 then str='ckp05'
	temp=strtrim(string(long(temp)),1)
	star_dir = 'Star/castelli/'+str+'/'+str+'_'+temp+'.fits'
	ckdata=mrdfits(star_dir,1,hdr,/silent)
	ckwave=ckdata.wavelength
	
	if lg eq 0.0 then ckflux=ckdata.g00
	if lg eq 0.5 then ckflux=ckdata.g05
	if lg eq 1.0 then ckflux=ckdata.g10
	if lg eq 1.5 then ckflux=ckdata.g15
	if lg eq 2.0 then ckflux=ckdata.g20
	if lg eq 2.5 then ckflux=ckdata.g25
	if lg eq 3.0 then ckflux=ckdata.g30
	if lg eq 3.5 then ckflux=ckdata.g35
	if lg eq 4.0 then ckflux=ckdata.g40
	if lg eq 4.5 then ckflux=ckdata.g45
	if lg eq 5.0 then ckflux=ckdata.g50
	
	ckdata=replicate({wave: 0., flux: 0.},n_elements(ckwave))
	
	for i=0,n_elements(ckwave)-1 do begin
		ckdata[i].wave=ckwave[i]
		ckdata[i].flux=ckflux[i]
	endfor
return,ckdata
end
