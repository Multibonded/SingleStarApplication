function readdet,file,sp=sp
    ;---------------------------------------------------------------------------
    ; name:
    ;   readdet
    ; purpose:
    ;   read detailed output parameters
    ; comments:
    ;   none yet
    ;---------------------------------------------------------------------------
    @const
    
    openr,lun,file,/get_lun

    ; read frequency indices
    nnu=0l
    readu,lun,nnu
    off=lonarr(nnu)
    readu,lun,off
    if n_elements(sp) gt 0l then wl=sp.wl[off-1l] else wl=-1
    freq = { nnu:nnu, off:off, wl:wl }

    ; read lines
    nline=0l
    readu,lun,nline
    if nline gt 0l then iatoml=lonarr(nline) else iatoml=-1l
    if nline gt 0l then wl0=lonarr(nline) else wl0=-1l
    for iline=0l,nline-1l do begin
        jatoml=0l
        readu,lun,jatoml
        iatoml[iline]=jatoml
        lam=0d
        readu,lun,lam
        wl0[iline]=lam*1d-2 ; m
    endfor

    ; read continua
    ncont=0l
    readu,lun,ncont
    if ncont gt 0l then iatomc=lonarr(ncont) else iatomc=-1l
    if ncont gt 0l then nu0=lonarr(ncont) else nu0=-1l
    for icont=0l,ncont-1l do begin
        jatomc=0l
        readu,lun,jatomc
        iatomc[icont]=jatomc
        lam=0d
        readu,lun,lam
        nu0[icont]=lam ; s^-1
    endfor

    ; read angles
    nangle=0l
    readu,lun,nangle
    if nangle gt 0l then begin
        mux=dblarr(nangle)
        muy=dblarr(nangle)
        muz=dblarr(nangle)
        wmu=dblarr(nangle)
        readu,lun,mux
        readu,lun,muy
        readu,lun,muz
        readu,lun,wmu
    endif else begin
        mux=-1l
        muy=-1l
        muz=-1l
        wmu=-1l
    endelse

    nmuav=0l
    readu,lun,nmuav
    if nmuav gt 0l then begin
        muav=dblarr(nmuav)
        wmuav=dblarr(nmuav)
        nphiav=lonarr(nmuav)
        readu,lun,muav
        readu,lun,wmuav
        readu,lun,nphiav
    endif else begin
        muav=-1l
        wmuav=-1l
        nphiav=-1l
    endelse
    
    ; read geo
    nx=0l
    ny=0l
    nz=0l
    ngeo=0l
    readu,lun,nx
    readu,lun,ny
    readu,lun,nz
    readu,lun,ngeo
    if ngeo gt 0l then begin
        geo_mux=dblarr(ngeo)
        geo_muy=dblarr(ngeo)
        geo_muz=dblarr(ngeo)
        geo_wmu=dblarr(ngeo)
        readu,lun,geo_mux
        readu,lun,geo_muy
        readu,lun,geo_muz
        readu,lun,geo_wmu
    endif else begin
        geo_mux=-1l
        geo_muy=-1l
        geo_muz=-1l
        geo_wmu=-1l
    endelse

    free_lun,lun

    out={nnu:nnu,off:off,wl:wl,nline:nline,iatoml:iatoml,wl0:wl0,$
        ncont:ncont,iatomc:iatomc,nu0:nu0,$
        nangle:nangle,mux:mux,muy:muy,muz:muz,wmu:wmu,$
        nmuav:nmuav,muav:muav,wmuav:wmuav,nphiav:nphiav,$
        nx:nx,ny:ny,nz:nz,ngeo:ngeo,$
        geo_mux:geo_mux,geo_muy:geo_muy,geo_muz:geo_muz,geo_wmu:geo_wmu}

    return,out
end
