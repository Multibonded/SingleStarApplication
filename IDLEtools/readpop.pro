function readpop,file,nx=nx,ny=ny,nz=nz,nlevel=nlevel
    ;---------------------------------------------------------------------------
    ; name:
    ;   readpop
    ; purpose:
    ;   reads populations
    ; comments:
    ;   none yet
    ;---------------------------------------------------------------------------

    openr,lun,file,/get_lun
            
    b=fltarr(nx,ny,nz,nlevel)
    nstar=fltarr(nx,ny,nz,nlevel)
    totn=fltarr(nx,ny,nz)
    readu,lun,b
    readu,lun,nstar
    readu,lun,totn
    free_lun,lun

    b=double(b)
    nstar=10d0^double(nstar)
    totn=10d0^double(totn)
    n=nstar*b

    pop = { n:n, nstar:nstar, totn:totn, b:b }
    return,pop

end

