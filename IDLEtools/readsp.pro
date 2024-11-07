function readsp,file
    ;---------------------------------------------------------------------------
    ; name:
    ;   readsp
    ; purpose:
    ;   read sp file
    ; comments:
    ;   none yet
    ;---------------------------------------------------------------------------
    @const

    openr,lun,file,/get_lun
    nnu=0l
    maxac=0l
    maxal=0l
    readu,lun,nnu
    readu,lun,maxac
    readu,lun,maxal

    nu=dblarr(nnu)
    wnu=dblarr(nnu)
    readu,lun,nu
    readu,lun,wnu

    if (maxac gt 0) then ac=lonarr(maxac,nnu) else ac=-1
    if (maxac gt 0) then readu,lun,ac

    if (maxal gt 0) then al=lonarr(maxal,nnu) else al=-1
    if (maxal gt 0) then readu,lun,al

    nac=lonarr(nnu)
    nal=lonarr(nnu)
    readu,lun,nac
    readu,lun,nal

    free_lun,lun

    wl=const.cc/nu
    sp = { nnu:nnu, maxac:maxac, maxal:maxal, nu:nu, wnu:wnu,$
        ac:ac, al:al, nac:nac, nal:nal, wl:wl}

    return,sp
end
