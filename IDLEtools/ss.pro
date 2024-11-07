pro ss,ion=ion
    star='marcs_sun'
    DIR='/crex/proj/snic2020-16-23/private/Jack/template_balder/project/HD140283/TiRun_Nq140k_k0001_250Abund/output'
    sp=readsp(dir+'/sp2')
    OUT=readdet(dir+'/detailed_'+star,sp=sp)
    IET=readie(dir+'/iet_'+star,out=out)
    IEC=readie(dir+'/iec_'+star,out=out)
    IETL=readie(dir+'/ietl_'+star,out=out)
    IECL=readie(dir+'/iecl_'+star,out=out)

    FLUXT=dblarr(out.nnu)
    FLUXC=dblarr(out.nnu)
    FLUXTL=dblarr(out.nnu)
    FLUXCL=dblarr(out.nnu)

    for inu=0l,out.nnu-1l do begin
        FLUXT[inu]=total(IET[*,inu]*out.muz*out.wmu)
        FLUXC[inu]=total(IEC[*,inu]*out.muz*out.wmu)
        FLUXTL[inu]=total(IETL[*,inu]*out.muz*out.wmu)
        FLUXCL[inu]=total(IECL[*,inu]*out.muz*out.wmu)
    endfor

    s=sort(out.wl)
    wl=out.wl[s]*1d9
    IET=IET[*,s]
    IEC=IEC[*,s]
    IETL=IETL[*,s]
    IECL=IECL[*,s]
    FLUXT=FLUXT[s]
    FLUXC=FLUXC[s]
    FLUXTL=FLUXTL[s]
    FLUXCL=FLUXCL[s]
    print, total(FLUXT)
    print, total(FLUXC)

    print, total(IET[-1, *])
    print, total(IEC[-1, *])
    print, "Ratio of flux", total(FLUXT)/total(FLUXC)
    print, "Ratio of intensity", total(IET[-1, *])/total(IEC[-1, *])
    plot,wl,FLUXT/FLUXC,xr=[350,450]
    oplot,wl,FLUXTL/FLUXCL,linestyle=1
    
end