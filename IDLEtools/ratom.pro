function nextline,unit
    text=''
    while (~EOF(unit)) do begin
        readf,unit,text
        if (strmid(text,0,1) ne '*') then break
    endwhile
    return,text
end

function ratom,atomfile
    ;---------------------------------------------------------------------------
    ; 2015 ama51
    ; name:
    ;   ratom
    ; purpose:
    ;   reads in m3d model atom 
    ; comments:
    ;   na
    ;---------------------------------------------------------------------------
    @const

    max_nq=10000l

    atomfiles=[atomfile,'atom.'+atomfile]
    atomfiles=file_search(atomfiles,count=count)
    if(count ne 1) then begin
        print,atomfile+' not found'
        stop
    endif

    text=''
    openr,lu1,atomfiles[0],/get_lun
        text=nextline(lu1)
        species=''
        reads,text,species
        species=strcompress(species)

        text=nextline(lu1)
        abund=0d
        awgt=0d
        reads,text,abund,awgt

        text=nextline(lu1)
        nlevel=0l
        nline=0l
        ncont=0l
        nfix=0l
        reads,text,nlevel,nline,ncont,nfix

        ev=0d
        g=0d
        label=''
        ion=0l
        s=0d
        l=0d
        p=-1l
        j=0d
        landeg=0d
        ;dum={ev:ev,g:g,label:label,ion:ion,j:j}
        dum={ev:ev,g:g,label:label,ion:ion,s:s,l:l,p:p,j:j,landeg:landeg}
        level=replicate(dum,nlevel)
        for ilevel=0l,nlevel-1l do begin
            text=nextline(lu1)

            str=strsplit(text,"'",/extract)
            level[ilevel].label=strtrim(str[1],2)

            reads,str[0],ev,g
            level[ilevel].ev=cmm1_to_ev(ev)
            level[ilevel].g=g

            str=strsplit(str[2],/extract)
            nl=0l
            s=-1d
            l=-1d
            p=-1l
            j=-1d
            landeg=-1d
            if n_elements(str) eq 1l then begin
                reads,str,ion
            endif else if n_elements(str) eq 2l then begin
                reads,str,ion,nl
            endif else if n_elements(str) eq 3l then begin
                reads,str,ion,nl,j
            endif else if n_elements(str) eq 7l then begin
                reads,str,ion,s,l,p,j,landeg,nl
            endif else stop
            level[ilevel].ion=ion
            level[ilevel].j=j
            level[ilevel].s=s
            level[ilevel].l=l
            level[ilevel].p=p
            level[ilevel].landeg=landeg
        endfor

        if (nline gt 0) then begin
            ul=0l
            ll=0l
            f=0d
            nq=0l
            qmax=0d
            q0=0d
            ncomp=0l
            ga=0d
            gw=0d
            gs=0d
            wl=0d
            mcomp=15l
            dum={ul:ul,ll:ll,f:f,nq:nq,qmax:qmax,q0:q0,ncomp:ncomp,$
                ga:ga,gw:gw,gs:gs,wl:wl,denergy:dblarr(mcomp),weight:dblarr(mcomp)}
            line=replicate(dum,nline)
            for iline=0l,nline-1l do begin
                text=nextline(lu1)
                reads,text,ul,ll,f,nq,qmax,q0,ncomp,ga,gw,gs
                line[iline].ul=ul
                line[iline].ll=ll
                line[iline].f=f
                line[iline].nq=nq
                line[iline].qmax=qmax
                line[iline].q0=q0
                line[iline].ncomp=ncomp
                line[iline].ga=ga
                line[iline].gw=gw
                line[iline].gs=gs
                line[iline].wl=const.hh*const.cc/$
                    (const.ee*(level[ul-1l].ev-level[ll-1l].ev))
                if ncomp gt 1l then begin
                    for icomp=0l,ncomp-1l do begin
                        denergy=0d
                        weight=0d
                        text=nextline(lu1)
                        ;reads,text,denergy,weight
                        str=strsplit(text,/extract)
                        denergy=double(str[0])
                        weight=double(str[1])
                        line[iline].denergy[icomp]=denergy
                        line[iline].weight[icomp]=weight
                    endfor
                endif

            endfor
        endif else line=-1l

        if (ncont gt 0) then begin
            ul=0l
            ll=0l
            f=0d
            nq=0l
            qmax=0d
            wl=0d
            al=0d
            dum={ul:ul,ll:ll,f:f,nq:nq,qmax:qmax,wavelen:dblarr(max_nq),alpha:dblarr(max_nq)}
            cont=replicate(dum,ncont)
            for icont=0l,ncont-1l do begin
                text=nextline(lu1)
                reads,text,ul,ll,f,nq,qmax
                cont[icont].ul=ul
                cont[icont].ll=ll
                cont[icont].f=f
                cont[icont].nq=nq
                cont[icont].qmax=qmax
                if qmax lt 0d then begin
                    for iq=0l,nq-1l do begin
                        text=nextline(lu1)
                        reads,text,wl,al
                        cont[icont].wavelen[iq]=wl
                        cont[icont].alpha[iq]=al
                    endfor
                endif

            endfor
        endif else cont=-1

        text=nextline(lu1)
        colroutine=text

        mcol=file_lines(atomfiles[0])
        coldata=strarr(mcol)
        icol=0l
        while (~EOF(lu1)) do begin
            text=nextline(lu1)
            coldata[icol]=text
            icol+=1l
            if (text eq 'END') then break
        endwhile
        ncoldata=icol
        if ncoldata gt 0l then coldata=coldata[0l:ncoldata-1l] else coldata=-1
    
    ; close 
    free_lun,lu1

    ; atom struct
    atom={species:species,$
        abund:abund, awgt:awgt,$
        nlevel:nlevel, nline:nline, ncont:ncont, nfix:nfix,$
        level:level,$
        line:line,$
        cont:cont,$
        colroutine:colroutine,$
        ncoldata:ncoldata,$
        coldata:coldata}

    ; return atom struct
    return,atom

end
