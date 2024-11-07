function readie,file,out=out,mode=mode
    ;---------------------------------------------------------------------------
    ; 2014 ama51:
    ; name:
    ;   readie
    ; purpose:
    ; comments:
    ;---------------------------------------------------------------------------
    if n_elements(mode) eq 0l then mode='avrg_xy'

    openr,lun,file,/get_lun
    case mode of
    'resolved':begin
        ie=dblarr(out.nx,out.ny,out.nangle,out.nnu)
        readu,lun,ie
    end
    'avrg_xy':begin
        ie=dblarr(out.nangle,out.nnu)
        readu,lun,ie
    end
    'avrg_xy_phi':begin
        ie=dblarr(out.nmuav,out.nnu)
        readu,lun,ie
    end
    endcase

    free_lun,lun

    return,ie
end

