PRO psfull,type,name=name,eps=eps

;colors
!p.charsize=1.0
!p.font=3
!x.thick=3
!y.thick=3
!p.thick=1

ind=keyword_set(name)
if ind eq 0 then name='idl.ps'
print,'Starting plot to: '+name

if type eq 'l' then begin
if keyword_set(eps) then begin 
set_plot,'ps'
device,filename=name,xsize=26,ysize=20,xoffset=0.0,/landscape,/color,/encapsulated,bits_per_pixel = 8,/ISOLATIN1
endif else begin 
set_plot,'ps'
device,filename=name,xsize=26,ysize=20,xoffset=0.0,/landscape,/color,bits_per_pixel = 8,/ISOLATIN1
endelse
endif

if type eq 'p' then begin
if keyword_set(eps) then begin
set_plot,'ps'
device,filename=name,xsize=20,ysize=26,xoffset=0.0,yoffset=0.0,/color,/encapsulated,bits_per_pixel = 8,/ISOLATIN1,/portrait
endif else begin
set_plot,'ps'
device,filename=name,xsize=20,ysize=26,xoffset=0.0,yoffset=0.0,/color,bits_per_pixel = 8,/ISOLATIN1,/portrait
endelse
endif

END
