PRO termdiag,file,ip,ev=ev,yr=yr,maxl=maxl,minf=minf,bb=bb,bf=bf,cid=cid,conf=conf,cut=cut,ps=ps,xyouts=xyouts,ytitle=ytitle

rtc=1.0973731569d5 ; Conversion factor between Rydberg and cm^(-1)
rte=13.605698d0    ; Conversion factor between Rydberg and eV 

L=['S','P','D','F','G','H','I','K','L','M','N','O','Q','R','T','U']

fits_info,file,n_ext=n_ext
print, "next", n_ext
ext1 = mrdfits(file,1)
if bb then print, "okdokes"
if n_ext ge 2 then ext2 = mrdfits(file,2)
if n_ext ge 3 then ext3 = mrdfits(file,3) 
if n_ext ge 4 then ext4 = mrdfits(file,4)

;ext3=ext3[where(fs(ext3.ref) eq 'Bautista')]
;uconf=sunique(strmid(ext1.conf,0,8))

nk = n_elements(ext1.ec)
if keyword_set(ev) then ext1.ec = ext1.ec/rtc*rte
if keyword_set(ev) and ip[0] gt 1e4 then ip=ip/rtc*rte

lc   = lonarr(nk) 
lf   = lonarr(nk)
par  = lonarr(nk)
 
for i=0,nk-1 do begin 

   if strlen(fs(ext1[i].term)) eq 1 then begin 
      lc[i]=0
      lf[i]=1 
   endif else begin
      lc[i]=where(L eq strmid(ext1[i].term,strpos(ext1[i].term,fs(ext1[i].mult))+1,1),ic)
      if ic eq 0 then print,'Error in term: '+ext1[i].term
   endelse
   if ext1[i].par eq 'E' then par[i]=0 else par[i]=1

endfor

if n_elements(ip) ge 2 then ext1.mult = ext1.mult-(ext1.ion-1)
;if n_elements(ip) ge 2 then ext1.mult = ext1.mult+(ext1.ion-1)

if keyword_set(cut) then i=where(ext1.ec lt cut,nk) else i=lindgen(nk)
if not keyword_set(maxl) then maxl=max(lc[i])
L = L[0:maxl]

j = where(lc[i] le maxl,nk) & i=i[j]
ext1 = ext1[i]
ind  = i+1
lc   = lc[i]
lf   = lf[i]
par  = par[i]

un=unique(ext1.mult) & un=un[sort(un)]
terml=''
termu=''
for i=0,n_elements(un)-1 do begin
   terml=[terml,'!13^'+fs(un[i])+L]
   termu=[termu,'!13^'+fs(un[i]+1)+L]
   lc[where(ext1.mult eq un[i])]=lc[where(ext1.mult eq un[i])]+n_elements(L)*i
endfor
terml=terml[1:*]
termu=termu[1:*]

maxt  = n_elements(terml)
xr = [-0.5,maxt+0.5]
xtickname = ['!13 ',terml,' ']
xticknamo = ['!13 ',termu,' ']
xtick = maxt+1
if not keyword_set(yr) then yr=[-max(ext1.ec)*0.02,max(ext1.ec)*1.05]
if not keyword_set(ytitle) then ytitle='!13Energy [eV]'

if not keyword_set(ps) then ps='termdiag.ps'

psfull,'l',name=ps
!p.charsize=1.2
!p.thick=1

plot,[0],[0],yr=yr,ys=1,xr=xr,xs=1,xticks=xtick,xtickname=textoidl(xtickname),ytitle=ytitle,xminor=2
;If upper axis is higher ionisation state
if n_elements(ip) ge 2 then axis,xaxis=1,xr=xr,xticks=xtick,xtickname=textoidl(xticknamo),xs=1,xminor=2 

;axis,xaxis=1,xr=xr,xticks=xtick,xtickname=textoidl(space(100L/xtick)+xticknamo),xs=1,xminor=2 ;If upper axis is odd

if not keyword_set(xyouts) then xyouts=''
xyouts,xr[0]+0.8*(xr[1]-xr[0]),yr[0]+0.1*(yr[1]-yr[0]),xyouts

loadct,0
if keyword_set(ext2) and keyword_set(bb) then begin
   
   print, "success"
   ext2.f = alog10(ext2.f*ext1[ext2.irad-1].g)
   ext2 = ext2[sort(ext2.f)]
   ext2.f = ext2.f-min(ext2.f)
   maxf = max(ext2.f)

   if not keyword_set(minf) then minf=-100.

   for i=0,n_elements(ext2)-1 do begin

      il = where(ind eq ext2[i].irad)
      iu = where(ind eq ext2[i].jrad)
      if min([il,iu]) ne -1 then begin 

         if lf[il] eq 1 then x1=lc[il]+maxl/2. else x1=lc[il]+par[il]/2.
         if lf[iu] eq 1 then x2=lc[iu]+maxl/2. else x2=lc[iu]+par[iu]/2.
   
         x=[x1+0.5,x2+0.5]
         y=[ext1[il].ec,ext1[iu].ec]

         if ext2[i].f gt minf and y[0] ge yr[0] and y[1] le yr[1] and x[0] ge xr[0] and x[1] le xr[1] then begin
            oplot,x,y,col=230.-ext2[i].f*230./maxf
         endif
      endif
   endfor
endif

if keyword_set(ext3) and keyword_set(bf) then begin 

   for i=0,n_elements(ext3)-1 do begin

      il = where(ind eq ext3[i].irad)
      iu = where(ind eq ext3[i].jrad)
      if min([il,iu]) ne -1 then begin 

         if lf[il] eq 1 then x1=lc[il]+maxl/2. else x1=lc[il]+par[il]/2.
         if lf[iu] eq 1 then x2=lc[iu]+maxl/2. else x2=lc[iu]+par[iu]/2.
   
         x=[x1+0.5,x2+0.5]
         y=[ext1[il].ec,ext1[iu].ec]

         if y[0] gt yr[0] and y[1] lt yr[1] and x[0] gt xr[0] and x[1] lt xr[1] then begin
            oplot,x,y,col=100.
         endif
      endif
   endfor
endif

if keyword_set(ext4) and keyword_set(cid) then begin 

   ext4 = ext4[where(fs(ext4.cid) eq cid)]
   ext4.rcoll[0] = alog10(ext4.rcoll[0])
   ext4 = ext4[sort(ext4.rcoll[0])]
   ext4.rcoll[0] = ext4.rcoll[0]-min(ext4.rcoll[0])
   maxc = max(ext4.rcoll[0])

   for i=0,n_elements(ext4)-1 do begin

;      if strmid(ext1[ext4[i].icol-1].conf,0,8) ne uconf[conf] then goto,next4
      il = where(ind eq ext4[i].icol)
      iu = where(ind eq ext4[i].jcol)

      if min([il,iu]) ne -1 then begin 

         if lf[il] eq 1 then x1=lc[il]+maxl/2. else x1=lc[il]+par[il]/2.
         if lf[iu] eq 1 then x2=lc[iu]+maxl/2. else x2=lc[iu]+par[iu]/2.

         x=[x1+0.5,x2+0.5]
         y=[ext1[il].ec,ext1[iu].ec]

         if y[0] gt yr[0] and y[1] lt yr[1] and x[0] gt xr[0] and x[1] lt xr[1] then begin
            oplot,x,y,col=0 ; col=230.-ext4[i].rcoll[0]*230./maxc
         endif
      endif
      next4:
   endfor
endif

colors
for i=0,nk-1 do begin
;   if strmid(ext1[i].conf,0,8) ne uconf[conf] and ext1[i].ion ne 2 then goto,next1
   col = 2 
   x = [lc[i]+0.3,lc[i]+0.7]+par[i]/2.
   thi = 4
   if lf[i] eq 1 then begin
      col = 2
      x = [0.1+lc[i]+0.3,lc[i]+0.7+maxl+0.4] 
      thi = 1
   endif
;   oplot,x,[ext1[i].ec,ext1[i].ec],thick=thi,col=col+par[i]*2,linestyle=sty
   if fs(ext1[i].note) eq 'Pr' then col = 2 else col = 4
   if fs(ext1[i].note) eq 'SUPER Ti I' then col = 2
   oplot,x,[ext1[i].ec,ext1[i].ec],thick=thi,col=col
   next1:
endfor

for i=0,n_elements(ip)-1 do print, ip[i]
for i=0,n_elements(ip)-1 do oplot,[-100,100],[ip[i],ip[i]],linestyle=2

psoff

END
