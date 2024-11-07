FUNCTION unique,v,nr=nr

u=v[0]

for i=1,n_elements(v)-1 do begin
   j=where(u eq v[i])
   if j[0] eq -1 then u=[u,v[i]]
endfor

nr=n_elements(u)
return,u
END
