	common/block2/list,head,map
	integer m,ncell,mapsiz
	parameter(m=15,ncell=m*m) !2<m<=3*NU/10, NEEDS TO BE CHANGED WHEN NU IS CHANGED
        parameter(mapsiz=4*ncell)
        integer list(n),head(ncell),map(mapsiz)

