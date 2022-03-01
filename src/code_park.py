

def dis_2ndmoment_forodespy(w,t):
    """
    ODESPY
    This function countains the model for boulder dislodgement. The equation of
    motion is adjusted to the Tower program. Fully and partially submerged
    boulders are possible. No functions are called within this functions.
    """
    x1,y1 = w
    aa,bb,cc,alpha,uu, rhos, cd, cl, mu, eta, CmA, x0,y0 = p
    rho = 1032.0
    gra = 9.81
    bx = [0.0,bb,bb,0.0,0.0]
    by = [0.0,0.0,cc,cc,0.0]
    water_line_x = [-10,10]
    water_line_y = [eta, eta]
#----------------Rotating Boulder--------------------------
    newbx = zeros(len(bx))
    newby = zeros(len(by))
    angle = x1 - arctan(bb/cc)
    if angle == 0.0: angle=0.001
    for i in range(len(bx)):
        newbx[i] = bx[i]*cos(angle) - by[i] * sin(angle)
        newby[i] = bx[i]*sin(angle) + by[i] * cos(angle)  
##    newbx,newby = rotate_boulder(bx,by,x1-arctan(bb/cc))
#------------------------------------------------------------

#---------------COM full body---------------------------------
    com_f_x,com_f_y = center_of_mass_xy(newbx,newby)
    mi_f = mi_polygon(newbx,newby)
#------------------------------------------------------------

#-------------COM of submerged part of body------------------
##-------Calculate coordinates of submerged part first--------------------
    adjbx,adjby = find_new_cords(newbx,newby,water_line_x,water_line_y)
    mi_a = mi_polygon(adjbx,adjby)
##------Calculate (x,) coordinates of submerged part-------------------
    com_a_x,com_a_y = center_of_mass_xy(adjbx,adjby)
#----------------------------------------------------------------------------

#----------------Body volumes-------------------------------------------
    area_f = find_area(newbx,newby)
    area_adj = find_area(adjbx,adjby)
    vol_full = aa * area_f
    vol_adj = aa * area_adj 


    ad_x_min = amin(adjbx); ad_x_max = amax(adjbx)
    ad_y_min = amin(adjby); ad_y_max = amax(adjby)


    ld_v = com_a_y
    ll_h = com_a_x
    lb_h = com_a_x
    lw_h = com_f_x
    drag_c = float(abs(ad_y_min-ad_y_max))
    lift_b = float(abs(ad_x_min-ad_x_max))
    m = (mi_f * rhos + CmA * mi_a * rho)**(-1.0)
    fw=rhos*gra*vol_full
    fb=rho*gra*vol_adj
    fd=0.5*aa*drag_c*cd*rho*uu**2.0
    fl=0.5*lift_b*aa*cl*rho*uu**2.0
    fsum2=fd*ld_v+fl*ll_h+fb*lb_h-fw*lw_h
    term1 = fsum2 * m
    t1=t2=0.0
    if fd*ld_v+fl*ll_h+fb*lb_h-fw*lb_h >0.0:
        t1 = y1
        t2 =term1
    return t1,t2



def dislodgement_forodespy(w,t):
    """
    This function countains the model for boulder dislodgement and uses odespy to 
    solve the governing equations. The equation of motion is adjusted to the Tower
    program. Fully and partially submerged boulders are possible. No functions are
    called within this functions.
    """
    x1,y1 = w
    aa,bb,cc,alpha,uu, rhos, cd, cl, mu, eta, CmA, x0,y0 = p
    rho = 1032.0
    gra = 9.81
    bx = [0.0,bb,bb,0.0,0.0]
    by = [0.0,0.0,cc,cc,0.0]
    water_line_x = [-10,10]
    water_line_y = [eta, eta]
#----------------Rotating Boulder--------------------------
    newbx = zeros(len(bx))
    newby = zeros(len(by))
    angle = x1 - arctan(bb/cc)
    if angle == 0.0: angle=0.001
    for i in range(len(bx)):
        newbx[i] = bx[i]*cos(angle) - by[i] * sin(angle)
        newby[i] = bx[i]*sin(angle) + by[i] * cos(angle)  
##    newbx,newby = rotate_boulder(bx,by,x1-arctan(bb/cc))
#------------------------------------------------------------

#---------------COM full body---------------------------------
    com_f_x,com_f_y = center_of_mass_xy(newbx,newby)
#------------------------------------------------------------

#-------------COM of submerged part of body------------------
##-------Calculate coordinates of submerged part first--------------------
    adjbx,adjby = find_new_cords(newbx,newby,water_line_x,water_line_y)
##------Calculate (x,) coordinates of submerged part-------------------
    com_a_x,com_a_y = center_of_mass_xy(adjbx,adjby)
#----------------------------------------------------------------------------

#----------------Body volumes-------------------------------------------
    area_f = find_area(newbx,newby)
    area_adj = find_area(adjbx,adjby)
    vol_full = aa * area_f
    vol_adj = aa * area_adj 


    ad_x_min = amin(adjbx); ad_x_max = amax(adjbx)
    ad_y_min = amin(adjby); ad_y_max = amax(adjby)


    ld_v1 = com_a_y
    ll_h1 = com_a_x
    lb_h1 = com_a_x
    lw_h1 = com_f_x
    drag_c1 = float(abs(ad_y_min-ad_y_max))
    lift_b1 = float(abs(ad_x_min-ad_x_max))
    ld_v = ld_v1
    ll_h=ll_h1
    lb_h=lb_h1
    lw_h =lw_h1
    drag_c = drag_c1
    lift_b=lift_b1
    ll_f = (bb**2.0 + cc**2.0)
    if eta < 0.6*ad_y_max:
        ll_adj = (eta**2.0 + (vol_adj/eta)**2.0)
    else:
        ll_adj = ll_f
    m= ((1./3. * ll_f*area_f*aa *rhos + ll_adj *CmA*rho* vol_adj*aa))**(-1.0)
    fw=rhos*gra*vol_full
    fb=rho*gra*vol_adj
    fd=0.5*aa*drag_c*cd*rho*uu**2.0
    fl=0.5*lift_b*aa*cl*rho*uu**2.0 #JI, not currently using lift, TODO if lift added: confirm bc is characteristic area and not ab
    fsum2=fd*ld_v+fl*ll_h+fb*lb_h-fw*lw_h
    term1 = fsum2 * m
    t1=t2=0.0
    if fd*ld_v+fl*ll_h+fb*lb_h-fw*lb_h >0.0:
        t1 = y1
        t2 =term1
    return t1,t2

def sliding_forodespy(w,t):
    """
    This function countains the model for boulder dislodgement. The equation of
    motion is adjusted to the Tower program. Fully and partially submerged
    boulders are possible. No functions are called within this functions.
    """
    x1,y1 = w
    aa,bb,cc,alpha,uu, rhos, cd, cl, mu, eta, CmA, x0,y0 = p
    rho = 1032.0
    gra = 9.81
    bx = [0.0,bb,bb,0.0,0.0]
    by = [0.0,0.0,cc,cc,0.0]
    water_line_x = [-10,10]
    water_line_y = [eta, eta]
# #----------------Rotating Boulder--------------------------
    newbx = zeros(len(bx))
    newby = zeros(len(by))
    angle = x1 - arctan(bb/cc)
    angle=0.001
    for i in range(len(bx)):
        newbx[i] = bx[i]*cos(angle) - by[i] * sin(angle)
        newby[i] = bx[i]*sin(angle) + by[i] * cos(angle)  
##    newbx,newby = rotate_boulder(bx,by,x1-arctan(bb/cc))
# #------------------------------------------------------------
# 
#---------------COM full body---------------------------------
    com_f_x,com_f_y = center_of_mass_xy(newbx,newby)
#------------------------------------------------------------

#-------------COM of submerged part of body------------------
##-------Calculate coordinates of submerged part first--------------------
    adjbx,adjby = find_new_cords(newbx,newby,water_line_x,water_line_y)
##------Calculate (x,) coordinates of submerged part-------------------
    com_a_x,com_a_y = center_of_mass_xy(adjbx,adjby)
#----------------------------------------------------------------------------

#----------------Body volumes-------------------------------------------
    area_f = find_area(newbx,newby)
    area_adj = find_area(adjbx,adjby)
    vol_full = aa * area_f
    vol_adj = aa * area_adj 


    ad_x_min = amin(adjbx); ad_x_max = amax(adjbx)
    ad_y_min = amin(adjby); ad_y_max = amax(adjby)


    ld_v1 = com_a_y
    ll_h1 = com_a_x
    lb_h1 = com_a_x
    lw_h1 = com_f_x
    drag_c1 = float(abs(ad_y_min-ad_y_max))
    lift_b1 = float(abs(ad_x_min-ad_x_max))
    ld_v = ld_v1
    ll_h=ll_h1
    lb_h=lb_h1
    lw_h =lw_h1
    drag_c = drag_c1
    lift_b=lift_b1
    ll_f = (bb**2.0 + cc**2.0)
    if eta < 0.6*ad_y_max:
        ll_adj = (eta**2.0 + (vol_adj/eta)**2.0)
    else:
        ll_adj = ll_f
    m = ((1./3.*area_f*aa *rhos + CmA*rho* vol_adj*aa))**(-1.0)
    fw = rhos * gra * vol_full
    fb = rho * gra * vol_adj
    fd = 0.5 * aa * drag_c * cd * rho * (uu-y1)**2.0
    fsum2=fd*ld_v + mu * (fb*lb_h - fw*lw_h)
    term1 = fsum2 * m
    t1=t2=0.0
    if fd*ld_v + mu *(fb*lb_h -fw*lb_h) >0.0:
        t1 = y1
        t2 = term1
    return t1,t2
    
    
def terminated_rot(u,t,step_no):
    return False if u[step_no,0]<=1.05*thetacr else True

def terminated_sli(u,t,step_no):
    return False if u[step_no,0]<=0.5*bbb_cp else True

def run_depth_withodespy(a,b,c,froude,eta_v,rhos,cd,cl,mu,CmA):
    """ 
    Runnning the Monte-Carlo type simulations for different flow depths and
    Froude numbers. To solve the equations, odespy is emplpoyed.
    """
    from numpy import linspace,shape,arctan,zeros,pi,sqrt
    import odespy
    #rhos = float(rhos)
    roughness = 0.0 #JI set for Towers
    stoptime =10. #JI see Weiss & Diplas
    nt=100
    x1 = [-b,-b,0.0]
    y1 = [0.0,c,c]
    t=linspace(0,stoptime,nt)
    alpha = arctan(b/c) #JI 0.0 * pi/180.0
    global thetacr
    thetacr = 0.5*pi #JI 0.5*pi - arctan(c/b) - alpha
    global p
    p = [a,b,c,alpha,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA,x1,y1]
#----------------Rotation---------------------------------------------------------
    w0_rot = [arctan(b/c),0.0]
    abserr = 1.0e-8
    relerr = 1.0e-6
    #solver=odespy.RKFehlberg(dislodgement_forodespy)
    #solver=odespy.RK4(dislodgement_forodespy)
    solver=odespy.RK4(dis_2ndmoment_forodespy)
    solver.set_initial_condition(w0_rot)
    wsol_rot,ti_rot = solver.solve(t,terminated_rot)
    theta = zeros(len(wsol_rot[:,0]))
    theta[:] = wsol_rot[:,0]
#---------------------------------------------------------------------------------

#----------------Sliding----------------------------------------------------------
    w0_sli = [0.0,0.0]
    abserr = 1.0e-8
    relerr = 1.0e-6
    global bbb_cp
    bbb_cp = b
    #solver=odespy.RKFehlberg(dislodgement_forodespy)
    solver=odespy.RK4(sliding_forodespy)
    solver.set_initial_condition(w0_sli)
    wsol_sli,ti_sli = solver.solve(t,terminated_sli)
    sli_dist = zeros(len(wsol_sli[:,0]))
    sli_dist[:] = wsol_sli[:,0]
#----------------------------------------------------------------------------------
    return ti_rot[-1],theta[-1],ti_sli[-1],sli_dist[-1]

def adjust_block(polygon_body, flow_d):
#    from numpy import array,amin,amax
#    from shapely import wkt
#    from shapely.ops import linemerge, unary_union, polygonize
#    from shapely.geometry.polygon import LinearRing, Polygon

    #polygon1 = Polygon([(x_local[0], y_local[0]), (x_local[1], y_local[1]), (x_local[2], y_local[2]), (x_local[3], y_local[3]),(x_local[4], y_local[4])])
    polygon2 = LinearRing([(-10,flow_d),(0.0,flow_d),(10,flow_d)])

    merged = linemerge([polygon_body.boundary, polygon2])
    borders = unary_union(merged)
    polygons = polygonize(borders)
    #print(borders)
    list_poly = list(polygons)
    cent_x,cent_y = list_poly[0].centroid.coords.xy
    return_x, return_y = list_poly[0].exterior.coords.xy
    #print(y)
    return_x = array(return_x)
    return_y = array(return_y)
    cent_x = array(cent_x)
    cent_y = array(cent_y)
    #print(y)
    return amin(return_x),amax(return_x),amin(return_y),amax(return_y),float(cent_x[0]),float(cent_y[0])

def polygone2xy(polygone_local):
    return_x, return_y = polygone_local.exterior.coords.xy
    return return_x, return_y



def center_of_mass(polygone_in):
#    from numpy import array
#    from shapely.geometry.polygon import LinearRing, Polygon
#from shapely import *
    centr_x,centr_y = polygone_in.centroid.coords.xy
    centr_x = array(centr_x)
    centr_y = array(centr_y)
    return float(centr_x[0]),float(centr_y[0])
    
def ode45_step(f,x_l,t_l,dt_l,*args):
    print(x_l,t_l,dt_l)
    k = dt_l
    k1 = k * f(t_l,x_l,*args)
    k2 = k * f(t_l + 0.5*k, x_l + 0.5*k1, *args)
    k3 = k * f(t_l + 0.5*k, x_l + 0.5*k2, *args)
    k4 = k * f(t_l + dt_l, x_l + k3, *args)
    return x + 1./6. * (k1 + 2.0*k2 + 2.0*k3 + k4)

def ode45(f, t_l, x0_l, *args):
    n = len(t_l)
    x = zeros([n,len(x0_l)])
    x[0] = x0_l
    print(x0_l)
    print(n)
    for i in range(n-1):
        dt_l = t_l[i+1] - t_l[i]
        x[i+1] = ode45_step(f, x[i], t_l[i], dt_l, *args)
    return t, x


def dislodgement_rw(w,t,p): #old functiont that uses shapely polygons
    """
    This function countains the model for boulder dislodgement. The equation of
    motion is adjusted to the Tower program. It is assumed that the tower is
    completely submerged in water during the flip.
    """
#    from numpy import sqrt, pi,sin,cos,arctan,zeros,abs,amax,amin
#    from shapely.geometry.polygon import LinearRing, Polygon
    x1,y1 = w
    a,b,c,alpha,uu, rhos, cd, cl, mu, eta, CmA, x0,y0, thetacrit = p
    rho = 1032.0
    gra = 9.81
    bx = [0.0,b,b,0.0,0.0]
    by = [0.0,0.0,c,c,0.0]
#     newbx = zeros(len(bx))
#     newby = zeros(len(bx))
#     angle = x1 - arctan(b/c)
#     for i in range(len(by)):
#         newbx[i] = bx[i]*cos(angle) - by[i] * sin(angle)
#         newby[i] = bx[i]*sin(angle) + by[i] * cos(angle)  

    newbx,newby = rotate_boulder(bx,by,x1-arctan(b/c))
# 
#     #print('in')
    polygon1 = Polygon([(newbx[0], newby[0]), (newbx[1], newby[1]), (newbx[2],newby[2]), (newbx[3], newby[3]),(newbx[4], newby[4])])
    com_f_x,com_f_y = center_of_mass(polygon1)
#    adjust_block(polygon1,eta)
    ad_x_min,ad_x_max,ad_y_min,ad_y_max,com_a_x,com_a_y = adjust_block(polygon1,eta)
#     ad_x,ad_y = polygone2xy(polygon_adj)
#    com_a_x,com_a_y = center_of_mass(polygon_adj)
    ld_v1 = com_a_y    
    ll_h1 = com_a_x
    lb_h1 = com_a_x
    lw_h1 = com_f_y
    drag_c1 = float(abs(ad_y_min-ad_y_max))
    lift_b1 = float(abs(ad_x_min-ad_x_max))
#    print(type(lb_h))
    #print('above ',ld_v1,ll_h1,lb_h1,lw_h1,drag_c1,lift_b1)
    #print(drag_c1)
    #print(lift_b)
    	#print(alpha,x1,x1-alpha)
    #print(len(yP))
    ld_v = ld_v1
    ll_h=ll_h1
    lb_h=lb_h1
    lw_h =lw_h1
    drag_c = drag_c1
    lift_b=lift_b1
 #   print(type(ld_v))
#(3.0, 3.0, 1.5, 1.5, 6.0, 6.0)
#     if eta >= h_fd:
#         new_c = h_fd
#new_b = w_fl
    #print(xP)
    ll = (b**2.0 + c**2.0)
    lo = sqrt(0.25*(b**2.0)+(c**2.0))
#     alpha = 0.0
    lv = lo * cos(x1-arctan(b/c))
    lh = lo * sin(x1-arctan(b/c))
#    print(x1-arctan(b/c),thetacrit)
    #m = (ll* (1./3.0*rhos+rho*CmA)*a*b*c)**(-1.0) #JI to update equation. added_mass_coef~0.7
    m= (ll * (1./3. * a * b *c *rhos + CmA*rho* a*b* drag_c))**(-1.0)
    fw=rhos*gra*a*b*drag_c
    fb=rho*gra*a*b*lift_b
    fd=cd*rho*uu**2.0*0.5*(a*drag_c)
    fl=cl*rho*uu**2.0*0.5*(lift_b*a) #JI, not currently using lift, TODO if lift added: confirm bc is characteristic area and not ab
    fsum2=fd*ld_v+fl*ll_h+fb*lb_h-fw*lw_h
    term1 = fsum2 * m
#    print('down ',type(m),type(fw),type(fd),type(fl),type(fb),type(term1),type(ll),type(CmA))
    if fd*ld_v+fl*ll_h+fb*lb_h-fw*lb_h >0.0 and x1<thetacrit:
#    if fd +fl +fb -fw >0.0: 
        f = [y1, term1]
        return f
    else:
#    if fd*ld_v+fl*ll_h+fb*lb_h-fw*lb_h > 0.0:
        f = [0.0,0.0]
        return f
#    if x1 > thetacrit:
#        return [0.0,0.0]


def run_depth_rw(froude,eta_v,a,b,c,rhos,cd,cl,mu,CmA): #old code that uses scipy.odeint
    """
    Runnning the Monte-Carlo type simulations for different flow depths and
    Froude numbers.
    """
#    from numpy import linspace,shape,arctan,zeros,pi,sqrt
#    from scipy.integrate import odeint
#    import odespy
    #rhos = float(rhos)
    roughness = 0.0 #JI set for Towers
    stoptime =5.5 #JI see Weiss & Diplas
    nt=500
    x1 = [-b,-b,0.0]
    y1 = [0.0,c,c]
    t=linspace(0,stoptime,nt)
    alpha = arctan(b/c) #JI 0.0 * pi/180.0
    thetacr = 0.5*pi  #JI 0.5*pi - arctan(c/b) - alpha
    global p
    p = [a,b,c,alpha,froude*sqrt(9.81*eta_v),rhos,cd,cl,mu,eta_v,CmA,x1,y1,thetacr]

    w0 = [arctan(b/c),0.0]
#    w0 = [0.5,0.5]
    abserr = 1.0e-8
    relerr = 1.0e-6
#    global b_g
#    b_g = b
    #wsol = odeint(dislodgement, w0,t, args=(p,), atol = abserr, rtol = relerr)
    wsol = odeint(dislodgement_rw, w0,t, args=(p,))
#    solver=odespy.Fehlberg(dislodgement_rw)
#    time_points=linspace(0,endtime,n+1)
#    solver.set_initial_condition(w0)
#    wsol = solver.solve(t)#:,terminated) 
    
    
    theta = zeros(len(wsol[:,0]))
    theta[:] = wsol[:,0]

##JI CHECK THIS
#	print eta_v, froude, max(theta), min(theta), thetacr #JI check
# Robert commented this out. I think we should set a debug flag.
# 	theta[theta<thetacr] = 'NaN'
# 	if max(theta)> 0.0:
    if max(theta)> thetacr:
        return 1
    else:
        return 0


##JI initial code to calculate moments of inertia, etc. No longer used
##Calculate submerged volume and distance from rotation axis for block
def submergedArea(a,b,c,lo,lW,eta,th,thetacr):
	from numpy import sin, cos, tan, pi
	x2,z2,x3,z3,x4,z4 = blockcorners(b,c,lo,th,thetacr)
# 	Areasub = -1.
# 	lsub = 0.
	if z2 < eta and z4 >= eta:
		#h = b of block
		#b = eta/cos(th)
		#a = eta/cos(th) - b*tan(th) where b is b of block
		#c = 0.
		Areasub,cxrelc,cyrelc = trapezoidchar(eta/cos(th) - b*tan(th), eta/cos(th), 0., b)
		cx,cy = rotatecoord(cxrelc, cyrelc, th - pi/2.)
		lsub = -cx		
# 		print th,Areasub,cxrelc,cyrelc,lsub #JItest		
	elif z2 >= eta and z4 < eta:
		#h = c of block
		#b = eta/sin(th)
		#a = eta/sin(th) - c/tan(th) where c is c of block
		#c = 0.
		Areasub,cxrelc,cyrelc = trapezoidchar(eta/sin(th) - c/tan(th), eta/sin(th), 0., c)	
		cx,cy = rotatecoord(-cxrelc, cyrelc, th)
		lsub = -cx
	elif z2 < eta and z4 < eta:
		trib = (z3 - eta)/sin(th)
		trih =  (z3 - eta)/cos(th)
		Atri = 0.5*trib*trih
		cx,cy = rotatecoord(-trib/3.,-trih/3., th - pi/2.)
		ltri = -(x3 + cx)
		Areasub = b*c - Atri
		lsub = (b*c*lW - Atri*ltri)/Areasub
	else:
		trib = eta/cos(th)
		trih = eta/sin(th)
		Areasub = 0.5*trib*trih
		cx,cy = rotatecoord(trib/3.,trih/3., th - pi/2.)
		lsub = -cx
	Vsub = a*Areasub
	return Vsub,lsub 


##THIS IS AN OLD VERSION, DELETE? JI version for Towers
def dislodgement(w,t,p):
	"""
	This function countains the model for block overturning to investigate
	Caesarea Tower overturning. The function allows for partially submerged
	and completely submerged blocks.
	"""
	from numpy import sqrt, pi, sin, cos, min
	th,dthdt = w
	a,b,c,alpha,rough ,uu, rhos, cd, cl, mu, eta, CmA, thetacr = p
	rho = 1032.0
	gra = 9.81
	lo = sqrt(0.25*(c**2.+b**2.))
	btil = b*cos(th) + c*sin(th)
	ctil = min([c,eta])
	ceq = c #placeholder
	beq = b #placeholder
	#JI TODO: estimate ceq,beq from Vsub,b,c
	lW = lo*sin(thetacr-th) #W acts at the block's center of mass
	if eta < c: #emergent
		lD = 0.5*eta #FD should act at the center of pressure
		#JI TODO: confirm lD
		Vsub,lB = submergedArea(a,b,c,lo,lW,eta,th,thetacr)
	else: #submerged
		lD = lo*cos(thetacr-th) #FD should act at the center of pressure
		Vsub = a*b*c
		lB = lW
		#JI TODO: confirm lD
# 	lB = lW #JI placeholder
	lL = lW #JI placeholder. FL acts at center of pressure
	#JI TODO: confirm lL
# 	Vsub = a*b*ctil #placeholder
	#JI TODO: correct Vsub calculation
	Wt = rhos*a*b*c
# 	print rhos #JItest

# JItest - turns off Buoyancy fix
# 	Vsub = a*b*min([c,eta])
# 	lB = lW
# endJItest

	FB = rho*Vsub
	FL = 0.5*cl*rho*(uu**2.)*a*btil #CONFIRM PLATFORM AREA IS PROJECTED AREA IN LINE WITH FLOW (horizontal plane)
	FD = 0.5*cd*rho*(uu**2.)*a*(2*lD) #PLATFORM AREA IS PROJECTED AREA INTO FLOW (vertical plane)
	IA = Wt*(c**2.+b**2.)/3.
	IM = CmA*FB*(ceq**2.+beq**2.)/3.
	d2thdt2 = ( -Wt*lW + FB*lB + FL*lL + FD*lD ) / (IA + IM) #confirm term is negative, with CW sign convention for theta
	sumM0 = -Wt*lW + FB*lB + FL*lL + FD*lD  #sum of moments	
# 	print eta, uu/sqrt(gra*eta), th, Wt, lW, FB, lB, FL, lL, FD, lD#JItest	
	if sumM0 > 0.:
		f = [dthdt,d2thdt2]
	else:
		f = [0.,0.]
	return f
#-----------------------------End Jen---------------------------------------------------------
