from grid import*
from math import*
import cmath as cma
import scipy as sipy
#from predictionCorrection import*

def concaY(A,eps,y):
	adj_A = A.conj().T
	temp_A = -y*eye(len(A)) + complex(0,1)*adj_A
	temp_Abis = complex(0,1)*A +y*eye(len(A))
	epsI = eye(len(A))*eps
	res1 = np.concatenate((temp_A,epsI),axis = 0)
	res2 = np.concatenate((-epsI,temp_Abis),axis = 0)
	res = np.concatenate((res1,res2),axis = 1)
	return res


def concaX(A,eps,x):
	adj_A = A.conj().T
	temp_A = x*eye(len(A)) - adj_A
	temp_Abis = A -x*eye(len(A))
	epsI = eye(len(A))*eps
	res1 = np.concatenate((temp_A,epsI),axis = 0)
	res2 = np.concatenate((-epsI,temp_Abis),axis = 0)
	res = np.concatenate((res1,res2),axis = 1)
	return res

def concaTETA(A,eps,teta):
	i = complex(0,1)
	adj_A = A.conj().T
	temp_A = i*cma.exp(i*teta)*adj_A
	temp_Abis = A*i*cma.exp(-i*teta)
	epsI = eye(len(A))*eps
	res1 = np.concatenate((temp_A,epsI),axis = 0)
	res2 = np.concatenate((-epsI,temp_Abis),axis = 0)
	res = np.concatenate((res1,res2),axis = 1)
	return res

def concaG(A,eps,r):
	epsI = eye(len(A))*eps
	rI = eye(len(A))*r
	res1 = np.concatenate((-epsI,rI),axis = 0)
	res2 = np.concatenate((A,zeros((len(A),len(A)))),axis = 0)
	res = np.concatenate((res1,res2),axis = 1)
	return res 

def concaD(A,eps,r):
	adj_A = A.conj().T
	epsI = eye(len(A))*eps
	rI = eye(len(A))*r
	res1 = np.concatenate((zeros((len(A),len(A))),adj_A),axis = 0)
	res2 = np.concatenate((rI,-epsI),axis = 0)
	res = np.concatenate((res1,res2),axis = 1)
	return res 



def intersect(A,eps,r):
	gauche = concaG(A,eps,r)
	droite = concaD(A,eps,r)
	res = sipy.linalg.eig(gauche,droite,right =False)

	return [a for a in res if (abs(abs(a)-1)<0.00001)]



def largestPureIm(M,A,eps,y):
	i = complex(0,1)
	n = len(A)
	val_p ,w= linalg.eig(M)
	res = val_p[0]
	for v in val_p:
		
		if(abs(v.real)<0.001) and res.imag<v.imag:
			
			res = v
	return res


def surchVert(M,A,eps,x):
	
	i = complex(0,1)
	val_p ,w= linalg.eig(M)
	res = []
	n = len(A)
	for v in val_p:
		
		if(abs(v.real)<0.00001):
			list_sig = svdvals((x+v.imag*i)*eye(n)-A)
			minsig = list_sig[len(list_sig)-1]
			if((abs(minsig-eps))<0.0001):	
				res.append(v.imag)
	return res

def surchTeta(A,eps,teta):
	M = concaTETA(A,eps,teta)
	i = complex(0,1)
	n = len(A)
	val_p ,w= linalg.eig(M)
	res = val_p[0]
	for v in val_p:
		
		if(abs(v.real)<0.000001) and res.imag<v.imag and v.imag>=0:
			res = v
	return res.imag*cma.exp(i*teta)



def abscissae(M,eps,axi = None):
	"""
	 Calcule le pseudo abscisse du eps pseudospectre d'une matrice M
 	"""
	#step1
	A = np.array(M)
	n = len(A)
	val_p,vecteurs=linalg.eig(M)
	max_valp = val_p[0]
	list_point_temp = []

	for val in val_p :
		if val.real > max_valp.real :
			max_valp = val


	#step 2
	maxtemp = -1
	point_prec = complex(0,0)
	Maty = concaY(A,eps,max_valp.imag)
	point_cour =  complex(largestPureIm(Maty,A,eps,max_valp.imag).imag,max_valp.imag)
	cpt = 0
	plt.plot(point_cour.real,point_cour.imag,"b:.")
	
	while(abs(point_prec-point_cour)>0.000001 and cpt<100):
	
		Matx = concaX(A,eps,point_cour.real)
		pointempY_list = []
		pointempY_list = surchVert(Matx,A,eps,point_cour.real)
		if(len(pointempY_list)%2==0 and len(pointempY_list)!=0):
			j = 0
			h = pointempY_list[j]
			pointempY_list = sorted(pointempY_list )
			while(j<len(pointempY_list)):
				x = np.array([point_cour.real,point_cour.real])
				y = np.array([pointempY_list[j],pointempY_list[j+1]])
				plt.plot(x,y,color = "black")
				plt.plot(point_cour.real,pointempY_list[j],"y:.")
				plt.plot(point_cour.real,pointempY_list[j+1],"y:.")
				new_point = complex(point_cour.real,((max(pointempY_list[j],pointempY_list[j+1]))+(min(pointempY_list[j],pointempY_list[j+1])))/2)
				list_point_temp.append(new_point)
				j = j+2
			maxtemp = -1000
			#print(list_point_temp)
			for elt in list_point_temp :
				Maty = concaY(A,eps,elt.imag)
				absi = largestPureIm(Maty,A,eps,elt.imag).imag
				if(absi.imag>maxtemp):
					maxtemp = absi
					h = elt.imag
					x_aff = elt.real
			list_point_temp = []
			point_prec = point_cour
			point_cour = complex(maxtemp,h)
			x = np.array([maxtemp,x_aff])
			y = np.array([h,h])
			plt.plot(x,y,color = "black")
			cpt = cpt+1
		else :
			print("fin")
			j = 0
			while(j<len(pointempY_list)):
				j = j+1
			break
	#pointemp = complex(pointempX.imag,max_valp.imag)
	plt.plot(point_cour.real,point_cour.imag,"r:.",label ="abscisse max")
	#print(point_cour)
	if(axi != None):
		axi.legend()
	return point_cour
	


	
def graphe_abs(M,eps_min,eps_max,pas):
	y = []
	x = np.linspace(eps_min,eps_max,pas)
	i = 0
	for elt in x :
		print(i)
		y.append(abscissae(M,elt).real)
		i = i+1
	plt.plot(x, y,color = "c")
	plt.show()


def radii(M,eps,axi = None):
	"""
	 Calcule le pseudo rayon du eps pseudospectre d'une matrice M
 	"""
	#step1
	i = complex(0,1)
	A = np.array(M)
	n = len(A)
	val_p,vecteurs=linalg.eig(M)
	max_valp = val_p[0]
	list_point_temp = []

	for val in val_p :
		if abs(val) > abs(max_valp) :
			max_valp = val
	r,phi =cma.polar(max_valp)
	point_temp = surchTeta(A,eps,phi)
	#plt.plot(point_temp.real,point_temp.imag,"g:.")

	cpt = 0
	point_prec = complex(10000,10000)
	while(abs(point_temp-point_prec)>0.0000001 and cpt<100):
		point_prec = point_temp
		r,phi = cma.polar(point_temp)
		list_point = intersect(A,eps,r)
		list_tri = []
		for elt in list_point : 
			rtemp,phitemp = cma.polar(elt)
			list_tri.append(phitemp)
		list_tri = sorted(list_tri)
		list_point = [cma.rect(r,elt) for elt in list_tri]
		j = 0
		while(j<len(list_point)):
				rbis1,phi1 = cma.polar(list_point[j])
				point1 = cma.rect(r,phi1)	
				rbis2,phi2 = cma.polar(list_point[j+1])
				point2 = cma.rect(r,phi2)
				plt.plot(point1.real,point1.imag,"r:.")
				plt.plot(point2.real,point2.imag,"r:.")
				point_mid = complex((point1.real+point2.real)/2,(point1.imag+point2.imag)/2)
				rbis,phi = cma.polar(point_mid)
				point_candidat = surchTeta(A,eps,phi)
				plt.plot(point_candidat.real,point_candidat.imag,"b:.")
				if abs(point_candidat)>abs(point_temp):
					point_temp = point_candidat
				j = j+2
		plt.plot(point_temp.real,point_temp.imag,"k:.")
		cpt = cpt+1
	plt.plot(point_temp.real,point_temp.imag,"m:.",label ="rayon max r1")
	c2 = (plt.Circle((0,0),radius =abs(point_temp),color ='r',fill=False,label="C(0,r1)"))
	plt.gca().add_patch(c2)
	if(axi != None):
		axi.legend()
	return point_temp


def graphe_radii(M,eps_min,eps_max,pas):
	y = []
	x = np.linspace(eps_min,eps_max,pas)
	i = 0
	for elt in x :
		print(i)
		y.append(abs(radii(M,elt)))
		i = i+1
	plt.plot(x, y,color = "b")
	plt.show()