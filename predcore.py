
from grid import*
import matplotlib as mp
seed()


def limite(M,eps):

	n = len(M)
	borne_X = [0,0]
	borne_Y = [0,0]
	x = 0
	y = 0
	r = 0
	cercle = ((x,y),r)
	liste_cercle = []
	liste_c = []
	for i in range(n):
		x = M[i][i].real
		y = M[i][i].imag
		r = 0
		for j in range(n):
			if(i != j):
				r = r+norm(M[i][j])
		cercle = ((x,y),r)
		liste_c.append(cercle)
		r = r+sqrt(n)*eps
		cercle = ((x,y),r)
		liste_cercle.append(cercle)
		if(borne_X[0] > x-r):
			borne_X[0] = x-r
		if(borne_X[1] < x+r):
			borne_X[1] = x+r
		if(borne_Y[0] > y-r):
			borne_Y[0] = y-r
		if(borne_Y[1] < y+r):
			borne_Y[1] = y+r
	return (borne_X,borne_Y,liste_cercle,liste_c)

	
def cherche_z(M,val_p,eps,tol,n):
	"""
	fonction qui permet de trouver un point sur le pseudospectre à partir d'une valeur propre
	"""
	z = val_p+complex(0,1)*eps
	#svdvals((z*eye(n))-M)[n-1]
	#sig_min[n-1]
	u_min,sig_min,v_min=svd((z*eye(n))-M)
	i = 0
	while (abs((sig_min[n-1]-eps)/eps)>tol and i<1000):
		u_min,sig_min,v_min= svd((z*eye(n))-M)
		z=z-((sig_min[n-1]-eps)/(complex(0,-1)*(ny.dot((v_min[n-1,:]),(u_min[:,n-1])))).real)*complex(0,1)
		i= i+1
	return z


def step_2_3(z1,M,taille_pred,n,eps,val_p):
	"""
	effectue une itération de prédiction puis de correction
	"""
	u_min,sig_min,v_min = svd((z1*eye(n))-M)
	r = (complex(0,1)*(ny.dot((v_min[n-1,:]),(u_min[:,n-1]))))/(abs(ny.dot((v_min[n-1,:]),(u_min[:,n-1]))))
	tau = min(taille_pred,norm(z1-val_p)/2)
	z = z1 +r*tau

	u_min,sig_min,v_min= svd((z*eye(n))-M)
	z_res = z-(sig_min[n-1]-eps)/(ny.dot(ny.conjugate(u_min[:,n-1]),ny.conjugate(v_min[n-1,:])))
	return z_res

def predcore(M,list_eps,tol,taille_pred):
	"""
	implémentation de l'algorithme de prédiction-correction, avec gestion de la barre de chargement
	"""
	lc = ['cyan','blue','m','g','r','k']
	n = len(M) 
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.gca()
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')
	window = Tk()
	window.title("processing")
	window.minsize(480,160)
	window.maxsize(640,360)
	label_title1 = Label(window,text= "Calcul en cours :",font ='Times')
	label_title1.pack(expand = YES)
	label_subtitle = Label(window,text= "nombre de valeur propre parcouru : 0/"+str(len(val_p)),font ='Times')
	label_subtitle.pack(expand = YES)
	label_subtitle2 = Label(window,text = "nombre de prédiction/correction : 0/1000",font = 'Times')
	label_subtitle2.pack(expand = YES)
	can = Canvas(window)
	can.pack(expand=YES)
	can.create_rectangle(89,39,291,61,outline="black")
	c1 = can.create_rectangle(90,40,90,60,outline="lime green",fill="lime green")
	i = 0
	cpt = 0
	k= 0

	for eps in list_eps :

		for val in val_p :
			if(eps == list_eps[0]):
				plt.plot(val.real,val.imag,'k:.',)
			z = cherche_z(M,val,eps,tol,n)
			b = True
			z_temp = step_2_3(z,M,taille_pred,n,eps,val)
			lx = []
			ly = []
			cpt = 1
			while b and cpt < 1000:
				z_temp = step_2_3(z_temp,M,taille_pred,n,eps,val)
				lx.append(z_temp.real)
				ly.append(z_temp.imag)
				b = norm(z-z_temp)>(taille_pred*0.3)
				cpt = cpt+1
				label_subtitle2.config(text ="nombre de prédiction/correction : "+str(cpt)+"/1000",font = 'Times')
				window.update()
			x = ny.array(lx)
			y = ny.array(ly)
			
			if(val ==val_p[0]):
				plt.plot(x,y,color = lc[list_eps.index(eps)],label=str(eps))
			else : 
				plt.plot(x,y,color = lc[list_eps.index(eps)])

				
			i = i +1
			label_subtitle.config(text="nombre de valeur propre parcouru :" +str(i)+"/"+str(len(val_p)),font ='Times')
			c2 = can.create_rectangle(90,40,90+2*int(i/len(val_p)*100),60,outline="lime green",fill="lime green")
			can.delete(c1)
			c1 = c2
			can.update()
			window.update()
		i = 0
		k = 1
	window.destroy()

	title("Pseudospectre predcore")
	axis('equal')
	ax.legend()
	plt.show()







def gridtest(M,list_eps,pas,b1,b2,tol,taille_pred):
	"""
	fonction effectuant le calcul du pseudospectre via l'lagorithme GRID puis par prédiction correction et superpose les deux résultats
	"""

	#on définit un rectangle via les disques de Gerschgorin
	borne_X,borne_Y,liste_cercle,liste_c = limite(M,list_eps[len(list_eps)-1])
	i = complex(0,1)
	n = len(M)
	#on maille le plan
	x = ny.linspace(borne_X[0],borne_X[1],pas)
	y = ny.linspace(borne_Y[0],borne_Y[1],pas)
	sigmin = zeros((pas,pas))
	list_sig = []
	#gestion de la barre de chargement
	proc = 0
	window = Tk()
	window.title("processing")
	window.minsize(480,160)
	window.maxsize(480,160)
	#window.configure(bg = 'green')
	label_title1 = Label(window,text= "Calcul en cours :",font ='Times')
	label_title1.pack(expand = YES)
	label_subtitle = Label(window,text= "0%",font ='Times')
	label_subtitle.pack(expand = YES)
	can = Canvas(window)
	can.pack(expand=YES)
	can.create_rectangle(89,39,291,61,outline="black")
	c1 = can.create_rectangle(90,40,90,60,outline="lime green",fill="lime green")

	for l in range(pas):
		for j in range(pas):
			list_sig = svdvals((x[l]+y[j]*i)*eye(n)-M)
			sigmin[j][l] = list_sig[len(list_sig)-1]
			proc=proc+1
			label_subtitle.config(text=str(round(proc/(pas*pas)*100,1))+"%",font ='Times')
			c2 = can.create_rectangle(90,40,90+2*int(proc/(pas*pas)*100),60,outline="lime green",fill="lime green")
			can.delete(c1)
			c1 = c2
			can.update()
			window.update()
	window.destroy()

	if(len(list_eps)>1):
		normi = cm.colors.BoundaryNorm(boundaries = list_eps,ncolors = 256)
		contour(x,y,sigmin,list_eps,cmap = 'tab10',norm = normi) #Set1 tab10 tab20
		plt.colorbar()
	else : 
		contour(x,y,sigmin,list_eps,cmap = 'cool') #Set1 tab10 tab20		plt.colorbar()
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.gca()
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')

	if(b1 or b2):

		for i in range (len(liste_cercle)):
			if(b1):
				co,r1 = liste_cercle[i]
				ax.add_artist(plt.Circle(co,radius =r1,color ='k',fill=False,label="d2"))
			if(b2):
				co,r1 = liste_c[i]
				ax.add_artist(plt.Circle(co,ls = "--",radius =r1,color ='b',fill=False,label="d1"))
	#2ème partie
	lc = ['blue','cyan','m','g','r','k']
	n = len(M) 
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.gca()
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')
	window = Tk()
	window.title("processing")
	window.minsize(480,160)
	window.maxsize(480,160)
	label_title1 = Label(window,text= "Calcul en cours :",font ='Times')
	label_title1.pack(expand = YES)
	label_subtitle = Label(window,text= "nombre de valeur propre parcouru : 0/"+str(len(val_p)),font ='Times')
	label_subtitle.pack(expand = YES)
	label_subtitle2 = Label(window,text = "nombre de prédiction/correction : 0/1000",font = 'Times')
	label_subtitle2.pack(expand = YES)

	can = Canvas(window)
	can.pack(expand=YES)
	can.create_rectangle(89,39,291,61,outline="black")
	c1 = can.create_rectangle(90,40,90,60,outline="lime green",fill="lime green")
	i = 0
	cpt = 0
	k= 0

	for eps in list_eps :

		for val in val_p :
			if(eps == list_eps[0]):
				plt.plot(val.real,val.imag,'k:.',)
			z = cherche_z(M,val,eps,tol,n)
			b = True
			z_temp = step_2_3(z,M,taille_pred,n,eps,val)
			lx = []
			ly = []
			cpt = 0
			while b and cpt < 1000:
				z_temp = step_2_3(z_temp,M,taille_pred,n,eps,val)
				lx.append(z_temp.real)
				ly.append(z_temp.imag)
				b = norm(z-z_temp)>(taille_pred*0.3)
				cpt = cpt+1
				label_subtitle2.config(text ="nombre de prédiction/correction : "+str(cpt)+"/1000",font = 'Times')
				window.update()
			x = ny.array(lx)
			y = ny.array(ly)
			
			if(val ==val_p[0]):
				plt.plot(x,y,color = lc[list_eps.index(eps)],label=str(eps))
			else : 
				plt.plot(x,y,color = lc[list_eps.index(eps)])

				
			i = i +1
			label_subtitle.config(text="nombre de valeur propre parcouru :" +str(i)+"/"+str(len(val_p)),font ='Times')
			c2 = can.create_rectangle(90,40,90+2*int(i/len(val_p)*100),60,outline="lime green",fill="lime green")
			can.delete(c1)
			c1 = c2
			can.update()
			window.update()
		i = 0
		k = 1
	window.destroy()

	'''
	for eps in list_eps : 
		if(len(val_p ==64)):
			z = cherche_z(M,val_p[40],eps,tol,n)
			plt.plot(z.real,z.imag,'r:.')
	'''
	for i in range(len(val_p)):
		plt.plot(val_p[i].real,val_p[i].imag,'k:.',)
	title("Pseudospectre predcore")
	axis('equal')
	plt.show()

def gridtestbis(M,list_eps,pas,b1,b2,x_min,x_max,y_min,y_max,tol):
	"""
	fonction de zoom pour le la fonction gridtestbis
	"""
	x = ny.linspace(x_min,x_max,pas)
	y = ny.linspace(y_min,y_max,pas)
	i = complex(0,1)
	n = len(M)
	sigmin = zeros((pas,pas))
	list_sig = []
	proc = 0
	window = Tk()
	window.title("processing")
	window.minsize(480,160)
	window.maxsize(480,160)
	label_title1 = Label(window,text= "Calcul en cours :")
	label_title1.pack(expand = YES)
	label_subtitle = Label(window,text= "0%",font = 'Times')
	label_subtitle.pack(expand = YES)
	can = Canvas(window)
	can.pack(expand=YES)
	can.create_rectangle(89,39,291,61,outline="black")
	c1 = can.create_rectangle(90,40,90,60,outline="lime green",fill="lime green")
	for l in range(pas):
		for j in range(pas):
			list_sig = svdvals((x[l]+y[j]*i)*eye(n)-M)
			sigmin[j][l] = list_sig[len(list_sig)-1]
			proc=proc+1
			label_subtitle.config(text=str(round(proc/(pas*pas)*100,1))+"%",font ='Times')
			c2 = can.create_rectangle(90,40,90+2*int(proc/(pas*pas)*100),60,outline="lime green",fill="lime green")
			can.delete(c1)
			c1 = c2
			can.update()
			window.update()
	window.destroy()
			
	if(len(list_eps)>1):
		norm = 	cm.colors.BoundaryNorm(boundaries = list_eps,ncolors = 256)
		contour(x,y,sigmin,list_eps,cmap = 'cool',norm = norm) #Set1 tab10 tab20
		plt.colorbar()
	else : 
		contour(x,y,sigmin,list_eps,cmap = 'cool') #Set1 tab10 tab20		plt.colorbar()
	val_p,vecteurs=linalg.eig(M)
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.gca()
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')

	
	if(b1 or b2):
		for i in range (len(liste_cercle)):
			if(b1):
				co,r1 = liste_cercle[i]
				ax.add_artist(plt.Circle(co,radius =r1,color ='k',fill=False,label="d2"))
			if(b2):
				co,r1 = liste_c[i]
				ax.add_artist(plt.Circle(co,ls = "--",radius =r1,color ='b',fill=False,label="d1"))
	for eps in list_eps:
		z1 = cherche_z(M,val_p[0],eps,tol,n)
		if(z1.real>=x_min and z1.real<=x_max and z1.imag<= y_max and z1.imag>=y_min):
				plt.plot(z1.real,z1.imag,'g:.',)	

	for i in range(len(val_p)):
		if(val_p[i].real>=x_min and val_p[i].real<=x_max and val_p[i].imag<= y_max and val_p[i].imag>=y_min):
			if(i ==0):
				plt.plot(val_p[i].real,val_p[i].imag,'r:.',)
			else:
				plt.plot(val_p[i].real,val_p[i].imag,'k:.',)
	title("Pseudospectre GRID")
	axis('equal')

	plt.show()
	return list(ax.axis())


