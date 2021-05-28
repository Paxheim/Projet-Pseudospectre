from grid import*
from math import*
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import toeplitz
seed()

def top3(n):
	r = []
	for a in range(1,n+1):
		r.append(-complex(0,1)/a)
	c = [-complex(0,1)]
	c.append(pi)
	for a in range(1,n-1):
		c.append(complex(0,1)/a)

	return toeplitz(c,r)
def gearm(n, i=None, j=None):
#nice
    if i is None:
        i = n - 1
        j = -(n - 1)

    if not(abs(i) < n and abs(j) < n):
        raise ValueError('Invalid i and j parameters')

    a = np.diag(np.ones(n - 1), -1) + np.diag(np.ones(n - 1), 1)
    a[0, np.abs(i)] = np.sign(i)
    a[n - 1, n - 1 - np.abs(j)] = np.sign(j)

    return a

def frank(n, k=0):
    #nice
    f = np.ones((n, n)) * np.arange(n, 0, -1)
    f = np.minimum(f, f.T)
    #   take upper Hessenberg part.
    f = np.triu(f, -1)
    if k == 1:
        f = f[::-1, ::-1].T
    return f
class Higham(Exception):
    pass


def invhess(x, y=None):
#nice
    try:
        n = np.max(x.shape)
    except AttributeError:
        n = x
        x = np.arange(1, n + 1)
        y = -1 * x

    a = np.outer(np.ones(n), x)
    for j in range(1, n):
        a[:j - 1, j] = y[:j - 1]

    return a


def ohess(x):
#nice
    if type(x) == int:
        n = x
        x = np.random.uniform(size=n - 1) * 2 * np.pi
        h = np.eye(n)
        h[n - 1, n - 1] = np.sign(np.random.randn())
    elif type(x) == numpy.ndarray:
        if np.imag(x).any():
            raise Higham('Parameter must be real.')
        n = np.max(x.shape)
        h = np.eye(n)
        # Second term ensures h[n-1, n-1] nonzero.
        h[n - 1, n - 1] = np.sign(x[n - 1]) + float(x[n - 1] == 0)
    else:
        raise Higham('Unknown type in ohess')        
        
        
    for i in range(n - 1, 0, -1):
        # Apply Givens rotation through angle x[i - 1]
        theta = x[i - 1]
        c = np.cos(theta)
        s = np.sin(theta)
        h[i - 1:i + 1, :] = np.vstack((c * h[i - 1, :] + s * h[i, :],-s * h[i - 1, :] + c * h[i, :]))
    return h
def redheff(n):
#ok
    i = np.outer(np.arange(1, n + 1), np.ones(n))
    a = np.where(np.remainder(i.T, i) == 0, 1, 0)
    a[:, 0] = np.ones(n)

    return a

def triw(n, alpha=-1, k=0):
#bof
    try:
        m, n = n
    except TypeError:
        m = n

    if k == 0:
        k = n - 1

    if np.array(alpha).size != 1:
        raise Higham("Second Argument Must Be A Scalar.")

    t = np.tril(np.eye(m, n) + alpha * np.triu(np.ones((m, n)), 1), k)

    return t



def compwise(M,list_eps,pas,b1,b2):
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
	E = np.copy(M)
	for i1 in range(n):
		for j1 in range(n):
			E[i1][j1] = abs(E[i1][j1])

	for l in range(pas):
		for j in range(pas):
			A = M-(x[l]+y[j]*i)*eye(n)
			A = linalg.inv(A)
			for i1 in range(n):
				for j1 in range(n):
					A[i1][j1] = abs(A[i1][j1])
			Y = np.dot(A,E)
			val_p,vecteurs=linalg.eig(Y)
			sigmin[j][l] = log(max([abs(v) for v in val_p]))# max([abs(v) for v in val_p])
			proc=proc+1

			label_subtitle.config(text=str(round(proc/(pas*pas)*100,1))+"%",font ='Times')
			c2 = can.create_rectangle(90,40,90+2*int(proc/(pas*pas)*100),60,outline="lime green",fill="lime green")
			can.delete(c1)
			c1 = c2
			can.update()
			window.update()
	window.destroy()

	#norm = 	cm.colors.Normalize(vmax = list_eps[len(list_eps)-1],vmin = list_eps[0])
	sval_p,vecteurs=linalg.eig(M)
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.add_subplot(111, projection='3d')
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')


	X, Y = np.meshgrid(x, y)
	ax.plot_surface(X,Y,sigmin,cmap = 'magma_r',alpha = 0.9,)#twilight_shifted terrain
	#ax.contour(X,Y,sigmin,levels = [log(1/e) for e in list_eps].reverse(),colors = 'black',alpha = 1)
	#cset = ax.contourf(X, Y, sigmin, zdir='z',cmap = 'summer',norm = norm)
	#ax.invert_zaxis()

	if(b1 or b2):

		for i in range (len(liste_cercle)):
			if(b1):
				co,r1 = liste_cercle[i]
				ax.add_artist(plt.Circle(co,radius =r1,color ='k',fill=False,label="d2"))
			if(b2):
				co,r1 = liste_c[i]
				ax.add_artist(plt.Circle(co,ls = "--",radius =r1,color ='b',fill=False,label="d1"))
	'''	
	for i in range(len(val_p)):
		ax.scatter3D(val_p[i].real,val_p[i].imag,0,c = 'k',alpha = 1,depthshade = False)
	'''
	title("Pseudospectre GRID")
	#axis('equal')
	plt.show()
	return list(ax.axis())




def grid2(M,list_eps,pas,b1,b2):
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
			#sigmin[j][l] = log10(1/(list_sig[len(list_sig)-1]))
			sigmin[j][l] = log(list_sig[len(list_sig)-1])
			proc=proc+1
			label_subtitle.config(text=str(round(proc/(pas*pas)*100,1))+"%",font ='Times')
			c2 = can.create_rectangle(90,40,90+2*int(proc/(pas*pas)*100),60,outline="lime green",fill="lime green")
			can.delete(c1)
			c1 = c2
			can.update()
			window.update()
	window.destroy()
	#norm = 	cm.colors.Normalize(vmax = list_eps[len(list_eps)-1],vmin = list_eps[0])
	sval_p,vecteurs=linalg.eig(M)
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.add_subplot(111, projection='3d')
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')


	X, Y = np.meshgrid(x, y)
	ax.plot_surface(X,Y,sigmin,cmap = 'magma',shade = True)#twilight_shifted terrain magma
	#ax.contour(X,Y,sigmin,levels = [log(e) for e in list_eps],colors = 'black',alpha = 1) #[log10(1/e) for e in list_eps].reverse()
	#cset = ax.contourf(X, Y, sigmin, zdir='z',cmap = 'summer',norm = norm)
	ax.invert_zaxis()

	if(b1 or b2):

		for i in range (len(liste_cercle)):
			if(b1):
				co,r1 = liste_cercle[i]
				ax.add_artist(plt.Circle(co,radius =r1,color ='k',fill=False,label="d2"))
			if(b2):
				co,r1 = liste_c[i]
				ax.add_artist(plt.Circle(co,ls = "--",radius =r1,color ='b',fill=False,label="d1"))
	'''	
	for i in range(len(val_p)):
		ax.scatter3D(val_p[i].real,val_p[i].imag,0,c = 'k',alpha = 1,depthshade = False)
	'''
	title("Pseudospectre GRID")
	#axis('equal')
	plt.show()
	return list(ax.axis())



def gridbis2(M,list_eps,pas,b1,b2,x_min,x_max,y_min,y_max):
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
			sigmin[j][l] = log(list_sig[len(list_sig)-1])
			proc=proc+1
			label_subtitle.config(text=str(round(proc/(pas*pas)*100,1))+"%",font ='Times')
			c2 = can.create_rectangle(90,40,90+2*int(proc/(pas*pas)*100),60,outline="lime green",fill="lime green")
			can.delete(c1)
			c1 = c2
			can.update()
			window.update()
	window.destroy()
	

	sval_p,vecteurs=linalg.eig(M)
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.add_subplot(111, projection='3d')
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')


	X, Y = np.meshgrid(x, y)
	ax.plot_surface(X,Y,sigmin,cmap = 'magma',shade = True)
	#ax.contour(X,Y,sigmin,levels = [log(e) for e in list_eps],colors = 'black',alpha = 1)
	#cset = ax.contourf(X, Y, sigmin, zdir='z',cmap = 'summer',norm = norm)
	ax.invert_zaxis()


	if(b1 or b2):

		for i in range (len(liste_cercle)):
			if(b1):
				co,r1 = liste_cercle[i]
				ax.add_artist(plt.Circle(co,radius =r1,color ='k',fill=False,label="d2"))
			if(b2):
				co,r1 = liste_c[i]
				ax.add_artist(plt.Circle(co,ls = "--",radius =r1,color ='b',fill=False,label="d1"))
	'''		
	for i in range(len(val_p)):
		ax.scatter3D(val_p[i].real,val_p[i].imag,0,c = 'k',alpha = 1)
	'''
	title("Pseudospectre GRID")
	#axis('equal')
	plt.show()
	return list(ax.axis())



def compwisebis(M,list_eps,pas,b1,b2,x_min,x_max,y_min,y_max):
	#on définit un rectangle via les disques de Gerschgorin
	i = complex(0,1)
	n = len(M)
	#on maille le plan
	x = ny.linspace(x_min,x_max,pas)
	y = ny.linspace(y_min,y_max,pas)
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
	E = np.copy(M)
	for i1 in range(n):
		for j1 in range(n):
			E[i1][j1] = abs(E[i1][j1])

	for l in range(pas):
		for j in range(pas):
			A = M-(x[l]+y[j]*i)*eye(n)
			A = linalg.inv(A)
			for i1 in range(n):
				for j1 in range(n):
					A[i1][j1] = abs(A[i1][j1])
			Y = np.dot(A,E)
			val_p,vecteurs=linalg.eig(Y)
			sigmin[j][l] = log(max([abs(v) for v in val_p]))# max([abs(v) for v in val_p])
			proc=proc+1
			label_subtitle.config(text=str(round(proc/(pas*pas)*100,1))+"%",font ='Times')
			c2 = can.create_rectangle(90,40,90+2*int(proc/(pas*pas)*100),60,outline="lime green",fill="lime green")
			can.delete(c1)
			c1 = c2
			can.update()
			window.update()
	window.destroy()

	#norm = 	cm.colors.Normalize(vmax = list_eps[len(list_eps)-1],vmin = list_eps[0])
	sval_p,vecteurs=linalg.eig(M)
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.add_subplot(111, projection='3d')
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')


	X, Y = np.meshgrid(x, y)
	ax.plot_surface(X,Y,sigmin,cmap = 'magma_r',alpha = 0.9,)#twilight_shifted terrain
	#ax.contour(X,Y,sigmin,levels = [log(1/e) for e in list_eps].reverse(),colors = 'black',alpha = 1)
	#cset = ax.contourf(X, Y, sigmin, zdir='z',cmap = 'summer',norm = norm)
	#ax.invert_zaxis()

	if(b1 or b2):

		for i in range (len(liste_cercle)):
			if(b1):
				co,r1 = liste_cercle[i]
				ax.add_artist(plt.Circle(co,radius =r1,color ='k',fill=False,label="d2"))
			if(b2):
				co,r1 = liste_c[i]
				ax.add_artist(plt.Circle(co,ls = "--",radius =r1,color ='b',fill=False,label="d1"))
	'''	
	for i in range(len(val_p)):
		ax.scatter3D(val_p[i].real,val_p[i].imag,0,c = 'k',alpha = 1,depthshade = False)
	'''
	title("Pseudospectre GRID")
	#axis('equal')
	plt.show()
	return list(ax.axis())



