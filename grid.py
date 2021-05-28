from scipy.linalg import svdvals
from numpy import linalg
import numpy as ny
from random import*
from pylab import *
import matplotlib.pyplot as plt 
from matplotlib import cm 
import matplotlib.animation as animation 
from math import *
import sys
from tkinter import *
import subprocess as sp
from easygui import *
import os
from abscisse import *

seed()

class Higham(Exception):
    pass


def pdtoep(n, m=None, w=None, theta=None):
   
    if m is None:
        m = n

    if w is None:
        w = ny.random.rand(m)

    if theta is None:
        theta = ny.random.rand(m)

    try:
        if ny.max(w.shape) != m or ny.max(theta.shape) != m:
            raise Higham('Arguments w and theta must be vectors of length M.')
    except:
        raise Higham('Arguments w and theta must be vectors of length M.')

    t = ny.zeros(n)
    e = 2 * ny.pi * (ny.outer(ny.arange(1, n + 1), ny.ones(n)) -
                     ny.ones(n) * ny.arange(1, n + 1))

    for i in range(m):
        t = t + w[i] * ny.cos(theta[i] * e)

    return t



def chebspec(n, k=0):

    if k == 1:
        n += 1

    n = n - 1
    c = ny.zeros((n + 1, n + 1))

    one = ny.ones(n + 1)
    x = ny.cos(ny.arange(0, n + 1) * (ny.pi / float(n)))
    d = ny.ones(n + 1)
    d[0] = 2
    d[n] = 2

    # np.eye(n + 1) in next expression avoids div by zero.
    c = ny.outer(d, (one / d)) / (ny.outer(x, one) -
                                  ny.outer(one, x) + ny.eye(n + 1))

    #  Now fix diagonal and signs.
    c[0, 0] = (2 * n ** 2 + 1) / 6.0
    for i in range(1, n + 1):
        if ((i + 1) % 2) == 0:
            c[:, i] = -c[:, i]
            c[i, :] = -c[i, :]

        if i < n:
            c[i, i] = -x[i] / (2 * (1 - x[i] ** 2))
        else:
            c[n, n] = -c[0, 0]

    if k == 1:
        c = c[1::, 1::]

    return c



def matrice_random(n,borne):
	"""int*int->list(list(complex))"""

	res = [[0 for j in range(n)]for i in range(n)]
	for i in range(n) : 
		for j in range(n) : 
			res[i][j] = complex(randint(-borne,borne),randint(-borne,borne))
	return res

def matrice_top(n):
	a = complex(0,0)
	b = complex(0,0)
	c = complex(1/2,1/8)
	d = complex(1,0)
	e = complex(0,0)
	res = [[0 for j in range(n)]for i in range(n)]
	for i in range(n) : 
		for j in range(n) : 
			if(i==j ):
				res[i][j] = a
			elif(i == j-1):
				res[i][j] = b
			elif(i == j+1):
				res[i][j] = c
			elif(i == j-2):
				res[i][j] = d
			elif(i==j+2):
				res[i][j] = e
			else:
				res[i][j] = complex(0,0)
	return res


def matrice_top4(n):
	p = uniform(0,1)
	a =  complex(0,0)
	b =  complex(-p,0)
	c =  complex(-p,0)
	d =  complex(p,0)
	e =  complex(-p,0)
	f =  complex(p,0)
	g =  complex(p,0)
	res = [[0 for j in range(n)]for i in range(n)]
	for i in range(n) : 
		for j in range(n) : 
			if(i==j ):
				res[i][j] = a
			elif(i == j-1):
				res[i][j] = b
			elif(i == j+1):
				res[i][j] = c
			elif(i == j-2):
				res[i][j] = d
			elif(i==j+2):
				res[i][j] = e
			elif(i==j-3):
				res[i][j] = f
			elif(i==j+3):
				res[i][j] = g
			else:
				res[i][j] = complex(0,0)
	return res

def matrice_top5(n):
	a = complex(0,0)
	b = complex(0,0)
	c = complex(0,0)
	d = complex(0,0)
	e = complex(1,0)
	f = complex(1/2,0)
	g = complex(0,0)
	res = [[0 for j in range(n)]for i in range(n)]
	for i in range(n) : 
		for j in range(n) : 
			if(i==j ):
				res[i][j] = a
			elif(i == j-1):
				res[i][j] = b
			elif(i == j+1):
				res[i][j] = c
			elif(i == j-2):
				res[i][j] = d
			elif(i==j+2):
				res[i][j] = e
			elif(i==j-3):
				res[i][j] = f
			elif(i==j+3):
				res[i][j] = g
			else:
				res[i][j] = complex(0,0)
	return res

def matrice_top3(n):
	a = complex(0,0)
	b = complex(1,0)
	c = complex(0,0)
	d = complex(0,0)
	e = complex(1/2,0)
	res = [[0 for j in range(n)]for i in range(n)]
	for i in range(n) : 
		for j in range(n) : 
			if(i==j ):
				res[i][j] = a
			elif(i == j-1):
				res[i][j] = b
			elif(i == j+1):
				res[i][j] = c
			elif(i == j-2):
				res[i][j] = d
			elif(i==j+2):
				res[i][j] = e
			else:
				res[i][j] = complex(0,0)
	return res
def matrice_top2(n):
	a = complex(1,0)
	b = complex(-1,0)
	res = [[0 for j in range(n)]for i in range(n)]
	for i in range(n) : 
		for j in range(n) : 
			if(i==j or i==j-1 or i==j-2 or i==j-3 ):
				res[i][j] = a
			elif(i==j+1):
				res[i][j] = b
			else:
				res[i][j] = complex(0,0)
	return res


def smoke(n, k=0):

    w = ny.exp(2 * ny.pi * 1j / n)
    a = (ny.diag(np.hstack((w ** ny.arange(1, n), 1))) +ny.diag(ny.ones(n - 1), 1))
    if k == 0:
        a[n - 1, 0] = 1

    return a

def riemann(n):
	n = n + 1
	i = np.outer(np.arange(2, n + 1), np.ones(n - 1))
	j = i.T
	a = np.where(np.remainder(j, i) == 0, i, -1)
	return a

def hilb(n,m=0):
	if (n == 1 and (m == 0 or m == 1)):
		return np.array([[1]])
	elif m == 0:
		m = n

	v = np.arange(1, n + 1) + np.arange(0, m)[:, np.newaxis]
	return 1. / v

def matrice_fibo(n):
	res = [[0 for j in range(n)]for i in range(n)]
	for i in range(n) : 
		for j in range(n) : 
			if(i>=j):
				res[i][j]=complex(0,0)
			else:
				res[i][j] = complex(0,0)
				while(res[i][j].real == 0):
					res[i][j] = complex(randint(-1,2),0)
	while(res[n-1][0].real == 0):
		res[n-1][0] = complex(randint(-1,2),0)
	return res

def affiche_matrice(M,n):
	"""list(list(complex))*int->void"""
	for i in range(n):
		print()
		for j in range(n):
			if(M[i][j].real == 0):
				print("   ",end = '')
				print(M[i][j] , end = '')
				print("   ",end = '')
			elif(M[i][j].real == -1):
				print(M[i][j] , end = '')
				print(" ",end = '')
			else :
				print(M[i][j] , end = '')
				print("  ",end = '')
	print()


#détermine la portion du plan complexe sur laquelle appliquer l'algorithme GRID
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



def grid1(M,list_eps,pas,b1,b2):
	"""

	première implémentation de l'algorithme Grid, calcul le pseudospectre de la matrice M
	pour les valeurs de epsilon contenues dans list_eps(par ordre croissant en faisant pas*pas SVD. 
	affiche le pseudo absicce si b1 = True et le pseudo rayon si b2= True
	
	renvoie les coordonnées de la fenêtre matplotlib au moment de sa destruction (utile pour le zoom)

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

	#norm = 	cm.colors.Normalize(vmax = list_eps[len(list_eps)-1],vmin = list_eps[0])
 	
 	#affichage du pseudospectre
	if(len(list_eps)>1):
		norm = 	cm.colors.BoundaryNorm(boundaries = list_eps,ncolors = 256)
		contour(x,y,sigmin,levels = list_eps,cmap = 'cool',norm = norm) #Set1 tab10 tab20
		plt.colorbar()
	else : 
		contour(x,y,sigmin,list_eps,cmap = 'cool') #Set1 tab10 tab20		plt.colorbar()
	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.gca()
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')

	#affichage du pseudo abcisse
	if(b1):
		for eps in list_eps:
			p = abscissae(M,eps,axi = ax)
			#plt.plot(p.real,p.imag,"b:.")
			p = radii(M,eps,axi = ax)
			#plt.plot(p.real,p.imag,"r:.")
	#affichage du pseudo rayon
	if(b2):

		for i in range (len(liste_cercle)):
			if(b1):
				co,r1 = liste_cercle[i]
				ax.add_artist(plt.Circle(co,radius =r1,color ='k',fill=False,label="d2"))
			if(b2):
				co,r1 = liste_c[i]
				ax.add_artist(plt.Circle(co,ls = "--",radius =r1,color ='b',fill=False,label="d1"))
				
	for i in range(len(val_p)):
		plt.plot(val_p[i].real,val_p[i].imag,'k:.',)
	title("Pseudospectre GRID")
	axis('equal')
	plt.show()
	
	return list(ax.axis())

	


def gridbis(M,list_eps,pas,b1,b2,x_min,x_max,y_min,y_max):
	"""
	fonction de zoom similaire à grid1(), fonctionnement similaire à grid1.
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
		contour(x,y,sigmin,list_eps,cmap = 'cool') 

	val_p,vecteurs=linalg.eig(M)
	fig = plt.gcf()
	ax = fig.gca()
	fig.patch.set_facecolor('#c0c0c0')
	ax.patch.set_facecolor('#9d9d9d')

	if(b1):
		for eps in list_eps:
			p= abscissae(M,eps)
			#plt.plot(p.real,p.imag,"b:.")
	
	if(b2):
		for i in range (len(liste_cercle)):
			if(b1):
				co,r1 = liste_cercle[i]
				ax.add_artist(plt.Circle(co,radius =r1,color ='k',fill=False,label="d2"))
			if(b2):
				co,r1 = liste_c[i]
				ax.add_artist(plt.Circle(co,ls = "--",radius =r1,color ='b',fill=False,label="d1"))
				
	for i in range(len(val_p)):
		if(val_p[i].real>=x_min and val_p[i].real<=x_max and val_p[i].imag<= y_max and val_p[i].imag>=y_min):
			plt.plot(val_p[i].real,val_p[i].imag,'k:.',)
	title("Pseudospectre GRID")
	axis('equal')

	plt.show()
	#plt.savefig("fig0",format = "png")
	return list(ax.axis())





