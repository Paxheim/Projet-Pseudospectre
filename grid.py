from scipy.linalg import svdvals
from numpy import linalg
import numpy as ny
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
from pseudoAbscisseRayon import *



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
	fonction de zoom, fonctionnement similaire à grid1.
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
			p = abscissae(M,eps,axi = ax)
			#plt.plot(p.real,p.imag,"b:.")
			p = radii(M,eps,axi = ax)
			#plt.plot(p.real,p.imag,"r:.")
	
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





