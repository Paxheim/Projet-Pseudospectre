from grid import*
from compwise import*
from predcore import*
from abscisse import *
# différents tableaux de choix 
#image = ["fig1.png","fig3.png","fig5.png","fig6.png","fig4.png","fig2.png"],
f_mat = ["matrice aléatoire","matrice de Fibonacci","matrice de Toeplitz","matrice de Hilbert","matrice smoke","autre matrice de Toeplitz","spirale","matrice de Toeplitz aléatoire"]
choices = ["lancer l'Algorithme GRID","lancer l'Algorihtme de prédiction-correction", "Exit"]
fields=["Nouvelle matrice","changer d'algorihtme","Exit"]
k = ["x min","x max", "y min","y max","nombre de points"]
res =''
b1 = False
l = [] # liste contenant les coordonnées  de la fenêtre matplotlib juste avant sa destruction
reply=choicebox(msg=" Calcul de Pseudospectre",title="Projet Pima",choices=choices)
temp = 1
compwise_bool = False
while(temp==1):
	if(reply==choices[2]):
		temp = -1

	# partie pour l'algorithme GRID
	if (reply==choices[0]):
		if(res != fields[1] and res != fields[2]):
			n_mat = choicebox(msg = "choix de la matrice",title = "Projet Pima",choices = f_mat)
			compwise_bool = False

			if(n_mat == f_mat[0]):
			#cas particulier de la matrice aléatoire
				f=["Taille de la matrice","borne","Epsilon","nombres de points","tester l'algorihtme de prédiction-correction?(oui/non)","rendu en 3D?(oui/on)\nimpossible si l'on teste l'algo de prédiction correction","Calcul abscisse max et rayon max: (oui/non)"]
				n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice aléatoire",fields=f,values=(10,10,"1 2 3 4 5 6",10000,"non","non","non"))
				m = matrice_random(int(n[0]),int(n[1]))
				if(n[6]=='oui'):
					b1 = True
				else :
					b1 = False
				if(n[4]=="oui"):
					tol = enterbox(msg = "veuillez choisir une tolérance",title="Variables pour lancer le programme",default = "0.0001")
					taille_pred = enterbox(msg = "veuillez choisir une distance de prédiction ",title="Variables pour lancer le programme",default = "0.1")
					l = gridtest(m,[float(s) for s in list((n[2]).split())],int(sqrt(int(n[3]))),False,False,float(tol),float(taille_pred))
				elif(n[5]=='oui'):
					l = grid2(m,[float(s) for s in list((n[2]).split())],int(sqrt(int(n[3]))),False,False)
				else : 
					l = grid1(m,[float(s) for s in list((n[2]).split())],int(sqrt(int(n[3]))),b1,False)
			else : 
				f=["Taille de la matrice","Epsilon","nombres de points","tester l'algorithme de prédiction-correction?(oui/non)","rendu en 3D?(oui/on)\nimpossible si l'on teste l'algo de prédiction correction","Calcul abscisse max et rayon max: (oui/non)"]
				if(n_mat ==f_mat[1]):
					n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice de Fibonacci",fields=f,values=(15,"0.01 0.05 0.1 0.5 1",10000,"non","non","non"))
					m = matrice_fibo(int(n[0]))
				elif(n_mat ==f_mat[2]):
					n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice de Toeplitz",fields=f,values=(64,"0.0000001 0.000001 0.00001 0.0001 0.001 0.01 ",10000,"non","non","non"))
					m = matrice_top(int(n[0]))
				elif(n_mat ==f_mat[3]):
					n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice de Hilbert",fields=f,values=(20,"0.0000001 0.000001 0.00001 0.0001 0.001 0.01 1",10000,"non","non","non"))
					m = hilb(int(n[0]))
				elif(n_mat ==f_mat[4]):
					n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice smoke",fields=f,values=(64,"0.0000001 0.000001 0.00001 0.0001 0.001 0.01",10000,"non","non","non"))
					m = smoke(int(n[0]))
				elif(n_mat ==f_mat[5]):
					n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice smoke",fields=f,values=(64,"0.0000001 0.000001 0.00001 0.0001 0.001 0.01",10000,"non","non","non"))
					m = matrice_top2(int(n[0]))
				elif(n_mat ==f_mat[6]):
					n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice smoke",fields=f,values=(64,"0.0000001 0.000001 0.00001 0.0001 0.001 0.01",10000,"non","non","non"))
					m = top3(int(n[0]))
				elif(n_mat ==f_mat[7]):
					n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice smoke",fields=f,values=(32,"0.0000001 0.000001 0.00001 0.0001 0.001 0.01",10000,"non","oui","non"))
					m = matrice_top4(int(n[0]))
				if(n[5]=="oui"):
					b1 = True
				else : 
					b1 = False
				if(n[3]=="oui"):
					tol = enterbox(msg = "veuillez choisir une tolérance",title="Variables pour lancer le programme",default = "0.0001")
					taille_pred = enterbox(msg = "veuillez choisir une distance de prédiction ",title="Variables pour lancer le programme",default = "0.1")
					l = gridtest(m,[float(s) for s in list((n[1]).split())],int(sqrt(int(n[2]))),b1,False,float(tol),float(taille_pred))
				elif(n[4]=='oui'):
					l = grid2(m,[float(s) for s in list((n[1]).split())],int(sqrt(int(n[2]))),b1,False)
				else :
					l = grid1(m,[float(s) for s in list((n[1]).split())],int(sqrt(int(n[2]))),b1,False)	
			fields=["Nouvelle matrice","Zoomer","Pseudospectre par composante","changer d'algorihtme","Exit"]
			res=buttonbox(msg = "Et maintenant ?",title="GRID",choices=fields)

		if(res==fields[3]):
			l = []
			reply=choicebox(msg=" Calcul de Pseudospectre",title="Projet Pima",choices=choices)
		elif(res==fields[4]):
			temp=-1
		elif(res ==fields[2]):
			#pseudospectre par composante
			compwise_bool = True
			if(n_mat == f_mat[0]):
				l = compwise(m,[float(s) for s in list((n[2]).split())],int(sqrt(int(n[3]))),False,False)
			else : 
				l = compwise(m,[float(s) for s in list((n[1]).split())],int(sqrt(int(n[2]))),False,False)
			res=buttonbox(msg = "Et maintenant ?",title="GRID",choices=fields)
			if(res==fields[3]):
				l = []
				reply=choicebox(msg=" Calcul de Pseudospectre",title="Projet Pima",choices=choices)
			elif(res==fields[4]):
				temp=-1
		elif(res == fields[1]):
			#gestion du ZOOM
			if(compwise_bool):
				if(n_mat==f_mat[0]):
					s =multenterbox(msg= "veuillez précisez les nouvelles coordonnées souhaitées",title="Variables pour lancer le programme",fields=k,values=(l[0],l[1],l[2],l[3],int(n[3])))
					l = compwisebis(m,[float(s) for s in list((n[2]).split())],int(sqrt(int(s[4]))),False,False,float(s[0]),float(s[1]),float(s[2]),float(s[3]))
				else : 
					s =multenterbox(msg= "veuillez précisez les nouvelles coordonnées souhaitées",title="Variables pour lancer le programme",fields=k,values=(l[0],l[1],l[2],l[3],int(n[2])))
					l = compwisebis(m,[float(s) for s in list((n[1]).split())],int(sqrt(int(s[4]))),False,False,float(s[0]),float(s[1]),float(s[2]),float(s[3]))
				res=buttonbox(msg = "Et maintenant ?",title="GRID",choices=fields)
				if(res==fields[3]):
					l = []
					reply=choicebox(msg=" Calcul de Pseudospectre",title="Pseudospectre d'une matrice",choices=choices)
				elif(res==fields[4]):
					temp=-1
			else : 
				if(n_mat==f_mat[0]):
					s =multenterbox(msg= "veuillez précisez les nouvelles coordonnées souhaitées",title="Variables pour lancer le programme",fields=k,values=(l[0],l[1],l[2],l[3],int(n[3])))
					if(n[5]=='oui'):
						l = gridbis2(m,[float(s) for s in list((n[2]).split())],int(sqrt(int(s[4]))),b1,False,float(s[0]),float(s[1]),float(s[2]),float(s[3]))
					else :
						l = gridbis(m,[float(s) for s in list((n[2]).split())],int(sqrt(int(s[4]))),b1,False,float(s[0]),float(s[1]),float(s[2]),float(s[3]))
				else :
					s =multenterbox(msg= "veuillez précisez les nouvelles coordonnées souhaitées",title="Variables pour lancer le programme",fields=k,values=(l[0],l[1],l[2],l[3],int(n[2])))
					if(n[4]=='oui'):
						l = gridbis2(m,[float(s) for s in list((n[1]).split())],int(sqrt(int(s[4]))),False,False,float(s[0]),float(s[1]),float(s[2]),float(s[3]))
					else :
						l = gridbis(m,[float(s) for s in list((n[1]).split())],int(sqrt(int(s[4]))),b1,False,float(s[0]),float(s[1]),float(s[2]),float(s[3]))
				res=buttonbox(msg = "Et maintenant ?",title="GRID",choices=fields)
				if(res==fields[3]):
					l = []
					reply=buttonbox(msg=" Calcul de Pseudospectre",title="Pseudospectre d'une matrice",choices=choices)
				elif(res==fields[4]):
					temp=-1


	#algorithme de prédiction-correction
	if(reply == choices[1]):
		n_mat = choicebox(msg = "choix de la matrice",title = "Projet Pima",choices = f_mat)
		if(n_mat == f_mat[0]):
			f=["Taille de la matrice","borne","Epsilon","tolérance","pas de prédiction"]
			n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice aléatoire",fields=f,values=(10,10,"3",0.001,0.1))
			m = matrice_random(int(n[0]),int(n[1]))
			predcore(m,[float(s) for s in list((n[2]).split())],float(n[3]),float(n[4]))
		else : 
			f=["Taille de la matrice","Epsilon","tolérance","pas de prédiction"]
			if(n_mat ==f_mat[1]):
				n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice de Fibonacci",fields=f,values=(15,"0.05 0.1 0.5",0.001,0.1))
				m = matrice_fibo(int(n[0]))
			elif(n_mat ==f_mat[2]):
				n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice de Toeplitz",fields=f,values=(64,"0.01 ",0.001,0.1))
				m = matrice_top(int(n[0]))
			elif(n_mat ==f_mat[3]):
				n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice de Hilbert",fields=f,values=(20,"1",0.001,0.1))
				m = hilb(int(n[0]))
			elif(n_mat ==f_mat[4]):
				n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice smoke",fields=f,values=(32,"0.01",0.001,0.07))
				m = smoke(int(n[0]))
			elif(n_mat ==f_mat[5]):
				n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice smoke",fields=f,values=(32,"0.01",0.001,0.05))
				m = matrice_top2(int(n[0]))
			elif(n_mat ==f_mat[6]):
				n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice smoke",fields=f,values=(64,"0.0000001 0.000001 0.00001 0.0001 0.001 0.01",0.0001,0.1))
				m = top3(int(n[0]))
			elif(n_mat ==f_mat[7]):
				n=multenterbox(msg= "Veuillez rentrer les paramètres du programme\nInsérer les epsilon par ordre croissant en les séparant par un espace",title="matrice smoke",fields=f,values=(32,"0.0000001 0.000001 0.00001 0.0001 0.001 0.01",10000,0.0001,0.1))
				m = matrice_top4(int(n[0]))
			predcore(m,[float(s) for s in list((n[1]).split())],float(n[2]),float(n[3]))
		fields=["Nouvelle matrice","changer d'algorihtme","Exit"]
		res=buttonbox(msg = "Et maintenant ?",title="GRID",choices=fields)

		if(res==fields[1]):
				l = []
				reply=choicebox(msg=" Calcul de Pseudospectre",title="Projet Pima",choices=choices,)
		if(res==fields[2]):
			temp=-1
