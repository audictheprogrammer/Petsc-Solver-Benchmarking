import matplotlib.pyplot as plt
import sys


""" Manuel d'utilisation.
Le paramètre OPTION_Y correspond au choix de l'ordonnée.
case 1: Temps d'exécution total des solveurs.
case 2: Temps de MatSolve total des solveurs.
case 3: Temps de MatFactor total des solveurs.
case 4: Norme d'erreur des solveurs.
case 5: Nombre de non zero facto.
"""
OPTION_Y = 5


# Données
## Abscisse
nombre_non_zero = [2596, 16160,	219812, 199984, 547625,	406858]
taille_matrice = [1138, 4800, 10974, 50000, 65025, 102158]

# Ordonnées 1: Temps
temps_mumps_lu = [2.605e-01, 5.039e-01,	2.232e+00, 2.779e+00,	6.520e+00,	5.382e+00]
temps_mumps_cho = [2.603e-01,	3.843e-01,	2.233e+00, 2.168e+00,	6.201e+00,	5.484e+00]
temps_super_lu = [2.291e-01, 	2.698e-01, 	2.226e+00, 3.866e+00,	7.656e+00,	5.519e+00]
temps_strumpack_lu = [2.616e-01,	3.317e-01,	2.191e+00, 2.098e+00,	6.696e+00,	5.312e+00]
temps_petsc_lu = [1.840e-01,	3.375e-01,	1.955e+00,  1.57845e-09, 6.337e+00,	4.735e+00]
temps_petsc_cho = [2.765e-01,	2.939e-01,	2.326e+00,	1.725e+01,  5.989e+00,	5.046e+00]
temps_suitesparse_lu = [2.580e-01,	2.921e-01, 	2.151e+00, 3.931e+00, 	5.939e+00,	4.915e+00]
temps_suitesparse_cho = [2.765e-01,	2.825e-01,	1.998e+00, 2.242e+00, 	5.796e+00,	4.463e+00]
temps_suitesparse_qr = [2.148e-01,	3.340e-01,	3.024e+00, 6.017e+00,	1.547e+01,	7.560e+00]

# Ordonnées 2: MatSolve PAS OK, il manque la colonne 4 ou 3(si on compte a partir de 0)
matsolve_mumps_lu = [5.4342e-03,	8.7046e-03,	1.7164e-02,	3.8196e-02,	5.6500e-02]
matsolve_mumps_cho = [3.2611e-03,	5.2407e-03,	1.9493e-02,	3.5768e-02,	4.8872e-02]
matsolve_super_lu = [1.3120e-04,	4.2990e-04,	2.1611e-02,	8.6814e-02,	7.4688e-02]
matsolve_strumpack_lu = [1.0737e-03,	2.4930e-03,	9.3431e-03,	6.0511e-02,	8.1426e-02]
matsolve_petsc_lu = [5.2900e-05,	3.1060e-04,	6.8581e-03,	1.7044e-02,	2.4162e-02]
matsolve_petsc_cho = [3.4720e-03,	4.0672e-03,	1.6291e-02,	2.8998e-02,	4.3356e-02]
matsolve_suitesparse_lu = [5.6668e-03,	4.3212e-03,	1.5555e-02,	3.0045e-02,	3.6883e-02]
matsolve_suitesparse_cho = [3.4720e-03,	2.5126e-03,	1.0904e-02,	2.7880e-02,	3.4997e-02]
matsolve_suitesparse_qr = [8.5420e-04,	1.8147e-03,	7.2991e-02,	3.9876e-01,	3.6159e-01]

# Ordonnées 3: MatFactor PAS OK
matfactor_mumps_lu = [0.072903, 0.042564, 0.33525, 0.676156, 0.78995]
matfactor_mumps_cho = [0.075533, 0.079379, 0.31879, 0.93573, 1.04705]
matfactor_super_lu = [1.6334e-03,	3.5307e-03,	3.7970e-01,	2.3460e+00,	1.4125e+00]
matfactor_strumpack_lu = [4.7358e-02,	4.6739e-02,	3.1306e-01,	1.0609e+00,	1.1229e+00]
matfactor_petsc_lu = [0.0005317, 0.001367, 0.334067, 0.6112299999999999, 0.9485]
matfactor_petsc_cho = [0.0004017, 0.0008109, 0.50456, 0.5290900000000001, 1.04817]
matfactor_suitesparse_lu = [0.0073056, 0.0104796, 0.179685, 0.575442, 0.69797]
matfactor_suitesparse_cho = [3.0170e-04,  4.4446e-03,  7.7877e-02,  2.8221e-01,  3.8227e-01]
matfactor_suitesparse_qr =  [9.7154e-03,  1.3436e-02,  1.3634e+00,  9.7290e+00,  3.0638e+00]

# Ordonnées 4: Norme PAS OK
norme_mumps_lu = [6.09739e-11,	1.06878e-11,	5.44022e-11,	1.09771e-11,	1.2898e-10]
norme_mumps_cho = [7.64738e-11,	1.1899e-11,	9.12008e-11,	1.06573e-11,	1.2898e-10]
norme_super_lu = [1.81161e-10,	7.16011e-12,	8.66171e-11,	8.10018e-12,	1.28981e-10]
norme_strumpack_lu = [6.19828e-11,	1.15536e-11,	2.39435e-10,	1.56831e-11,	1.2898e-10]
norme_petsc_lu = [9.55287e-11,	8.43788e-12,	2.38356e-10,	8.63165e-12,	1.28981e-10]
norme_petsc_cho = [7.69126e-11,	9.07803e-12,	8.66133e-10,	1.08716e-11,	1.28983e-10]
norme_suitesparse_lu = [8.99435e-11,	1.14642e-11,	8.5205e-11,	7.72324e-12,	1.2898e-10]
norme_suitesparse_cho = [7.49487e-11,	9.30665e-12,	1.06056e-10,	6.9209e-12,	1.28976e-10]
norme_suitesparse_qr = [1.394e-10,	4.64286,	7.40947e-09,	7.34749e-12,	1.31537e-10]

# Ordonnées 5: NNZFacto
nnzfacto_mumps_lu = [24180,	72776, 	2456418, 	5293510,	7394695,	8883348]
nnzfacto_mumps_cho = [12659, 38788,	1238104,	2711187,	3731157,	4500436]
nnzfacto_super_lu = [0, 0, 0, 0, 0, 0]
nnzfacto_strumpack_lu = [0, 0, 0, 0, 0, 0]
nnzfacto_petsc_lu = [11730,	38284,	2966252,	19185240,	5517827,	8933756]
nnzfacto_petsc_cho = [6434,	21542,	1488613,	9617620,	2791426,	4517957]
nnzfacto_suitesparse_lu = [0, 0, 0, 0, 0, 0]
nnzfacto_suitesparse_cho = [3265, 16160, 1043601,	1915433,	2956869,	3681339]
nnzfacto_suitesparse_qr = [9002, 28076,	3095965,	7821888,	20290483,	12417253]


solveurs = ['MUMPS LU', 'MUMPS Cho', 'SuperLU', 'Strumpack LU', 'Petsc LU', 'Petsc Cho',
            'SuiteSparse LU', 'SuiteSparse Cho', 'SuiteSparse QR']

# plt.figure(figsize=(10, 6))
# for i, solveur in enumerate(solveurs):
#     plt.plot(nombre_non_zero,nnzfacto_mumps_lu[i], nnzfacto_mumps_cho[i], nnzfacto_super_lu[i],
#                                    nnzfacto_strumpack_lu[i], nnzfacto_petsc_lu[i], nnzfacto_petsc_cho[i],
#                                    nnzfacto_suitesparse_lu[i], nnzfacto_suitesparse_cho[i],
#                                    nnzfacto_suitesparse_qr[i]],
#                 label=solveur)

# # Personnalisation de l'axe des x et des y, des étiquettes et du titre
# plt.xlabel('Nombre de non-zéros après la factorisation')
# plt.ylabel('Nombre de non-zéros avant la factorisation')
# plt.title('Nombre de non-zéros avant vs après la factorisation pour chaque solveur')

# # Ajouter la légende
# plt.legend()

# # Afficher le graphique
# plt.grid(True)
# plt.show()

if len(sys.argv) != 2:
    print("Utilisation: python3 script.py <chemin_vers_fichier>")
    sys.exit(1)

Final_X = []            # Abscisse
Final_Y = []            # Ordonnée
Final_L = []            # Label
Final_Y_label = str()   # YLabel
Final_O = sys.argv[1]   # Savefig


match OPTION_Y:
    case 1:
        Final_Y_label = 'Temps d\'exécution'
        Final_Y.append(temps_mumps_lu)
        Final_L.append("MUMPS LU")
        Final_Y.append(temps_mumps_cho)
        Final_L.append("MUMPS CHO")
        Final_Y.append(temps_super_lu)
        Final_L.append("SUPER LU")
        Final_Y.append(temps_strumpack_lu)
        Final_L.append("STRUMPACK LU")
        Final_Y.append(temps_petsc_lu)
        Final_L.append("PETSC LU")
        Final_Y.append(temps_petsc_cho)
        Final_L.append("PETSC CHO")
        Final_Y.append(temps_suitesparse_lu)
        Final_L.append("SUITESPARSE LU")
        Final_Y.append(temps_suitesparse_cho)
        Final_L.append("SUITESPARSE CHO")
        Final_Y.append(temps_suitesparse_qr)
        Final_L.append("SUITESPARSE QR")

    case 2:
        Final_Y_label = 'Temps d\'exécution de MatSolve'
        Final_Y.append(matsolve_mumps_lu)
        Final_L.append("MUMPS LU")
        Final_Y.append(matsolve_mumps_cho)
        Final_L.append("MUMPS CHO")
        Final_Y.append(matsolve_super_lu)
        Final_L.append("SUPER LU")
        Final_Y.append(matsolve_strumpack_lu)
        Final_L.append("STRUMPACK LU")
        Final_Y.append(matsolve_petsc_lu)
        Final_L.append("PETSC LU")
        Final_Y.append(matsolve_petsc_cho)
        Final_L.append("PETSC CHO")
        Final_Y.append(matsolve_suitesparse_lu)
        Final_L.append("SUITESPARSE LU")
        Final_Y.append(matsolve_suitesparse_cho)
        Final_L.append("SUITESPARSE CHO")
        Final_Y.append(matsolve_suitesparse_qr)
        Final_L.append("SUITESPARSE QR")


    case 3:
        Final_Y_label = 'Temps d\'exécution de MatFactor'
        Final_Y.append(matfactor_mumps_lu)
        Final_L.append("MUMPS LU")
        Final_Y.append(matfactor_mumps_cho)
        Final_L.append("MUMPS CHO")
        Final_Y.append(matfactor_super_lu)
        Final_L.append("SUPER LU")
        Final_Y.append(matfactor_strumpack_lu)
        Final_L.append("STRUMPACK LU")
        Final_Y.append(matfactor_petsc_lu)
        Final_L.append("PETSC LU")
        Final_Y.append(matfactor_petsc_cho)
        Final_L.append("PETSC CHO")
        Final_Y.append(matfactor_suitesparse_lu)
        Final_L.append("SUITESPARSE LU")
        Final_Y.append(matfactor_suitesparse_cho)
        Final_L.append("SUITESPARSE CHO")
        Final_Y.append(matfactor_suitesparse_qr)
        Final_L.append("SUITESPARSE QR")
    case 4:
        Final_Y_label = 'Norme d\' erreur'
        Final_Y.append(norme_mumps_lu)
        Final_L.append("MUMPS LU")
        Final_Y.append(norme_mumps_cho)
        Final_L.append("MUMPS CHO")
        Final_Y.append(norme_super_lu)
        Final_L.append("SUPER LU")
        Final_Y.append(norme_strumpack_lu)
        Final_L.append("STRUMPACK LU")
        Final_Y.append(norme_petsc_lu)
        Final_L.append("PETSC LU")
        Final_Y.append(norme_petsc_cho)
        Final_L.append("PETSC CHO")
        Final_Y.append(norme_suitesparse_lu)
        Final_L.append("SUITESPARSE LU")
        Final_Y.append(norme_suitesparse_cho)
        Final_L.append("SUITESPARSE CHO")
        Final_Y.append(norme_suitesparse_qr)
        Final_L.append("SUITESPARSE QR")
    case 5:
        Final_Y_label = 'Nombre de non-zero Facto'
        Final_Y.append(nnzfacto_mumps_lu)
        Final_L.append("MUMPS LU")
        Final_Y.append(nnzfacto_mumps_cho)
        Final_L.append("MUMPS CHO")
        Final_Y.append(nnzfacto_super_lu)
        Final_L.append("SUPER LU = 0")
        Final_Y.append(nnzfacto_strumpack_lu)
        Final_L.append("STRUMPACK LU = 0")
        Final_Y.append(nnzfacto_petsc_lu)
        Final_L.append("PETSC LU")
        Final_Y.append(nnzfacto_petsc_cho)
        Final_L.append("PETSC CHO")
        Final_Y.append(nnzfacto_suitesparse_lu)
        Final_L.append("SUITESPARSE LU = 0")
        Final_Y.append(nnzfacto_suitesparse_cho)
        Final_L.append("SUITESPARSE CHO")
        Final_Y.append(nnzfacto_suitesparse_qr)
        Final_L.append("SUITESPARSE QR")


# Tracer les données
plt.figure(figsize=(10, 6))

plt.plot(taille_matrice, Final_Y[0], marker='o', label=Final_L[0])
plt.plot(taille_matrice, Final_Y[1], marker='o', label=Final_L[1])
plt.plot(taille_matrice, Final_Y[2], marker='o', label=Final_L[2])
plt.plot(taille_matrice, Final_Y[3], marker='o', label=Final_L[3])
plt.plot(taille_matrice, Final_Y[4], marker='o', label=Final_L[4])
plt.plot(taille_matrice, Final_Y[5], marker='o', label=Final_L[5])
plt.plot(taille_matrice, Final_Y[6], marker='o', label=Final_L[6])
plt.plot(taille_matrice, Final_Y[7], marker='o', label=Final_L[7])
plt.plot(taille_matrice, Final_Y[8], marker='o', label=Final_L[8])

plt.xlabel('Nombre de nnz initial')
plt.ylabel(Final_Y_label)
plt.title(Final_Y_label + ' en fonction du nnz itinial')
plt.legend()
plt.grid(True)
plt.savefig(Final_O)
print("Plot avec succès dans le fichier:", Final_O)
plt.show()

