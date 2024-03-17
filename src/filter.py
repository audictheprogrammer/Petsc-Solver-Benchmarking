import sys

if len(sys.argv) != 2:
    print("Utilisation: python3 script.py <chemin_vers_fichier>")
    sys.exit(1)

# Récupérer le chemin du fichier à lire à partir des arguments de la ligne de commande
chemin_fichier = sys.argv[1]

# Générer le nom du fichier de sortie
nom_fichier_sortie = chemin_fichier.replace('.txt', '_filtered.txt')

# Lignes à sélectionner
lignes_a_selectionner = ["MUMPS LU", "MUMPS Cho", "SuperLU", "Strumpack LU", "Petsc LU", "Petsc Cho",
                         "SuiteSparse LU", "SuiteSparse Cho", "SuiteSparse QR",
                        "FICHIER:", "nnz:","Norm of error","Time (sec):", 
                        "MatSolve", "MatLUFactor", 
                        "MatLUFactorSym", "MatLUFactorNum",
                        "MatCholFctrSym","MatCholFctrNum",
                        "MatQRFactorSym", "MatQRFactorNum",
                        "package used",
                        "total: nonzeros"]

# Liste pour stocker les lignes sélectionnées
lignes_selectionnees = []

# Lecture du fichier
with open(chemin_fichier, 'r') as f:
    # Parcourir chaque ligne du fichier
    for ligne in f:
        # Vérifier si la ligne doit être sélectionnée
        for motif in lignes_a_selectionner:
            if motif in ligne:
                # Ajouter la ligne sélectionnée à la liste
                lignes_selectionnees.append(ligne.strip())
                break  # Passer à la ligne suivante

# Écriture des lignes sélectionnées dans le fichier de sortie
with open(nom_fichier_sortie, 'w') as f:
    for ligne in lignes_selectionnees:
        f.write(ligne + '\n')
