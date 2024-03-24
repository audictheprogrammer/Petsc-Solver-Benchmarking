#!/bin/bash

# Définir la liste des fichiers

LIST_FICHIERS=("1138_bus.mtx")      # 1138 1138 2596
LIST_FICHIERS=("mhd4800b.mtx")      # 4800 4800 16160
LIST_FICHIERS=("bcsstk17.mtx")      # 10974 10974 219812
LIST_FICHIERS=("cvxbqp1.mtx")       # 50000 50000 199984
LIST_FICHIERS=("Dubcova2.mtx")      # 65025 65025 547625
LIST_FICHIERS=("thermomech_TC.mtx") # 102158 102158 406858
# LIST_FICHIERS=("offshore.mtx")      # 259789 259789 2251231 OUT OF MEMORY.



# Chemin vers le répertoire contenant les fichiers
wPETSC_DIR=/mnt/c/Users/audic/petsc/share/petsc/datafiles/matrices/MYMAT
# wPETSC_DIR=/mnt/c/Users/xu/petsc/share/petsc/datafiles/matrices/MYMAT

# Boucle à travers chaque fichier dans la liste
for fichier in "${LIST_FICHIERS[@]}"; do
    echo "FICHIER: $fichier"
    echo -e "#####################################"
    echo -e "############# MUMPS LU ##############"
    echo -e "#####################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type lu -pc_factor_mat_solver_type mumps -memory_view -log_view -ksp_view
    echo -e "\n"

    echo -e "#####################################"
    echo -e "############# MUMPS Cho #############"
    echo -e "#####################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type cholesky -pc_factor_mat_solver_type mumps -memory_view -log_view -ksp_view
    echo -e "\n"

    echo -e "#####################################"
    echo -e "############# SuperLU LU ############"
    echo -e "#####################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type lu -pc_factor_mat_solver_type superlu -memory_view -log_view -ksp_view
    echo -e "\n"

    echo -e "#####################################"
    echo -e "############ Strumpack LU ###########"
    echo -e "#####################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type lu -pc_factor_mat_solver_type strumpack -memory_view -log_view -ksp_view
    echo -e "\n\n"

    echo -e "#####################################"
    echo -e "########### SuiteSparse LU ##########"
    echo -e "#####################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type lu -pc_factor_mat_solver_type umfpack -memory_view -log_view -ksp_view
    echo -e "\n"

    echo -e "#####################################"
    echo -e "########### SuiteSparse Cho #########"
    echo -e "#####################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type cholesky -pc_factor_mat_solver_type cholmod -memory_view -log_view -ksp_view
    echo -e "\n"

    echo -e "#####################################"
    echo -e "########### SuiteSparse QR ##########"
    echo -e "#####################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type qr -pc_factor_mat_solver_type spqr -memory_view -log_view -ksp_view
    echo -e "\n"

    echo -e "#####################################"
    echo -e "############## Petsc LU #############"
    echo -e "#####################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type lu -pc_factor_mat_solver_type petsc -memory_view -log_view -ksp_view
    echo -e "\n"

    echo -e "######################################"
    echo -e "############## Petsc Cho #############"
    echo -e "######################################"
    ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -pc_type cholesky -pc_factor_mat_solver_type petsc -memory_view -log_view -ksp_view
    echo -e "\n"

done

#   ./ex72 -fin "${wPETSC_DIR}/${fichier}" -fout petscmat.aij -aij_only -ksp_monitor -pc_type lu -ksp_max_it 1 -pc_factor_mat_solver_type strumpack
