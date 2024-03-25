static char help[] = "Distributed: Loads a MTX matrix then solves a linear system with KSP.\n\n";

#include <petsc.h>
#include "mmloader.h"

/*
export wPETSC_DIR=/mnt/c/Users/audic/petsc/share/petsc/datafiles/matrices/MYMAT
export wPETSC_DIR=/mnt/c/Users/xu/petsc/share/petsc/datafiles/matrices/MYMAT

mpiexec -n 2 ./ex23 -fin ${wPETSC_DIR}/cvxbqp1.mtx petscmat.aij -aij_only -pc_type lu -pc_factor_mat_solver_type mumps -memory_view -log_view -ksp_view> TEST2.txt
*/

int main(int argc, char **args)
{
  MM_typecode matcode;
  FILE       *file;

  Vec         x, b, u;                                   /* approx solution, RHS, exact solution */
  Mat         A, A_local, A_MPI;                         /* linear system matrix */
  KSP         ksp;                                       /* linear solver context */
  PC          pc;                                        /* preconditioner context */
  PetscReal   norm;                                     /* norm of solution error */
  PetscInt    its, zero = 0, one = 1;

  PetscInt    n = PETSC_DETERMINE, M, N, nz;
  char        filein[PETSC_MAX_PATH_LEN];
  char        ordering[256] = MATORDERINGRCM;
  // PetscViewer view;
  PetscBool   flag, aijonly = PETSC_FALSE, permute = PETSC_FALSE;
  IS          rowperm = NULL, colperm = NULL;
  IS          row_is, col_is;
  PetscMPIInt size, rank;


  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, (char *)0, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL));


  /* Loading matrix. */
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Matrix Market example options", "");
  {
    PetscCall(PetscOptionsString("-fin", "Input Matrix Market file", "", filein, filein, sizeof(filein), &flag));
    PetscCheck(flag, PETSC_COMM_SELF, PETSC_ERR_USER_INPUT, "Please use -fin <filename> to specify the input file name!");
    PetscCall(PetscOptionsBool("-aij_only", "Use MATAIJ for all cases", "", aijonly, &aijonly, NULL));
    PetscCall(PetscOptionsFList("-permute", "Permute matrix and vector to solving in new ordering", "", MatOrderingList, ordering, ordering, sizeof(ordering), &permute));
  }
  PetscOptionsEnd();

  PetscCall(MatCreateFromMTX(&A, filein, aijonly));
  PetscCall(PetscFOpen(PETSC_COMM_SELF, filein, "r", &file));
  PetscCallExternal(mm_read_banner, file, &matcode);
  PetscCallExternal(mm_write_banner, stdout, matcode);
  PetscCallExternal(mm_read_mtx_crd_size, file, &M, &N, &nz);
  PetscCall(PetscFClose(PETSC_COMM_SELF, file));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "M: %d, N: %d, nnz: %d\n", M, N, nz));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "Reading matrix completes.\n"));
  if (permute) {
    Mat Aperm;
    PetscCall(MatGetOrdering(A, ordering, &rowperm, &colperm));
    PetscCall(MatPermute(A, rowperm, colperm, &Aperm));
    PetscCall(MatDestroy(&A));
    A = Aperm; /* Replace original operator with permuted version */
  }

  /* Matrix loaded. */


  /* Creating local matrices, vectors avec local vectors. */
  PetscCall(PetscSplitOwnership(PETSC_COMM_WORLD, &n, &N));
  /* Fully working, but ugly. */
  int starting_pos = 0;
  for (int i = 0; i < rank; i++){
    int n2 = N / size;

    int k = N % n2;
    if (n2 == 1){
      k = N - size;
    }

    if (k != 0 && i < k){
      n2++;
    }
    starting_pos += n2;
  }


  PetscCall(ISCreateStride(PETSC_COMM_WORLD, M, zero, one, &col_is));             // Full columns.
  PetscCall(ISCreateStride(PETSC_COMM_WORLD, n, starting_pos, one, &row_is));     // Partiel rows.

  PetscCall(MatCreateSubMatrix(A, row_is, col_is, PETSC_DECIDE, &A_local));
  PetscCall(MatCreateMPIMatConcatenateSeqMat(PETSC_COMM_WORLD, A_local, n, MAT_INITIAL_MATRIX, &A_MPI));

  // PetscCall(PetscPrintf(PETSC_COMM_SELF, "RANK = %d, n = %d.\n", rank, n));
  // PetscCall(MatView(A_MPI, PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(VecCreate(PETSC_COMM_WORLD, &x));
  PetscCall(PetscObjectSetName((PetscObject)x, "Solution"));
  PetscCall(VecSetSizes(x, n, N));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecDuplicate(x, &b));
  PetscCall(VecDuplicate(x, &u));

  // PetscCall(PetscPrintf(PETSC_COMM_SELF, "RANK = %d, n = %d.\n", rank, n));
  // PetscCall(VecView(x, PETSC_VIEWER_STDOUT_WORLD));


  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  Compute the matrix and right-hand-side vector that define
  the linear system, Ax = b.
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Assemble matrix.

     The linear system is distributed across the processors by
     chunks of contiguous rows, which correspond to contiguous
     sections of the mesh on which the problem is discretized.
     For matrix assembly, each processor contributes entries for
     the part that it owns locally.
  */


  /* Assemble the matrix */
  // USEFULL ??
  // PetscCall(MatAssemblyBegin(A_MPI, MAT_FINAL_ASSEMBLY));
  // PetscCall(MatAssemblyEnd(A_MPI, MAT_FINAL_ASSEMBLY));

  /*
     Set exact solution; then compute right-hand-side vector.
  */

  PetscCall(VecSet(u, 1.0));
  PetscCall(MatMult(A_MPI, u, b));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create linear solver context
  */
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
  */
  PetscCall(KSPSetOperators(ksp, A_MPI, A_MPI));

  /*
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following four statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions();
  */
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCJACOBI));
  PetscCall(KSPSetTolerances(ksp, 1.e-5, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));

  /*
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  PetscCall(KSPSetFromOptions(ksp));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Solve linear system
  */
  PetscCall(KSPSolve(ksp, b, x));

  /*
     View solver info; we could instead use the option -ksp_view to
     print this info to the screen at the conclusion of KSPSolve().
  */
  // PetscCall(KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Check the error
  */
  PetscCall(VecAXPY(x, -1.0, u));
  PetscCall(VecNorm(x, NORM_2, &norm));
  PetscCall(KSPGetIterationNumber(ksp, &its));
  // if (norm > tol) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Norm of error %g, Iterations %" PetscInt_FMT "\n", (double)norm, its));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "Norm of error %g, Iterations %" PetscInt_FMT "\n", (double)norm, its));

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(MatDestroy(&A_local));
  PetscCall(MatDestroy(&A_MPI));
  PetscCall(KSPDestroy(&ksp));

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_view).
  */
  PetscCall(PetscFinalize());
  return 0;
}

/*TEST

   build:
      requires: !complex !single

   test:
      args: -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always

   test:
      suffix: 2
      nsize: 3
      args: -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always

   test:
      suffix: 3
      nsize: 2
      args: -ksp_monitor_short -ksp_rtol 1e-6 -ksp_type pipefgmres

TEST*/
