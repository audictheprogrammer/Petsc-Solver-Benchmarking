static char help[] = "Loads a MTX matrix then solves a linear system with KSP.\n\n";

#include <petsc.h>
#include "mmloader.h"

/*
export wPETSC_DIR=/mnt/c/Users/audic/petsc/share/petsc/datafiles/matrices/MYMAT
export wPETSC_DIR=/mnt/c/Users/xu/petsc/share/petsc/datafiles/matrices/MYMAT
./ex72 -fin ${wPETSC_DIR}/1138_bus.mtx petscmat.aij -aij_only -pc_type lu -pc_factor_mat_solver_type mumps -memory_view -log_view -ksp_view
./ex72 -fin ${wPETSC_DIR}/cvxbqp1.mtx petscmat.aij -aij_only -pc_type lu -pc_factor_mat_solver_type mumps -memory_view -log_view -ksp_view> TEST1.txt
*/

int main(int argc, char **args)
{
  MM_typecode matcode;
  FILE       *file;

  Vec         x, b, u; /* approx solution, RHS, exact solution */
  Mat         A;       /* linear system matrix */
  KSP         ksp;     /* linear solver context */
  PC          pc;      /* preconditioner context */
  PetscReal   norm;    /* norm of solution error */
  PetscInt    n, its;

  PetscInt    M, N, nz;
  char        filein[PETSC_MAX_PATH_LEN];
  char        ordering[256] = MATORDERINGRCM;
  // PetscViewer view;
  PetscBool   flag, aijonly = PETSC_FALSE, permute = PETSC_FALSE;
  IS          rowperm = NULL, colperm = NULL;

  PetscMPIInt size;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &args, (char *)0, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "This is a uniprocessor example only!");

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

  n = N;

   /* Matrix loaded. */

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Ax = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create vectors.  Note that we form 1 vector from scratch and
     then duplicate as needed.
  */
  PetscCall(VecCreate(PETSC_COMM_SELF, &x));
  PetscCall(PetscObjectSetName((PetscObject)x, "Solution"));
  PetscCall(VecSetSizes(x, PETSC_DECIDE, n));
  PetscCall(VecSetFromOptions(x));
  PetscCall(VecDuplicate(x, &b));
  PetscCall(VecDuplicate(x, &u));

  /*
     Set exact solution; then compute right-hand-side vector.
  */
  PetscCall(VecSet(u, 1.0));
  PetscCall(MatMult(A, u, b));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the linear solver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));

  /*
     Set operators. Here the matrix that defines the linear system
     also serves as the matrix that defines the preconditioner.
  */
  PetscCall(KSPSetOperators(ksp, A, A));

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

  /* check that KSP automatically handles the fact that the the new non-zero values in the matrix are propagated to the KSP solver */
  // PetscCall(MatShift(A, 2.0));

  PetscCall(KSPSolve(ksp, b, x));

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check the solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  PetscCall(VecAXPY(x, -1.0, u));
  PetscCall(VecNorm(x, NORM_2, &norm));
  PetscCall(KSPGetIterationNumber(ksp, &its));
  PetscCall(PetscPrintf(PETSC_COMM_SELF, "Norm of error %g, Iterations %" PetscInt_FMT "\n", (double)norm, its));

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  PetscCall(KSPDestroy(&ksp));

  /* test if prefixes properly propagate to PCMPI objects */
  if (PCMPIServerActive) {
    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
    PetscCall(KSPSetOptionsPrefix(ksp, "prefix_test_"));
    PetscCall(MatSetOptionsPrefix(A, "prefix_test_"));
    PetscCall(KSPSetOperators(ksp, A, A));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp, b, x));
    PetscCall(KSPDestroy(&ksp));
  }

  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&u));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));

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

   test:
      args: -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always

   test:
      suffix: 2
      args: -pc_type sor -pc_sor_symmetric -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always

   test:
      suffix: 2_aijcusparse
      requires: cuda
      args: -pc_type sor -pc_sor_symmetric -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always -mat_type aijcusparse -vec_type cuda
      args: -ksp_view

   test:
      suffix: 3
      args: -pc_type eisenstat -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always

   test:
      suffix: 3_aijcusparse
      requires: cuda
      args: -pc_type eisenstat -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always -mat_type aijcusparse -vec_type cuda -ksp_view

   test:
      suffix: aijcusparse
      requires: cuda
      args: -ksp_monitor_short -ksp_gmres_cgs_refinement_type refine_always -mat_type aijcusparse -vec_type cuda -ksp_view
      output_file: output/ex1_1_aijcusparse.out

   test:
      requires: defined(PETSC_USE_SINGLE_LIBRARY)
      suffix: mpi_linear_solver_server_1
      nsize: 3
      filter: sed 's?ATOL?RTOL?g' | grep -v HERMITIAN
      # use the MPI Linear Solver Server
      args: -mpi_linear_solver_server -mpi_linear_solver_server_view
      # controls for the use of PCMPI on a particular system
      args: -mpi_linear_solver_server_minimum_count_per_rank 5 -mpi_linear_solver_server_ksp_view -mpi_linear_solver_server_mat_view
      # the usual options for the linear solver (in this case using the server)
      args: -ksp_monitor -ksp_converged_reason -mat_view -ksp_view  -ksp_type cg -pc_type none
      # the options for the prefixed objects
      args: -prefix_test_mpi_linear_solver_server_mat_view -prefix_test_ksp_monitor -prefix_test_mpi_linear_solver_server_minimum_count_per_rank 5

   test:
      requires: !__float128
      suffix: minit
      args: -ksp_monitor -pc_type none -ksp_min_it 8

TEST*/


/*TEST

   build:
      requires:  !complex double !defined(PETSC_USE_64BIT_INDICES)
      depends: mmloader.c mmio.c

   test:
      suffix: 1
      args: -fin ${wPETSC_DIR}/share/petsc/datafiles/matrices/amesos2_test_mat0.mtx

   test:
      suffix: 2
      args: -fin ${wPETSC_DIR}/share/petsc/datafiles/matrices/LFAT5.mtx petscmat.sbaij

   test:
      suffix: 3
      args: -fin ${wPETSC_DIR}/share/petsc/datafiles/matrices/m_05_05_crk.mtx petscmat2.aij

   test:
      suffix: 4
      args: -fin ${wPETSC_DIR}/share/petsc/datafiles/matrices/amesos2_test_mat0.mtx petscmat.aij -permute rcm
TEST*/
