#include "mmloader.h"

PetscErrorCode MatCreateFromMTX(Mat *A, const char *filein, PetscBool aijonly, PetscInt *m, PetscInt *n)
{
  MM_typecode matcode;
  FILE        *file;

  PetscInt    M, N, NZ;                     // M: Number of rows, N: Number of cols, NZ: Number of non zeros.
  PetscInt    nz = 0;                       // nz: local number of non zeros.
  PetscInt    *ia, *ja;                     // Indices of rows, Indices of cols.
  PetscInt    i, j, ninput, *rownz;
  PetscMPIInt starting_row = 0, starting_nz = 0;

  PetscScalar *val;
  PetscBool   sametype, symmetric = PETSC_FALSE, skew = PETSC_FALSE;
  PetscMPIInt size, rank;

  /* Read in matrix */
  PetscFunctionBeginUser;

  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  PetscCall(PetscFOpen(PETSC_COMM_SELF, filein, "r", &file));
  PetscCheck(mm_read_banner(file, &matcode) == 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Could not process Matrix Market banner.");
  /*  This is how one can screen matrix types if their application */
  /*  only supports a subset of the Matrix Market data types.      */
  PetscCheck(mm_is_matrix(matcode) && mm_is_sparse(matcode) && (mm_is_real(matcode) || mm_is_integer(matcode)), PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Input must be a sparse real or integer matrix. Market Market type: [%s]", mm_typecode_to_str(matcode));

  if (mm_is_symmetric(matcode)) symmetric = PETSC_TRUE;
  if (mm_is_skew(matcode)) skew = PETSC_TRUE;

  /* Find out size of sparse matrix .... */
  PetscCheck(mm_read_mtx_crd_size(file, &M, &N, &NZ) == 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Size of sparse matrix is wrong.");

  /* Find out local number of rows. */
  PetscCall(PetscSplitOwnership(PETSC_COMM_WORLD, m, &M));
  PetscCall(PetscSplitOwnership(PETSC_COMM_WORLD, n, &N));

  /* Starting row = Sum of previous local n, myself excluded. */
  PetscCallMPI(MPI_Exscan(m, &starting_row, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD));

  /* Reserve memory for matrices, can be more restricted. */
  PetscCall(PetscMalloc4(NZ, &ia, NZ, &ja, NZ, &val, M, &rownz));
  for (i = 0; i < M; i++) rownz[i] = 0;

  /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
  /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
  /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

  for (i = 0; i < NZ; i++) {
    ninput = fscanf(file, "%d %d %lg\n", &ia[nz], &ja[nz], &val[nz]);
    PetscCheck(ninput >= 3, PETSC_COMM_WORLD, PETSC_ERR_FILE_UNEXPECTED, "Badly formatted input file");
    ia[nz]--;
    ja[nz]--;                                 /* adjust from 1-based to 0-based */

    // Read only in my work zone.
    if ((ia[nz] < starting_row) || ((starting_row + *m) <= ia[nz])) continue;

    if ((symmetric && aijonly) || skew) {     /* transpose */
      rownz[ia[nz] - starting_row]++;
      if (ja[nz] != ia[nz]) {
        rownz[ja[nz] - starting_row]++;
      }
    } else {
      if (symmetric) {
        rownz[ja[nz] - starting_row]++;
      }
      else {
        rownz[ia[nz] - starting_row]++;
      }
    }
    nz++;                                     /* Incrementing the local number of non zeros. */
  }

  PetscCallMPI(MPI_Exscan(&nz, &starting_nz, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD));
  PetscCall(PetscFClose(PETSC_COMM_SELF, file));

  /* Create, preallocate, and then assemble the matrix */
  PetscCall(MatCreate(PETSC_COMM_WORLD, A));
  PetscCall(MatSetSizes(*A, *m, PETSC_DECIDE, M, N));

  if (symmetric && !aijonly) {
    PetscCall(MatSetType(*A, MATSBAIJ));
    PetscCall(MatSetFromOptions(*A));
    PetscCall(PetscObjectTypeCompareAny((PetscObject)*A, &sametype, MATSEQSBAIJ, MATMPISBAIJ, ""));
    PetscCall(PetscObjectTypeCompare((PetscObject)*A, MATSBAIJ, &sametype));
    PetscCheck(sametype, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Only AIJ and SBAIJ are supported. Your mattype is not supported");
  } else {
    PetscCall(MatSetType(*A, MATAIJ));
    PetscCall(MatSetFromOptions(*A));
    PetscCall(PetscObjectTypeCompareAny((PetscObject)*A, &sametype, MATSEQAIJ, MATMPIAIJ, ""));
    PetscCheck(sametype, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Only AIJ and SBAIJ are supported. Your mattype is not supported");
  }
  /* Add values to the matrix, these correspond to lower triangular part for symmetric or skew matrices */
  if (!(symmetric && !aijonly)){
      for (j = 0; j < nz; j++) {
        PetscCall(MatSetValues(*A, 1, &ia[j], 1, &ja[j], &val[j], INSERT_VALUES));
      }
    }

  /* Add values to upper triangular part for some cases */
  if (symmetric && !aijonly) {
    /* MatrixMarket matrix stores symm matrix in lower triangular part. Take its transpose */
    for (j = 0; j < nz; j++) {
      PetscCall(MatSetValues(*A, 1, &ja[j], 1, &ia[j], &val[j], INSERT_VALUES));
    }
  }
  if (skew) {
    for (j = 0; j < nz; j++) {
      val[j] = -val[j];
      PetscCall(MatSetValues(*A, 1, &ja[j], 1, &ia[j], &val[j], INSERT_VALUES));
    }
  }
  PetscCall(MatAssemblyBegin(*A, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(*A, MAT_FINAL_ASSEMBLY));
  PetscCall(PetscFree4(ia, ja, val, rownz));
  PetscFunctionReturn(PETSC_SUCCESS);
}
