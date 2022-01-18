#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

const TMatrixD multiply_matrix(const TMatrixD& a, const TMatrixD& b) {
  
  //Code to multiply two matrices as follows: Cij = sum_k [ Aik * Bkj ] 

  cout << "Matrix A:" << endl;
  a.Print();
  cout << "Matrix B:" << endl;
  b.Print();
  const int a_ncols = a.GetNcols();
  const int a_nrows = a.GetNrows();
  const int b_ncols = b.GetNcols();
  const int b_nrows = b.GetNrows();

  if(a_ncols != b_nrows)
    {
      cout << Form("The number of columns in Matrix A (%i columns) does NOT match the number of rows in Matrix B (%i rows)", a_ncols, b_nrows ) << endl;
      cout << "Unable to Multiply incompatible dimensions in Matrices. Exiting Now . . ." << endl;
      gSystem->Exit(1);
    }
  
  cout << "Multiplying Matrix AxB: " << endl;

  //Result Matrix TMatrixD
  TMatrixD R(a_nrows,b_ncols);
  R.Zero();

  for(int i=0; i<a_nrows; i++){
    
    for(int j=0; j<b_ncols; j++){

      R[i][j] = 0;

      for(int k=0; k<b_nrows; k++ ){

	R[i][j] += a[i][k]*b[k][j];
	
      }
    
    }

  }

  R.Print();
  

  return R;
}
