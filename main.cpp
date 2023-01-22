#include <iostream>
#include <vector>
 
using namespace std;

//Make a class Matrix that will help us to do matrix math operations
class Matrix {
  // private variables: rows_, cols_, data_ (internal values of the matrix)
  // we use _ suffix to denote inner variables
  // for data we use vector in order to avoid pointer to pointer
  private:
    unsigned rows_, cols_;
    vector<vector<double>> data_;

  public:
    // constructor method

    // syntax sugar, we use inline asignation, intead of manual asignation
    // Matrix(int rows, int cols) {
      // rows_ = rows;
      // cols_ = cols;

    Matrix(unsigned rows, unsigned cols, double initial=0.0): rows_(rows), cols_(cols) {
      // generates a new vector of len rows
      data_.resize(rows);
      for (int i = 0; i < rows; i++) {
        // generates a new vector of len cols for each row
        data_[i].resize(cols, initial);
      }
    }
    
    // desconstructor method
    // clean vector in memory
    // it is important due to memory leaks
    ~Matrix() {
      data_.clear();
    }

    unsigned rows() { return rows_; }
    unsigned cols() { return cols_; }

    // global getter, to obtain the private value in data_[row][col]
    double get(unsigned row, unsigned col) {
      return data_[row][col];
    }

    // global setter, to change the private value in data_[row][col]
    void set(unsigned row, unsigned col, double value) {
      data_[row][col] = value;
    }

    // metho to print matrix
    void display() {
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          std::cout << data_[i][j] << " ";
        }
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    // method to sum self matrix to an external scalar
    // this method will create a new matrix
    Matrix scalarSum(double value) {
      Matrix result(rows_, cols_);

      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j] + value);
        }
      }

      return result;
    }

    // method to substract self matrix to an external scalar
    // this method will create a new matrix
    Matrix scalarSub(double value) {
      Matrix result(rows_, cols_);

      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j] - value);
        }
      }

      return result;
    }

    // method to multiply self matrix to an external scalar
    // this method will create a new matrix
    Matrix scalarDot(double value) {
      Matrix result(rows_, cols_);

      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j] * value);
        }
      }

      return result;
    }

    // method to devide self matrix to an external scalar
    // this method will create a new matrix
    Matrix scalarDiv(double value) {
      Matrix result(rows_, cols_);

      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j] / value);
        }
      }

      return result;
    }

    // method to sum self matrix to an external matrix
    // this method will create a new matrix
    Matrix matrixSum(Matrix &matrix) {
      Matrix result(rows_, cols_);

      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j] + matrix.get(i, j));
        }
      }

      return result;
    }

    // method to substract self matrix to an external matrix
    // this method will create a new matrix
    Matrix matrixSub(Matrix &matrix) {
      Matrix result(rows_, cols_);

      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j] - matrix.get(i, j));
        }
      }

      return result;
    }

    // method to miultiply self matrix to an external matrix
    // this method will create a new matrix
    Matrix matrixDot(Matrix &matrix) {
      Matrix result(rows_, cols_);

      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < matrix.cols(); j++) {
          double tmp = 0.0;
          for (int k = 0; k < cols_; k++) {
            tmp += data_[i][k] * matrix.get(k, j);
          }
          result.set(i, j, tmp);
        }
      }

      return result;
    }

    // method to miultiply self matrix to an external matrix by the Strassen method
    // this method will create a new matrix
    Matrix strassen(Matrix &A, Matrix &B) {
      int n = A.rows();

      // If a unit matrix is entered, then simply return their product 
      if (n == 1) {
        Matrix res(1, 1, A.get(0, 0) * B.get(0, 0));
      }

      int half = n / 2;

      Matrix A11(half, half), A12(half, half), A21(half, half), A22(half, half);

      Matrix B11(half, half), B12(half, half), B21(half, half), B22(half, half);

      for (int i = 0; i < half; i++) {
        for (int j = 0; j < half; j++) {
          A11.set(i, j, A.get(i, j));
          A12.set(i, j, A.get(i, j + half));
          A21.set(i, j, A.get(i + half, j));
          A22.set(i, j, A.get(i + half, j + half));

          B11.set(i, j, B.get(i, j));
          B12.set(i, j, B.get(i, j + half));
          B21.set(i, j, B.get(i + half, j));
          B22.set(i, j, B.get(i + half, j + half));
        }
      }

      Matrix S1 = A11.matrixSum(A22);
      Matrix S2 = B11.matrixSum(B22);
      Matrix P = strassen(S1, S2);

      Matrix S3 = A21.matrixSum(A22);
      Matrix Q = strassen(S3, B11);

      Matrix S4 = B12.matrixSub(B22);
      Matrix R = strassen(A11, S4);
      
      Matrix S5 = B21.matrixSub(B11);
      Matrix S = strassen(A22, S5);

      Matrix S6 = A11.matrixSum(A12);
      Matrix T = strassen(S6, B22);

      Matrix S7 = A21.matrixSub(A11);
      Matrix S8 = B11.matrixSum(B12);
      Matrix U = strassen(S7, S8);

      Matrix S9 = A12.matrixSub(A22);
      Matrix S10 = B21.matrixSum(B22);
      Matrix V = strassen(S9, S10);


      // to avoid extra code, use class concatenation
      // Matrix Tmp1 = P.matrixSum(S);
      // Matrix Tmp2 = Tmp1.matrixSub(T);
      // Matrix C11 = Temp2.matrixSum(V);
      Matrix C11 = P.matrixSum(S).matrixSub(T).matrixSum(V);
      Matrix C12 = R.matrixSum(T);
      Matrix C21 = Q.matrixSum(S);
      Matrix C22 = P.matrixSum(R).matrixSub(Q).matrixSum(U);

      Matrix Result(n, n);
      
      for (int i = 0; i < half; i++) {
        for (int j = 0; j < half; j++) {
          Result.set(i, j, C11.get(i, j));
          Result.set(i, j + half, C12.get(i, j));
          Result.set(i + half, j, C21.get(i, j));
          Result.set(i + half, j + half, C22.get(i, j));
        }
      }

      return Result;
    }

    // private method
    // this method will create a new matrix
    Matrix cofactor(int p, int q) {
      Matrix result(rows_-1, cols_-1);

      int i = 0, j = 0;
  
      for (int row = 0; row < rows_; row++) {
        for (int col = 0; col < cols_; col++) {
          if (row != p && col != q) {
            result.set(i, j++, data_[row][col]);

            if (j == rows_ - 1) {
              j = 0;
              i++;
            }
          }
        }
      }

      return result;
    }

    // method to find the determinant of the matrix
    // this method will create a new matrix
    double determinant(int n) {
      int D = 0; // Initialize result

      //  Base case : if matrix contains single element
      if (n == 1) {
        return data_[0][0];
      }

      Matrix result(rows_, cols_);

      // Looping for each element of the matrix
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j]);
        }
      }

      int sign = 1; // To store sign multiplier

      // Iterate for each element of first row
      for (int f = 0; f < n; f++) {
        // Getting Cofactor of A[0][f]
        Matrix cofact = result.cofactor(0, f);
        D += sign * result.get(0, f) * cofact.determinant(n - 1);

        // terms are to be added with alternate sign
        sign = -sign;
      }

      return D;
    }

    // method to find the transpose of the matrix
    // this method will create a new matrix
    Matrix transpose() {
      Matrix result(cols_, rows_);
      
      for (int i = 0; i < cols_; i++) {
        for (int j = 0; j < rows_; j++) {
          result.set(i, j, data_[j][i]);
        }
      }

      return result;
    }

    // method to find the adjoint of the matrix
    // this method will create a new matrix
    Matrix adjoint() {
      Matrix result(rows_, cols_);

      if (rows_ == 1 && cols_ == 1) {
        result.set(0, 0, 1);
        return result;
      }

      Matrix temp(rows_, cols_);

      // Looping for each element of the matrix
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j]);
          temp.set(i, j, data_[i][j]);
        }
      }

      int sign = 1;
      int n = result.rows();

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Matrix cof = temp.cofactor(i, j);
          sign = ((i + j) % 2 == 0) ? 1 : -1;
          int det = cof.determinant(n - 1);
          result.set(j, i, sign * det);
        }
      }

      return result;
    }

    // method to miultiply self matrix to an external matrix
    // this method will create a new matrix
    Matrix inverse() {
      Matrix result(rows_, cols_);

      // Looping for each element of the matrix
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, data_[i][j]);
        }
      }
      
      // Find determinant of the given Matrix
      int det = result.determinant(result.cols());
      if (det == 0) {
        cout << "Unable to find the inverse of a singular matrix. Returning -1 error code." << endl;
        Matrix error(1, 1, -1);
        return error;
      }
  
      // Find adjoint
      Matrix adjacent = result.adjoint();
  
      // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
      for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_; j++) {
          result.set(i, j, adjacent.get(i, j) / double(det));
        }
      }
  
      return result;
    }

    // method to find the pseiudo inverse matrix
    // this method will create a new matrix
    Matrix psdinverse() {
      if(cols_ >= rows_) {
        cout << "cols_ ≥ rows_" << endl;
        Matrix origin(rows_, cols_);

        for (int i = 0; i < rows_; i++) {
          for (int j = 0; j < cols_; j++) {
            origin.set(i, j, data_[i][j]);
          }
        }

        Matrix Transp = origin.transpose();
        Matrix RightInvertible = (origin.matrixDot(Transp));
        Matrix RI_Inv(rows_, cols_);
        Matrix Pseudo(rows_, cols_);

        RI_Inv = RightInvertible.inverse();
        Pseudo = RightInvertible.matrixDot(RI_Inv);

        return Pseudo;
      } else if(cols_ <= rows_) {
        Matrix pending(1, 1);
        return pending; // pending to develop this part
      } else {
        Matrix error(1, 1, -1);
        cout << "No valid alternative found to solve this problem. Returning error code -1" << endl;
        return error;
      }
    }
};

int main() {
  int N = 4;

  Matrix A(N, N);
  int A_data[N][N] = { {1, 2, 3, 4}, {5, 6, 7, 8}, {9, 0, 1, 2}, {3, 4, 5, 6} };
  for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        A.set(i, j, A_data[i][j]);
      }
    }

  Matrix B(N, N);
  int B_data[N][N] = { {0, 9, 8, 7}, {6, 5, 4, 3}, {2, 1, 0, 9}, {8, 7, 6, 5} };
  for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        B.set(i, j, B_data[i][j]);
      }
    }

    
  std::cout << "Matrix A: " << std::endl;
  A.display();

  std::cout << "Matrix B: " << std::endl;
  B.display();

  Matrix scalarSum = A.scalarSum(5.0);
  std::cout << "Scalar sum: " << std::endl;
  scalarSum.display();

  Matrix scalarSub = A.scalarSub(5.0);
  std::cout << "Scalar sub: " << std::endl;
  scalarSub.display();

  Matrix scalarDot = A.scalarDot(5.0);
  std::cout << "Scalar dot: " << std::endl;
  scalarDot.display();

  Matrix scalarDiv = A.scalarDiv(5.0);
  std::cout << "Scalar div: " << std::endl;
  scalarDiv.display();

  Matrix sum = A.matrixSum(B);
  std::cout << "Matrix sum: " << std::endl;
  sum.display();

  Matrix subs = A.matrixSub(B);
  std::cout << "Matrix subs: " << std::endl;
  subs.display();

  Matrix dot = A.matrixDot(B);
  std::cout << "Matrix dot: " << std::endl;
  dot.display();

  // Matrix strassen = A.strassen(A, B);
  // std::cout << "Matrix strassen: " << std::endl;
  // strassen.display();

  double determinant = A.determinant(N);
  std::cout << "Matrix determinant: " << determinant << std::endl << std::endl;

  Matrix transpose = A.transpose();
  std::cout << "Matrix transpose: " << std::endl;
  transpose.display();

  Matrix inverse = A.inverse();
  std::cout << "Matrix inverse: " << std::endl;
  inverse.display();

  Matrix adjoint = A.adjoint();
  std::cout << "Matrix adjoint: " << std::endl;
  adjoint.display();

  return 0;
}