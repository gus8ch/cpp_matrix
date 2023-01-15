#include <iostream>

class Matrix {
  // private variables: rows_, cols_, data_
  private:
    int rows_, cols_;
    int** data_;

  public:
    // constructor method

    // syntax sugar, we use inline asignation, intead of manual asignation
    // Matrix(int rows, int cols) {
      // rows_ = rows;
      // cols_ = cols;

    Matrix(int rows, int cols): rows_(rows), cols_(cols) {
      // generates a new pointer to array of len rows
      data_ = new int*[rows];
      for (int i = 0; i < rows; i++) {
        // generates a new array of len clos for each cols
        data_[i] = new int[cols];
      }
    }
    // desconstructor method
    ~Matrix() {
      for (int i = 0; i < rows_; i++) {
        delete[] data_[i];
      }
      delete[] data_;
    }

  int rows() { return rows_; }
  int cols() { return cols_; }

  // global getter, to obtain the private value in data_[row][col]
  int get(int row, int col) {
    return data_[row][col];
  }

  // global setter, to change the private value in data_[row][col]
  void set(int row, int col, int value) {
    data_[row][col] = value;
  }

  // setter, to change values from array
  // void set_from_array(int data[][3]) {
  //   for (int i = 0; i < rows_; i++) {
  //     for (int j = 0; j < cols_; j++) {
  //       data_[i][j] = data[i][j];
  //     }
  //   }
  // }

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

  void sum(Matrix &matrix) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        int tmp = data_[i][j] + matrix.get(i, j);
        std::cout << tmp << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  void dot(Matrix &matrix) {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        int tmp = data_[i][j] * matrix.get(j, i);
        std::cout << tmp << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
};

int main() {
  Matrix A(3, 3);
  int A_data[3][3] = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
  // A.set_from_array(A_data);
  for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        A.set(i, j, A_data[i][j]);
      }
    }
  A.display();

  Matrix B(3, 3);
  int B_data[3][3] = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
  // A.set_from_array(B_data);
  for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        B.set(i, j, B_data[i][j]);
      }
    }
  B.display();

  A.sum(B);

  A.dot(B);

  return 0;
}