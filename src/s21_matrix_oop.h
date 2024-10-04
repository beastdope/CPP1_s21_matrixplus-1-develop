#ifndef S21_MATRIX_HEADER
#define SRC_S21_MATRIX_H_

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

class S21Matrix {
 private:
  // attributes
  int rows_, cols_;  // rows and columns
  double *matrix_;   // pointer to the memory where the matrix is allocated

 public:
  S21Matrix();  // default constructor
  S21Matrix(int rows, int cols);
  S21Matrix(const double matrix[], int rows, int cols);
  S21Matrix(const S21Matrix &other);  // copy constructor
  S21Matrix(S21Matrix &&other);       // move constructor
  ~S21Matrix();                       // destructor

  // operators
  S21Matrix &operator+=(const S21Matrix &other);
  S21Matrix &operator-=(const S21Matrix &other);
  S21Matrix &operator*=(const S21Matrix &other);
  S21Matrix &operator*=(double other);
  S21Matrix &operator=(const S21Matrix &other);
  S21Matrix &operator=(S21Matrix &&other) noexcept;

  S21Matrix operator+(const S21Matrix &other);
  S21Matrix operator-(const S21Matrix &other);
  S21Matrix operator*(const S21Matrix &other);
  S21Matrix operator*(double other);

  bool operator==(const S21Matrix &other) const;
  bool operator!=(const S21Matrix &other) const;

  double &operator()(int row, int col) const;

  double &at(int row, int col) const;

  void SumMatrix(const S21Matrix &other);
  bool EqMatrix(const S21Matrix &other) const;
  void SubMatrix(const S21Matrix &other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix &other);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  double Determinant() const;
  double Minor(int row, int col) const;
  S21Matrix InverseMatrix() const;
  void RemoveRowCol(const double *matrix, int matrix_size, double *target,
                    int row_skip, int col_skip) const;
  double Det(const double *mat, int e) const;
  int GetRows() const;
  int GetCols() const;
  void SetRows(int new_rows);
  void SetCols(int new_col);
  double *GetMatrix() const;

  void SetElement(int row, int col, double value);
  float GetElement(int row, int col) const;
};

#endif