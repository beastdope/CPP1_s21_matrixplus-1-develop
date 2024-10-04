#include <gtest/gtest.h>

#include "../s21_matrix_oop.h"
constexpr double kEpsilon = 1e-7;
/**
 * Функция matrixCheck выполняет итерацию по матрице и проверяет, равен ли
 * каждый элемент примерно заданному значению в пределах указанного допуска
 * эпсилон.
 *
 * @param matrix Параметр «matrix» имеет тип «S21Matrix», который является
 * пользовательским классом матрицы. Функция CheckMatrix выполняет итерацию по
 * каждому элементу матрицы, чтобы проверить, равно ли значение в этой позиции
 * примерно значению параметра value в пределах допуска kEpsilon.
 * @param value Параметр value в функции CheckMatrix представляет ожидаемое
 * значение, с которым следует сравнивать каждый элемент матрицы. Функция
 * перебирает каждый элемент матрицы и проверяет, приблизительно ли значение
 * этого элемента равно указанному «значению» в пределах допуска, определенного
 * параметром «
 */

bool almost_equal(double a, double b) { return std::fabs(a - b) < kEpsilon; }

bool matrixEqual(const S21Matrix &matrix, const double *otherMatrix) {
  const double *matrix1 = matrix.GetMatrix();
  long element = static_cast<long>(matrix.GetRows()) * matrix.GetCols();

  for (long i = 0; i < element; ++i) {
    if (!almost_equal(matrix1[i], otherMatrix[i])) {
      delete[] matrix1;
    }
  }
  return true;
}

void matrixCheck(const S21Matrix &matrix, double value) {
  for (int i = 0; i < matrix.GetRows(); ++i) {
    for (int j = 0; j < matrix.GetCols(); ++j) {
      ASSERT_NEAR(matrix(i, j), value, kEpsilon);
    }
  }
}

void FillMatrix(S21Matrix &matrix, double value) {
  for (int i = 0; i < matrix.GetRows(); ++i) {
    for (int j = 0; j < matrix.GetCols(); ++j) {
      matrix(i, j) = value;
    }
  }
}

const double EPSILON = 1e-7;
const double identity_2_by_3[6] = {1, 0, 0, 0, 1, 0};

/**
 * Функция тестирует конструкцию матричного объекта и проверяет его исходное
 * состояние.
 *
 * @param  Похоже, вы пишете тестовый пример для класса S21Matrix. В этом
 * тестовом примере вы проверяете исходное состояние матричного объекта после
 * построения. Тест проверяет, что количество строк и столбцов равно 0, а
 * указатель данных матрицы имеет значение nullptr.
 */
TEST(S21Matrix, construct) {
  S21Matrix matrix(3, 3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      matrix(i, j) = i * 3 + j;
    }
  }

  EXPECT_EQ(matrix(0, 0), 0);
  EXPECT_EQ(matrix(1, 1), 4);
  EXPECT_EQ(matrix(2, 2), 8);

  EXPECT_THROW(matrix(-1, 0), std::out_of_range);
  EXPECT_THROW(matrix(0, 3), std::out_of_range);
  EXPECT_THROW(matrix(3, 0), std::out_of_range);
}
TEST(S21Matrix, constructException) {
  EXPECT_ANY_THROW(S21Matrix matrix(-1, 2));
  EXPECT_ANY_THROW(S21Matrix matrix(2, -3));
}
TEST(S21Matrix, constructExceptionZero) {
  EXPECT_ANY_THROW(S21Matrix matrix(0, 2));
  EXPECT_ANY_THROW(S21Matrix matrix(2, 0));
}

TEST(S21Matrix, constructRowCorrect) {
  S21Matrix matrix;
  matrix.SetRows(1);
  EXPECT_TRUE(matrix.GetRows() == 1);
}
TEST(S21Matrix, constructRows) {
  S21Matrix matrix(2, 3);
  S21Matrix matrix_last = matrix;
  EXPECT_ANY_THROW(matrix_last.SetRows(-1));
  EXPECT_TRUE(matrix_last == matrix);
}
TEST(S21Matrix, constructCols) {
  S21Matrix matrix(2, 3);
  S21Matrix matrix_last = matrix;
  EXPECT_ANY_THROW(matrix_last.SetCols(-1));
  EXPECT_TRUE(matrix_last == matrix);
}

TEST(S21Matrix, initialiseParametersNullPtr) {
  EXPECT_ANY_THROW(S21Matrix matrix(nullptr, 2, 3));
}
TEST(S21Matrix, initialiseParametersException) {
  EXPECT_ANY_THROW(S21Matrix matrix(identity_2_by_3, 0, 3));
}
TEST(S21Matrix, initialiseParametersExceptionLessZero) {
  EXPECT_ANY_THROW(S21Matrix matrix(identity_2_by_3, -1, 3));
}

TEST(S21Matrix, initialiseParameters) {
  S21Matrix matrix(2, 3);
  EXPECT_EQ(matrix.GetRows(), 2);
  EXPECT_EQ(matrix.GetCols(), 3);
}
TEST(S21Matrix, initialiseParameters2) {
  int rows = 3;
  int cols = 2;
  double matrixIdent[rows * cols]{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  S21Matrix matrix(matrixIdent, rows, cols);

  EXPECT_EQ(matrix.GetRows(), rows);
  EXPECT_EQ(matrix.GetCols(), cols);
  double *matrixData = matrix.GetMatrix();
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      EXPECT_EQ(matrixData[i * cols + j], matrixIdent[i * cols + j]);
    }
  }
}
TEST(S21Matrix, matrixMove) {
  S21Matrix matrix;

  matrix.SetCols(2);
  matrix.SetRows(3);
  EXPECT_EQ(matrix.GetCols(), 2);
  EXPECT_EQ(matrix.GetRows(), 3);

  S21Matrix matrixCopy(matrix);
  EXPECT_EQ(matrixCopy.GetCols(), matrix.GetCols());
  EXPECT_EQ(matrixCopy.GetRows(), matrix.GetRows());
  EXPECT_TRUE(matrixCopy == matrix);

  S21Matrix matrixMoved(std::move(matrix));

  EXPECT_EQ(matrixMoved.GetCols(), 2);
  EXPECT_EQ(matrixMoved.GetRows(), 3);

  EXPECT_TRUE(matrixMoved == matrixCopy);
}

TEST(S21Matrix, copy) {
  S21Matrix matrix(2, 2);
  S21Matrix matrixCopy;
  matrix = matrixCopy;
  EXPECT_TRUE(matrixCopy == matrix);
}
TEST(S21Matrix, matrixSumException) {
  S21Matrix matrix(1, 2);
  S21Matrix matrixCopy(2, 3);

  S21Matrix matrixToWork = matrix;
  S21Matrix matrixToWork2 = matrixCopy;

  EXPECT_ANY_THROW(matrixToWork.SumMatrix(matrixToWork2));
  EXPECT_ANY_THROW(matrixToWork2.SumMatrix(matrixToWork));
  EXPECT_TRUE(matrix == matrixToWork);
  EXPECT_TRUE(matrixCopy == matrixToWork2);
}
/**
 * Функция определяет тестовый пример для операции суммирования матриц в C++ с
 использованием
 * специального класса матрицы.
 *
 * @param  Похоже, вы пишете тестовый пример для матричного класса S21Matrix.
 В этом тестовом примере
 * вы создаете две матрицы, устанавливаете их размеры 2x2, а затем
 инициализируете их значения
 * определенными числами.
 */
TEST(S21Matrix, matrixSum1) {
  S21Matrix matrix(2, 3);
  S21Matrix matrix2(2, 3);
  matrix.SumMatrix(matrix2);
  double matrixSum[6] = {2, 0, 0, 0, 2, 0};
  EXPECT_TRUE(matrixEqual(matrix, matrixSum));
  S21Matrix matrix3(3, 3);
  EXPECT_THROW(matrix.SumMatrix(matrix3), std::invalid_argument);
}

TEST(S21MatrixOperator, plusOperator) {
  double *matrix1_ptr = new double[4]{1, 2, 3, 4};
  double *matrix2_ptr = new double[4]{41, 42, 43, 44};
  S21Matrix matrix1(matrix1_ptr, 2, 2);
  S21Matrix matrix2(matrix2_ptr, 2, 2);
  S21Matrix matrix3 = matrix1 + matrix2;
  double other[4] = {42, 44, 46, 48};
  EXPECT_TRUE(matrixEqual(matrix3, other));
  matrix1 += matrix2;
  EXPECT_TRUE(matrixEqual(matrix1, other));
  matrix1.SetRows(3);
  EXPECT_THROW({ matrix1 += matrix2; }, std::invalid_argument);
  EXPECT_THROW({ S21Matrix matrix4 = matrix1 + matrix2; },
               std::invalid_argument);
  delete[] matrix1_ptr;
  delete[] matrix2_ptr;
}
TEST(S21MatrixOperator, minusOperator) {
  double *matrix1_ptr = new double[4]{41, 42, 43, 44};
  double *matrix2_ptr = new double[4]{1, 2, 3, 4};
  S21Matrix matrix1(matrix1_ptr, 2, 2);
  S21Matrix matrix2(matrix2_ptr, 2, 2);
  S21Matrix matrix3 = matrix1 - matrix2;
  double other[4] = {40, 40, 40, 40};
  EXPECT_TRUE(matrixEqual(matrix3, other));
  matrix1 -= matrix2;
  EXPECT_TRUE(matrixEqual(matrix1, other));
  matrix1.SetRows(3);
  EXPECT_THROW({ matrix1 -= matrix2; }, std::invalid_argument);
  EXPECT_THROW({ S21Matrix matrix4 = matrix1 - matrix2; },
               std::invalid_argument);
  delete[] matrix1_ptr;
  delete[] matrix2_ptr;
}
TEST(S21MatrixOperator, mulOperator) {
  double *matrix_initialise = new double[4]{1, -2, 3, -4};

  S21Matrix matrix1(matrix_initialise, 2, 2);

  S21Matrix matrix2 = matrix1 * 2;
  double matrixSum[4] = {2, -4, 6, -8};
  EXPECT_TRUE(matrixEqual(matrix2, matrixSum));

  matrix1 *= 2;
  EXPECT_TRUE(matrixEqual(matrix1, matrixSum));
  delete[] matrix_initialise;
}

/*умножаются элементы матрицы верхней строчки matrix_initialise1 на столбец
 * matrix_initialise2 схема выглядит так:
 * matrix_initialise1 = |2, 3, 4| ||  matrix_initialise2 = |2, 3|
 * matrix_initialise1 = |5, 6, 7| ||  matrix_initialise2 = |15, 11|
 *                                ||  matrix_initialise2 = |25, 0|
 * (2 * 2) + (3 * 15) + (4 * 25) = 149
 * (2 * 3) + (3 * 11) + (4 * 0) = 39
 *
 * (5 * 2) + (6 * 15) + (7 * 25) = 275
 * (5 * 3) + (6 * 11) + (7 * 0) = 81
//  */
TEST(S21MatrixTest, MulMatrixTest) {
  S21Matrix mat1(2, 3);
  mat1.SetElement(0, 0, 1.0);
  mat1.SetElement(0, 1, 2.0);
  mat1.SetElement(0, 2, 3.0);
  mat1.SetElement(1, 0, 4.0);
  mat1.SetElement(1, 1, 5.0);
  mat1.SetElement(1, 2, 6.0);

  S21Matrix mat2(3, 2);
  mat2.SetElement(0, 0, 7.0);
  mat2.SetElement(0, 1, 8.0);
  mat2.SetElement(1, 0, 9.0);
  mat2.SetElement(1, 1, 10.0);
  mat2.SetElement(2, 0, 11.0);
  mat2.SetElement(2, 1, 12.0);

  S21Matrix expected_result(2, 2);
  expected_result.SetElement(0, 0, 58.0);
  expected_result.SetElement(0, 1, 64.0);
  expected_result.SetElement(1, 0, 139.0);
  expected_result.SetElement(1, 1, 154.0);

  mat1.MulMatrix(mat2);
  for (int i = 0; i < expected_result.GetRows(); ++i) {
    for (int j = 0; j < expected_result.GetCols(); ++j) {
      ASSERT_DOUBLE_EQ(mat1.GetElement(i, j), expected_result.GetElement(i, j));
    }
  }
}

TEST(S21Matrix, multiOperator) {
  double *matrix_initialise1 = new double[4]{1, 2, 3, 4};
  double *matrix_initialise2 = new double[4]{-3, 4, -5, 6};

  S21Matrix matrix1(matrix_initialise1, 2, 2);
  S21Matrix matrix2(matrix_initialise2, 2, 2);

  S21Matrix matrix3 = matrix1 * matrix2;
  double matrixSum[4] = {-13, 16, -29, 36};
  EXPECT_TRUE(matrixEqual(matrix3, matrixSum));

  matrix1 *= matrix2;
  EXPECT_TRUE(matrixEqual(matrix1, matrixSum));

  matrix1.SetCols(3);
  EXPECT_THROW({ matrix1 *= matrix2; }, std::invalid_argument);
  EXPECT_THROW({ S21Matrix matrix4 = matrix1 * matrix2; },
               std::invalid_argument);

  delete[] matrix_initialise1;
  delete[] matrix_initialise2;
}

TEST(S21Matrix, matrixNotEqual) {
  S21Matrix matrix(1, 2);
  S21Matrix matrixCopy(2, 3);
  S21Matrix matrixThird(2, 3);

  S21Matrix matrixToWork = matrix;
  S21Matrix matrixToWork2 = matrixCopy;
  S21Matrix matrixToWork3 = matrixThird;

  EXPECT_ANY_THROW(matrixToWork.SumMatrix(matrixToWork2));
  EXPECT_ANY_THROW(matrixToWork2.SumMatrix(matrixToWork));

  EXPECT_TRUE(matrix != matrixToWork3);
}

TEST(S21Matrix, matrixExcept) {
  S21Matrix matrix(3, 3);
  EXPECT_THROW(matrix(-1, 0), std::out_of_range);
}
TEST(S21Matrix, matrixAt) {
  S21Matrix matrix(3, 3);
  EXPECT_THROW(matrix.at(-1, 0), std::out_of_range);
}

// проверка неравных элементов матрицы
TEST(S21Matrix, matrixNotEqualExcept) {
  S21Matrix matrix(2, 2);
  S21Matrix matrixCopy(2, 2);
  matrix.SetElement(0, 0, 1);
  matrix.SetElement(0, 1, 1);
  matrix.SetElement(1, 0, 1);
  matrix.SetElement(1, 1, 1);

  matrixCopy.SetElement(0, 0, 1);
  matrixCopy.SetElement(0, 1, 1);
  matrixCopy.SetElement(1, 0, 2);
  matrixCopy.SetElement(1, 1, 1);

  EXPECT_FALSE(matrix.EqMatrix(matrixCopy));
}
/** transpose делает из строчной матрицы столбцовую например:
 * | 1, 2, 3 | -> | 1, 4 |
 * | 4, 5, 6 | -> | 2, 5 |
 *                | 3, 6 |
 */
TEST(S21Matrix, transpose) {
  S21Matrix matrix1(2, 2);
  EXPECT_EQ(matrix1.GetRows(), 2);
  EXPECT_EQ(matrix1.GetCols(), 2);

  matrix1.SetElement(0, 0, 1.1);
  matrix1.SetElement(1, 0, 1.2);
  matrix1.SetElement(0, 1, 1.3);
  matrix1.SetElement(1, 1, 1.4);

  S21Matrix result_check{2, 2};
  result_check.SetElement(0, 0, 1.1);
  result_check.SetElement(0, 1, 1.2);
  result_check.SetElement(1, 0, 1.3);
  result_check.SetElement(1, 1, 1.4);

  S21Matrix matrix_before = matrix1;
  S21Matrix result = matrix1.Transpose();

  for (int i = 0; i < result.GetRows(); ++i) {
    for (int j = 0; j < result.GetCols(); ++j) {
      EXPECT_NEAR(result.GetElement(i, j), result_check.GetElement(i, j),
                  1e-6);  // 1e-6 с точностью до 6 знаков после запятой т.е
                          // погрешность не должна быть больше 0.000001
    }
  }
  EXPECT_TRUE(matrix1.EqMatrix(matrix_before));
  ;
}
TEST(S21Matrix, minorNotSquare) {
  S21Matrix matrix(3, 2);
  EXPECT_THROW(matrix.Minor(1, 1), std::logic_error);
}
TEST(S21Matrix, minorOutOfRange) {
  S21Matrix matrix(2, 2);
  EXPECT_THROW(matrix.Minor(-1, 1), std::out_of_range);
}
TEST(S21Matrix, minor) {
  S21Matrix matrix(3, 3);
  matrix.SetElement(0, 0, 1.0);
  matrix.SetElement(0, 1, 2.0);
  matrix.SetElement(0, 2, 3.0);
  matrix.SetElement(1, 0, 4.0);
  matrix.SetElement(1, 1, 5.0);
  matrix.SetElement(1, 2, 6.0);
  matrix.SetElement(2, 0, 7.0);
  matrix.SetElement(2, 1, 8.0);
  matrix.SetElement(2, 2, 9.0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      std::cout << "Element (" << i << ", " << j
                << ") = " << matrix.GetElement(i, j) << std::endl;
    }
  }
  EXPECT_EQ(matrix.Minor(0, 0), -3);
}

TEST(S21Matrix, det) {
  S21Matrix matrix2x2(2, 2);
  matrix2x2.SetElement(0, 0, 1.0);
  matrix2x2.SetElement(0, 1, 2.0);
  matrix2x2.SetElement(1, 0, 3.0);
  matrix2x2.SetElement(1, 1, 4.0);
  EXPECT_DOUBLE_EQ(matrix2x2.Determinant(), -2.0);

  S21Matrix matrix3x3(3, 3);
  matrix3x3.SetElement(0, 0, 1.0);
  matrix3x3.SetElement(0, 1, 2.0);
  matrix3x3.SetElement(0, 2, 3.0);
  matrix3x3.SetElement(1, 0, 4.0);
  matrix3x3.SetElement(1, 1, 5.0);
  matrix3x3.SetElement(1, 2, 6.0);
  matrix3x3.SetElement(2, 0, 7.0);
  matrix3x3.SetElement(2, 1, 8.0);
  matrix3x3.SetElement(2, 2, 9.0);
  EXPECT_DOUBLE_EQ(matrix3x3.Determinant(), 0.0);

  S21Matrix matrix3x4{3, 4};
  EXPECT_THROW(matrix3x4.Determinant(), std::logic_error);
}

TEST(S21Matrix, setElementExcept) {
  S21Matrix matrix(2, 2);
  EXPECT_THROW(matrix.SetElement(-1, 2, 3), std::out_of_range);
}
TEST(S21Matrix, getElementExcept) {
  S21Matrix matrix(2, 2);
  EXPECT_THROW(matrix.GetElement(-1, 2), std::out_of_range);
}

TEST(S21MatrixComplements, calcComp) {
  S21Matrix matrix(3, 3);
  matrix.SetElement(0, 0, 1);
  matrix.SetElement(0, 1, 2);
  matrix.SetElement(0, 2, 3);
  matrix.SetElement(1, 0, 4);
  matrix.SetElement(1, 1, 5);
  matrix.SetElement(1, 2, 6);
  matrix.SetElement(2, 0, 7);
  matrix.SetElement(2, 1, 8);
  matrix.SetElement(2, 2, 9);
  S21Matrix matrixCompl = matrix.CalcComplements();

  EXPECT_EQ(matrixCompl.GetElement(0, 0), -3);
  EXPECT_EQ(matrixCompl.GetElement(0, 1), 6);
  EXPECT_EQ(matrixCompl.GetElement(0, 2), -3);
  EXPECT_EQ(matrixCompl.GetElement(1, 0), 6);
  EXPECT_EQ(matrixCompl.GetElement(1, 1), -12);
  EXPECT_EQ(matrixCompl.GetElement(1, 2), 6);
  EXPECT_EQ(matrixCompl.GetElement(2, 0), -3);
  EXPECT_EQ(matrixCompl.GetElement(2, 1), 6);
  EXPECT_EQ(matrixCompl.GetElement(2, 2), -3);
}
TEST(S21MatrixComplements, calcCompExcept) {
  S21Matrix matrix(4, 3);
  EXPECT_THROW(S21Matrix matrixCompl = matrix.CalcComplements(),
               std::logic_error);
}

TEST(S21MatrixInverse, inverse) {
  S21Matrix matrix(2, 2);
  matrix.SetElement(0, 0, 2);
  matrix.SetElement(0, 1, 3);
  matrix.SetElement(1, 0, 1);
  matrix.SetElement(1, 1, 4);
  S21Matrix inverse = matrix.InverseMatrix();
  double inverseToCheck[4] = {0.8, -0.6, -0.2, 0.4};
  EXPECT_TRUE(matrixEqual(inverse, inverseToCheck));
}
TEST(S21MatrixInverse, inverseExcept) {
  S21Matrix matrix(3, 2);
  EXPECT_THROW(matrix.InverseMatrix(), std::logic_error);
}
TEST(S21MatrixInverse, inverseExcept2) {
  S21Matrix matrix(2, 2);
  matrix.SetElement(0, 0, 0);
  matrix.SetElement(0, 1, 0);
  matrix.SetElement(1, 0, 0);
  matrix.SetElement(1, 1, 0);
  EXPECT_THROW(matrix.InverseMatrix(), std::logic_error);
}