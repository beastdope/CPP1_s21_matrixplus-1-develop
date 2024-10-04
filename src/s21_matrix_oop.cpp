#include "s21_matrix_oop.h"

#include <iostream>

S21Matrix::S21Matrix() {
  rows_ = 4;
  cols_ = 4;
  matrix_ = new double[16]{};

  for (int i = 0; i < 4; ++i) {
    matrix_[i * cols_ + i] = 1;  // Устанавливаем единицы на главной диагонали
  }
}
/**
 * Конструктор S21Matrix инициализирует матрицу с указанными размерами,
 * устанавливает диагональные элементы равными 1 и выводит сообщение, если
 * создается новая матрица.
 *
 * @param rows Параметр rows в конструкторе S21Matrix представляет количество
 * строк в создаваемой матрице. Он определяет вертикальный размер матрицы.
 */
S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 0 || cols < 0) {
    throw std::invalid_argument("matix is not square");
  }
  if (rows == 0 || cols == 0) {
    throw std::invalid_argument("matix is 0");
  }

  rows_ = rows;
  cols_ = cols;

  matrix_ = new double[static_cast<long>(rows_) * cols_]{};

  int min_dim = std::min(rows_, cols_);
  for (int i = 0; i != min_dim; ++i) {
    matrix_[i * cols_ + i] = 1;
  }
  bool news_log = true;
  if (news_log) std::cerr << ": matrix is created" << std::endl;
}

/**
 * Конструктор S21Matrix инициализирует матрицу с предоставленными данными и
 * размерами и выполняет проверку ошибок на наличие недопустимых аргументов.
 *
 * @param matrix Параметр «matrix» — это указатель на постоянный массив двойных
 * значений.
 * @param rows Строки представляют количество строк в матрице. Он указывает
 * горизонтальный размер матрицы, определяя, сколько горизонтальных строк
 * элементов присутствует в матрице.
 * @param cols Параметр cols в конструкторе S21Matrix представляет количество
 * столбцов в матрице, которая передается в качестве входных данных. Он
 * определяет ширину матрицы.
 */
S21Matrix::S21Matrix(const double matrix[], int rows, int cols) {
  if (matrix == nullptr) {
    throw std::invalid_argument("null ptr");
  }
  if (rows < 0 || cols < 0) {
    throw std::invalid_argument("cols and rows are 0");
  }
  if (rows == 0 || cols == 0) {
    throw std::invalid_argument("rows and cols are 0");
  }
  rows_ = rows;
  cols_ = cols;
  matrix_ = new double[static_cast<long>(rows_) * cols_];
  std::copy_n(matrix, static_cast<long>(rows_) * cols_, matrix_);
  bool news_log = true;
  if (news_log) std::cerr << ": matris is created";
}

/**
 * Конструктор копирования S21Matrix инициализирует новую матрицу, копируя
 * данные из другой матрицы.
 *
 * @param other Параметр `other` в конструкторе копирования
 * `S21Matrix::S21Matrix` является ссылкой на другой объект `S21Matrix`, из
 * которого необходимо скопировать данные. Предоставленный вами фрагмент кода
 * показывает реализацию конструктора копирования, в котором элементы данных
 * (`rows_`, `cols
 */
S21Matrix::S21Matrix(const S21Matrix &other) {  // copy

  rows_ = other.rows_;
  cols_ = other.cols_;

  matrix_ = new double[static_cast<long>(rows_) * cols_]{};
  std::copy_n(other.matrix_, rows_ * cols_, matrix_);
  bool news_log = true;
  if (news_log) std::cerr << ": matrix is copied";
}

/**
 * Конструктор перемещения S21Matrix передает право владения данными от другого
 * объекта к текущему объекту.
 *
 * @param other В данном фрагменте кода «other» — это параметр типа «S21Matrix»,
 * из которого выполняется перемещение. Конструктор перемещения передает ресурсы
 * (например, владение памятью) из «другого» объекта вновь созданному объекту.
 */
S21Matrix::S21Matrix(S21Matrix &&other) {  // move
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;
  other.matrix_ = nullptr;
  bool news_log = true;
  if (news_log) std::cerr << ": matrix is moved)";
}

/**
 * Функция перегружает оператор «+=", чтобы добавить еще один S21Matrix к
 * текущему.
 *
 * @return Возвращается указатель `*this`.
 */
S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}
/**
 * Функция перегружает оператор вычитания-присваивания для класса S21Matrix,
 * вычитая другой объект S21Matrix из текущего объекта.
 *
 * @return Возвращается указатель `*this`.
 */
S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}
S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}
/**
 * Функция перегружает оператор умножения-присваивания для класса S21Matrix,
 * умножая другой объект S21Matrix из текущего объекта.
 *
 * @return Возвращается указатель `*this`.
 */
S21Matrix &S21Matrix::operator*=(double other) {
  MulNumber(other);
  return *this;
}

/**
 * Функция перегружает оператор прибавления для класса S21Matrix,
 * прибавляя другой объект S21Matrix из текущего объекта.
 *
 * @return Возвращается указатель `*this`.
 */
S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix sum = *this;
  sum += other;
  return sum;
}
/**
 * Функция перегружает оператор вычитания для класса S21Matrix,
 * вычитая другой объект S21Matrix из текущего объекта.
 *
 * @return Возвращается указатель `*this`.
 */
S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix sum = *this;
  sum -= other;
  return sum;
}
/**
 * Функция перегружает оператор умножения для класса S21Matrix,
 * умножая другой объект S21Matrix из текущего объекта.
 *
 * @return Возвращается указатель `*this`.
 */
S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix sum = *this;
  sum *= other;
  return sum;
}
S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (this != &other) {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(matrix_, other.matrix_);
  }

  return *this;
}
S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  S21Matrix copy{other};
  *this = std::move(copy);
  return *this;
}
bool S21Matrix::operator==(const S21Matrix &other) const {
  return EqMatrix(other);
}
bool S21Matrix::operator!=(const S21Matrix &other) const {
  return !EqMatrix(other);
}

double &S21Matrix::operator()(int row, int col) const {
  if (row < 0 || col < 0 || row >= rows_ || col >= cols_) {
    throw std::out_of_range("indexes out of range");
  }
  return matrix_[row * cols_ + col];
}

/**
 * Функция перегружает оператор умножения для класса S21Matrix,
 * умножая другой объект S21Matrix из текущего объекта.
 *
 * @return Возвращается указатель `*this`.
 */
S21Matrix S21Matrix::operator*(double other) {
  S21Matrix sum = *this;
  sum *= other;
  return sum;
}

double &S21Matrix::at(int row, int col) const {
  return (*this)(row, col);
}  // возвращает ссылку на элемент по заданному индексу

bool S21Matrix::EqMatrix(const S21Matrix &other) const {
  if (rows_ != other.rows_ || cols_ != other.cols_) return false;

  long n_elements = static_cast<long>(rows_) * cols_;
  for (long i = 0; i < n_elements; ++i) {
    if (matrix_[i] != other.matrix_[i]) {
      return false;
    }
  }
  return true;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("S21Matrix: wrong dimensions");
  }
  long index = static_cast<long>(rows_) * cols_;
  for (long i = 0; i < index; ++i) {
    matrix_[i] += other.matrix_[i];
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("S21Matrix: wrong dimensions");
  }
  long index = static_cast<long>(rows_) * cols_;
  for (int i = 0; i < index; ++i) {
    matrix_[i] -= other.matrix_[i];
  }
}

void S21Matrix::MulNumber(const double num) {
  long index = static_cast<long>(rows_) * cols_;
  for (int i = 0; i < index; ++i) {
    matrix_[i] *= num;
  }
}
/**
 * Функция MulMatrix умножает две матрицы и сохраняет результат в вызывающем
 * объекте.
 *
 * @param other Похоже, что предоставленный вами фрагмент кода представляет
 * собой функцию MulMatrix в классе S21Matrix, которая умножает две матрицы.
 * Функция принимает другой объект `S21Matrix` `other` в качестве параметра и
 * выполняет умножение матрицы на текущий матричный объект.
 */

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_) {
    throw std::invalid_argument("can't multiply matrix");
  }
  double *result = new double[static_cast<long>(rows_) * other.cols_]{};
  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < other.cols_; ++j) {
      for (int k = 0; k < cols_; ++k) {
        result[i * other.cols_ + j] +=
            matrix_[i * cols_ + k] * other.matrix_[k * other.cols_ + j];
      }
    }
  }

  delete[] matrix_;
  matrix_ = result;
  cols_ = other.cols_;
}

// из строчной матрицы в столбец переводит
S21Matrix S21Matrix::Transpose() const {
  S21Matrix transpose(cols_, rows_);

  for (int i = 0; i < rows_; ++i) {
    for (int j = 0; j < cols_; ++j) {
      transpose.matrix_[j * transpose.cols_ + i] = matrix_[i * cols_ + j];
    }
  }
  return transpose;
}

void S21Matrix::RemoveRowCol(const double *matrix, int matrix_size,
                             double *target, int row_skip, int col_skip) const {
  int i = 0;
  for (int row = 0; row < matrix_size; ++row) {
    if (row == row_skip) continue;
    for (int col = 0; col < matrix_size; ++col) {
      if (col != col_skip) {
        target[i] = matrix[row * matrix_size + col];
        i++;
      }
    }
  }
}
/** Методом разложения по первой строке
 * Если матрица 2х2 то вычисляется по формуле mat[0] * mat[3] - mat[1] * mat[2]
 * Если матрица 1х1 то её детерминант всегда равен первому элементу этой матрицы
 * т.е | n | - это и будет её детерминант
 * Минор матрицы вычисляется путем удаления определенной строки и столбца в
 * матрице затем нахождение его определителя т.е чтобы найти определитель
 * сначала надо найти минор элемента матрицы, если матрица больше 2х2 так как в
 * 2х2 определитель равен умножая (1 * 4) - (2 * 3) = 4 - 6 = определитель -2
 *                  |1, 2|
 *                  |3, 4|
 *
 * если матрица 3х3, то сначала мы вычисляем минор элемента например если
 * матрица  | 1, 2, 3 |
 *          | 4, 5, 6 |
 *          | 7, 8, 9 |
 *
 * нам необходимо найти минор элемента матрицы и затем по формуле нахождения
 * определителя разложим определитель по первой строке
 * det(A) = 1 * | 5, 6 | - 2 *| 4, 6 | + 3 * | 4, 5 |
 *              | 8, 9 |      | 7, 9 |       | 7, 8 |
 * затем находим определитель миноров
 * | 5, 6 | = (5 * 9) - (6 * 8) = -3
 * | 8, 9 |
 *
 * | 4, 6 | = (4 * 9) - (6 * 7) = -6
 * | 7, 9 |
 *
 * | 4, 5 | = (4 * 8) - (5 * 7) = -3
 * | 7, 8 |
 * затем всё подставляем и получается
 * det(A) = 1 * (-3) - 2 *(-6) + 3 * (-3) =
 * 1 * (-3) = -3
 * 2 * (-6) = -12
 * 3 * (-3) = -9
 * -3 + 12 - 9 =
 * det = 9 - 9 = 0
 * определитель матрицы 3 на 3 равен 0
 */
/**
 * Функция вычисляет определитель квадратной матрицы, используя рекурсивные
 * вызовы и второстепенные матрицы.
 *
 * @param mat Параметр mat в функции Det представляет собой указатель на
 * элементы матрицы, хранящиеся в порядке строк. Предполагается, что матрица
 * представляет собой квадратную матрицу размера «e x e», где «e» — количество
 * строк (или столбцов) в матрице. Элементы
 * @param e Параметр e в функции Det представляет размер квадратной матрицы, для
 * которой вы хотите вычислить определитель. Он используется для определения
 * размера матрицы и выполнения необходимых вычислений для нахождения
 * определителя.
 *
 * @return Функция Determinant() в классе S21Matrix возвращает определитель
 * матрицы, хранящийся в переменной-члене matrix_. Сначала он проверяет,
 * является ли матрица квадратной (т. е. количество строк равно количеству
 * столбцов), а затем вызывает функцию Det(), передавая матрицу и ее размер для
 * вычисления определителя. Затем возвращается вычисленное значение
 * определителя.
 */
double S21Matrix::Det(const double *mat, int e) const {
  double determinant = 0;
  if (e == 1) return mat[0];
  if (e == 2) return mat[0] * mat[3] - mat[1] * mat[2];

  double *tmp = new double[static_cast<unsigned long>(e - 1) * (e - 1)]{};
  int sign = 1;
  for (int i = 0; i < e; ++i) {
    RemoveRowCol(mat, e, tmp, 0, i);
    std::cout << "Intermediate determinant for element (" << 0 << ", " << i
              << ") is: " << sign * mat[i] * Det(tmp, e - 1) << std::endl;

    determinant += sign * mat[i] * Det(tmp, e - 1);
    sign = -sign;
  }
  delete[] tmp;
  return determinant;
}

double S21Matrix::Determinant() const {
  if (rows_ != cols_) {
    throw std::logic_error("S21Matrix is not a square");
  }
  return Det(matrix_, rows_);
}

double S21Matrix::Minor(int row, int col) const {
  if (row < 0 || col < 0 || row > rows_ - 1 || col > cols_ - 1) {
    throw std::out_of_range("S21Matrix: is out of range");
  }
  if (rows_ != cols_) {
    throw std::logic_error("S21Matrix: is not a square");
  }
  int e = rows_;

  double *tmp = new double[static_cast<unsigned long>(e - 1) * (e - 1)]{};
  RemoveRowCol(matrix_, e, tmp, row, col);

  double minor = Det(tmp, e - 1);
  delete[] tmp;
  return minor;
}
/*
 * Вычисление алгебраических дополнение это число которое используется в
 * операциях с матрицами, в нашем случае это необходимо для нахождения обратной
 * матрицы, алгебраическое дополнение элемента (A)ij матрицы А определяется так
 * Aij = (-1)i+j * Mij
 * M - минор элемента aij
 * Например если матрица A = | 1, 2, 3 |
 * Например если матрица A = | 4, 5, 6 |
 * Например если матрица A = | 7, 8, 9 |
 * То алгебраические дополнения будут вычисляться следующим образом:
 * A = Ai+j * Mij
 * A00 = M00 = | 5, 6 |
 * A00 = M00 = | 8, 9 |
 * находим определитель минора = (5*9) - (6 * 8) = 45 - 48 = -3
 * Затем находим Алгебраическое дополнение элемента
 * A00 = (-1)0 + 0 * (-3) = 1 *-3 = 3
 *
 * A00 = -3 // 00
 *
 *
 * M01 = | 4, 6 |
 * M01 = | 7, 9 |
 * det01 = (4 * 9) - (6 * 7) = 36 - 42 = -6
 * A01 = (-1)0+1 * (-6) = -1 * -6 = 6
 *
 * A01 = 6 // 01
 *
 * M02 = | 4, 5 |
 * M02 = | 7, 8 |
 * det02 = (4 * 8) - (5 * 7) = 32 - 35 = -3
 * A02 = (-1)0+2 * (-3) = 1 * -3 = -3
 *
 * A02 = -3 // 02

 *
 * M10 = | 2, 3 |
 * M10 = | 8, 9 |
 * det10 = 18 - 24 = -6
 *
 * A10 = (-1)1+0 * (-6) = -1 * -6 = 6
 *
 * A10 = 6 // 10
 *
 * M11 = | 1, 3 |
 * M11 = | 7, 9 |
 * det11 = 9 - 21 = -12
 * A11 = (-1)1+1 * (-12) = 1 * -12 = -12
 *
 * A11 = -12 // 11
 *
 *
 * M12 = | 1, 2 |
 * M12 = | 7, 8 |
 * det12 = 8 - 14 = -6
 * Alg = (-1)1+2 *(-6) = -1 * -6 = 6
 *
 * A12 == 6 // 12
 *
 * M20 = | 2, 3 |
 * M20 = | 5, 6 |
 * det20 = 12 - 15 = -3
 * Alg = (-1)2+0 * (-3) = 1 * -3 = -3
 *
 * A20 = -3 // 20
 *
 *
 * A21 = | 1, 3 |
 * A21 = | 4, 6 |
 * M21 = 6 - 12 = -6
 * Alg = (-1)2+1 * (-6) = -1 * -6 = 6
 *
 * A21 = 6 // 21
 *
 * Alg = | -3, 6, -3 |
 * Alg = | 6,-12, 6 |
 * Alg = | -3, 6,-3 |
 *
 *
 * A22
 * M22 = | 1, 2 |
 * M22 = | 4, 5 |
 * det22 = 5 - 8 = -3
 * A22 = (-1)2+2 * (-3) = 1 * -3 = -3
 *
 */
S21Matrix S21Matrix::CalcComplements() const {
  if (rows_ != cols_) {
    throw std::logic_error("S21Matrix: not a square");
  }
  int e = rows_;
  double *complements = new double[static_cast<unsigned long>(e) * e];
  for (int row = 0; row < e; ++row) {
    for (int col = 0; col < e; ++col) {
      int sign = (row + col) % 2 == 0 ? 1 : -1;
      complements[row * e + col] = sign * Minor(row, col);
    }
  }
  S21Matrix result(complements, e, e);
  delete[] complements;
  return result;
}
/**
 *  Обратная матрица работает по схеме следующей:
 * Сначала проверяется крадратная ли матрица если нет - ошибка
 * затем вычисляется детерминант матрицы, если близок к нулю = ошибка
 * если все условия соблюдены сначала вычисляем алгебраические дополнения, затем
 * транспонируем затем находим обратное значение определителя матрицы и затем
 * умножаем его на каждый элемент матрицы
 * A = | 2, 3 |
 * A = | 1, 4 |
 * det(A) = (2 * 4) - ( 3 * 1) = 5
 * После находим алгебраические дополнения(кофакторы)
 * minor(A11) = 4
 * compl = (-1)1+1 * 4 = 1 * 4 = 4
 * minor(A12) = 1
 * compl = (-1)1+2 * 1 = -1 * 1 = -1
 * minor(A21) = 3
 * compl = (-1)2+1 * 3 = -1 * 3 = -3
 * minor(A22) = 2
 * compl = (-1)2+2 * 2 = 1 * 2 = 2
 * Затем транспонируем полученную матрицу
 * A before transpose = | 4,-1 |
 * A before transpose = |-3, 2 |
 *
 * A after transpose =  | 4,-3 |
 * A after transpose =  |-1, 2 |
 * Затем находим обратный определитель путем деления на 1 а затем умножаем
 * каждый элемент на это значение Обратный Определитель = 1 / 5 = 0.2
 * теперь умножаем на каждый элемент
 * A11 = 0.2 * 4 = 0.8
 * A12 = 0.2 * (-3) = -0.6
 * A21 = 0.2 * (-1) = -0.2
 * A22 = 0.2 * (2) = 0.4
 * Итого:
 * Обратная матрица А = |0.8, -0.6 |
 * Обратная матрица А = |-0.2, 0.4 |
 */
S21Matrix S21Matrix::InverseMatrix() const {
  if (rows_ != cols_) {
    throw std::logic_error("Is not square");
  }
  double determinant = Determinant();
  if (std::fabs(determinant) < 1e-7) {  // abs (абсолютное число != -)
    throw std::logic_error("singular matrix");
  }
  S21Matrix inverse = CalcComplements().Transpose();
  inverse.MulNumber(
      1 /
      determinant);  // находим обратное значение определителя путем умножения
                     // обратного определителя на каждый элемент матрицы
  return inverse;
}

int S21Matrix::GetRows() const { return rows_; }
int S21Matrix::GetCols() const { return cols_; }
void S21Matrix::SetRows(int new_rows) {
  if (new_rows < 1) {
    throw std::out_of_range("S21Matrix: rows are < 1");
  }
  double *new_matrix = new double[static_cast<long>(new_rows) * cols_]{};
  std::copy_n(matrix_, std::min(new_rows, rows_) * cols_, new_matrix);
  rows_ = new_rows;
  delete[] matrix_;
  matrix_ = new_matrix;
}

void S21Matrix::SetCols(int new_col) {
  if (new_col < 1) {
    throw std::out_of_range("Invalid col");
  }

  double *new_matrix = new double[static_cast<long>(rows_) * new_col]{};
  int col_copy = std::min(new_col, cols_);
  for (int i = 0; i != rows_; ++i) {
    for (int j = 0; j != col_copy; ++j) {
      new_matrix[i * new_col + j] = matrix_[i * cols_ + j];
    }
  }

  cols_ = new_col;
  delete[] matrix_;
  matrix_ = new_matrix;
}

double *S21Matrix::GetMatrix() const { return matrix_; }

void S21Matrix::SetElement(int row, int col, double value) {
  if (rows_ < 0 || col < 0 || row >= rows_ || col >= cols_) {
    throw std::out_of_range("Indexes are out of range");
  }
  matrix_[row * cols_ + col] = value;
}

float S21Matrix::GetElement(int row, int col) const {
  if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
    // Обработка ошибки: выбросить исключение или вернуть значение по умолчанию
    throw std::out_of_range("Индекс строки или столбца вне диапазона");
  }
  return matrix_[row * cols_ + col];
}
/**
 * Деструктор класса S21Matrix удаляет матрицу и записывает сообщение, если
 * значение news_log истинно.
 */
S21Matrix::~S21Matrix() {
  if (matrix_ != nullptr) {
    delete[] matrix_;
    matrix_ = nullptr;
  }
  bool news_log = true;
  if (news_log) std::cerr << (": matrix is deleted");
}
