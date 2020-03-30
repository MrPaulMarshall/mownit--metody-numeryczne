#include "aghMatrix.h"


// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(const std::vector<std::vector<T>>& mat) 
{
  matrix.resize(mat.size());
  for (unsigned i = 0; i < mat.size(); i++) 
  {
    matrix[i].resize(mat[i].size());
    for(unsigned j = 0; j < mat[i].size(); j++)
    {
      matrix[i][j] = mat[i][j];
    }
  }
  rows = matrix.size();
  cols = matrix[1].size();
}

// Parameter Constructor                                                                                                                                                      
template<typename T>
AGHMatrix<T>::AGHMatrix(unsigned _rows, unsigned _cols, const T& _initial) 
{
  matrix.resize(_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(_cols, _initial);
  }
  rows = _rows;
  cols = _cols;
}

// Copy Constructor                                                                                                                                                           
template<typename T>
AGHMatrix<T>::AGHMatrix(const AGHMatrix<T>& rhs) 
{
  matrix = rhs.matrix;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// Get the number of rows of the matrix                                                                                                                                       
template<typename T>
unsigned AGHMatrix<T>::get_rows() const 
{
  return this->rows;
}

// Get the number of columns of the matrix                                                                                                                                    
template<typename T>
unsigned AGHMatrix<T>::get_cols() const 
{
  return this->cols;
}

// Assignment Operator                                                                                                                                                        
template<typename T>
AGHMatrix<T>& AGHMatrix<T>::operator=(const AGHMatrix<T>& rhs) 
{
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();

  matrix.resize(new_rows);
  for (unsigned i=0; i<matrix.size(); i++) 
  {
    matrix[i].resize(new_cols);
  }

  for (unsigned i=0; i<new_rows; i++) 
  {
    for (unsigned j=0; j<new_cols; j++) 
    {
      matrix[i][j] = rhs(i, j);
    }
  }
  rows = new_rows;
  cols = new_cols;

  return *this;
}

// compares two matrices
template<typename T>
bool AGHMatrix<T>::operator==(const AGHMatrix<T>& rhs)
{
    if (&rhs == this)
        return true;

    if (this->get_rows() != rhs.get_rows() || this->get_cols() != rhs.get_cols())
        return false;

    for (int i = 0; i < rhs.get_rows(); i++) {
        for (int j = 0; j < rhs.get_cols(); j++) {
            if (abs(this->matrix[i][j] - rhs(i, j)) > cEPS)
                return false;
        }
    }
    
    return true;
}

// Access the individual elements                                                                                                                                             
template<typename T>
T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) 
{
  return this->matrix[row][col];
}

// Access the individual elements (const)                                                                                                                                     
template<typename T>
const T& AGHMatrix<T>::operator()(const unsigned& row, const unsigned& col) const 
{
  return this->matrix[row][col];
}

// Addition of two matrices                                                                                                                                                   
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator+(const AGHMatrix<T>& rhs) 
{
  // Task 1 - implement addition of two matrices
    if (this->get_rows() != rhs.get_rows() || this->get_cols() != rhs.get_cols()) {
        std::cout << "Error: attempt to add matrices with different sizes\n";
        exit(1);
    }
    
    AGHMatrix newMatrix = *this;
    for (int i = 0; i < newMatrix.get_rows(); i++) {
        for (int j = 0; j < newMatrix.get_cols(); j++) {
            newMatrix(i, j) += rhs(i, j);
        }
    }

    return newMatrix;
}

// Left multiplication of this matrix and another                                                                                                                              
template<typename T>
AGHMatrix<T> AGHMatrix<T>::operator*(const AGHMatrix<T>& rhs) 
{
  // Task 1 - implement multiplication of two matrices
    if (this->get_cols() != rhs.get_rows()) {
        std::cout << "Error: attempt to multiplicate matrices when number of cols in A is not equal to number of rows in B\n";
        exit(1);
    }

    AGHMatrix newMatrix = AGHMatrix(this->get_rows(), rhs.get_cols(), 0);
    for (int i = 0; i < newMatrix.get_rows(); i++) {
        for (int j = 0; j < newMatrix.get_cols(); j++) {
            for (int k = 0; k < rhs.get_rows(); k++) {
                newMatrix(i, j) += this->matrix[i][k] * rhs(k, j);
            }
        }
    }

    return newMatrix;
}

// Check if matrix is symmetric
template<typename T>
bool AGHMatrix<T>::isSymmetric()
{
    if (this->get_rows() != this->get_cols()) {
        return false;
    }
    for (int i = 0; i < this->get_rows(); i++) {
        for (int j = i; j < this->get_cols(); j++) {
            if ((this->matrix[i][j] - this->matrix[j][i]) > cEPS)
                return false;
        }
    }

    return true;
}

// Calculate daterminant of a matrix
template<typename T>
T AGHMatrix<T>::det()
{
    if (this->get_rows() != this->get_cols())
        return 0;

    unsigned n = this->get_rows();
    // if matrix is 1x1, then determinant is just that element
    if (n == 1)
        return this->matrix[0][0];
    
    T det = 0;
    int sign;
    for (int k = 0; k < n; k++) {
        // i use recursive definition of determinant
        //       so i need the smaller matrix in each step
        AGHMatrix<T> innerMatrix(n - 1, n - 1, 0);
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < n - 1; j++)
                innerMatrix(i, j) = this->matrix[i][j];
        }
        for (int i = k + 1; i < n; i++) {
            for (int j = 0; j < n - 1; j++)
                innerMatrix(i, j) = this->matrix[i][j];
        }
        if ((n - 1 + k) % 2 == 0)
            sign = 1;
        else
            sign = -1;
        det += sign * innerMatrix.det() * this->matrix[k][n - 1];
    }
    return det;
}

// Transpose matrix
template<typename T>
AGHMatrix<T> AGHMatrix<T>::transpose()
{
    AGHMatrix<T> transposedMatrix(this->get_cols(), this->get_rows(), 0);
    for (int i = 0; i < this->get_rows(); i++) {
        for (int j = 0; j < this->get_cols(); j++) {
            transposedMatrix(j, i) = this->matrix[i][j];
        }
    }

    return transposedMatrix;
}

// Printing matrix                                                                                                                        
template<typename T>
std::ostream& operator<<(std::ostream& stream, const AGHMatrix<T>& matrix) 
{
  for (int i = 0; i < matrix.get_rows(); i++) 
  { 
    for (int j = 0; j < matrix.get_cols(); j++) 
    {
        stream << matrix(i,j) << ", ";
    }
    stream << std::endl;
  }
  
  stream << std::endl;

  return stream;
}
