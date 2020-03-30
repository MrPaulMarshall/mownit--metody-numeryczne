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
        return NULL;
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
        return NULL;
    }

    AGHMatrix newMatrix = new AGHMatrix(this->get_rows(), rhs.get_cols(), 0);
    for (int i = 0; i < newMatrix.get_rows(); i++) {
        for (int j = 0; j < newMatrix.get_cols(); j++) {
            for (int k = 0; k < rhs.get_rows(); k++) {
                newMatrix(i, j) += this(i, k) * rhs(k, j);
            }
        }
    }

    return newMatrix;
}

// Printing matrix                                                                                                                        
template<typename T>
std::ostream& operator<<(std::ostream& stream, const AGHMatrix<T>& matrix) 
{
  for (int i=0; i<matrix.get_rows(); i++) 
  { 
    for (int j=0; j<matrix.get_cols(); j++) 
    {
        stream << matrix(i,j) << ", ";
    }
    stream << std::endl;
  }
    stream << std::endl;
}