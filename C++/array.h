// include the global includes
#include "includes.h"

// make sure only included once
#ifndef array_H
#define array_H

/**
 * @brief Describes a Matrix which are used storing matrices more easily
 * 
 * @tparam T is the container type
 * @tparam U is the value type
 */
template <template <typename...> class T, typename U>
class Matrix {
    public:
        /**
         * @brief Constructor for a Matrix class
         */
        Matrix() {
            size_1 = 0; size_2 = 0;
        };
        ~Matrix() {};

    protected:
        // define the sizes of the matrix
        int size_1, size_2;

    public:
        /**
         * @brief Get the size of the matrix
         * 
         * @return the size of the matrix
         */
        tuple<int, int> get_size() const { return make_tuple(size_1, size_2); };

        /**
         * @brief Get the element of the matrix at (i, j)
         * 
         * @param i index
         * @param j index
         * @return the element at (i, j)
         */
        virtual U get_element(int i, int j) const { return 0; };
    
        /**
         * @brief stream operator for the Matrix to print it out more easily
         * 
         * @param os output stream to parse the object to 
         * @param matrix object in question
         * @return the output stream
         */
        friend ostream& operator<<(ostream& os, const Matrix& matrix) {

            // get the size of the matrix
            auto [size_1, size_2] = matrix.get_size();

            // loop over it's elements and print it out
            for (int i = 0; i < size_1; i++) {
                for (int j = 0; j < size_2; j++) {
                    os << matrix.get_element(i, j) << " ";
                }

                // use a linebreak after each row
                os << endl;
            }

            return os;
        };
};

/**
 * @brief Describes a nested array matrix that is based on the base matrix
 * 
 * @tparam T is the container type
 * @tparam U is the value type
 */
template <template <typename...> class T, typename U>
class Nested_Array : public Matrix<T, U> {

    public:
        /**
         * @brief Constructor for a Nested_Array class, creating a matrix from
         *        multiplying the values of 2 arrays.
         * 
         * @tparam V is the container type
         * @param arr1 first array;
         * @param arr2 second arrray;
         */
        template <template <typename...> class V>
        Nested_Array(V<U> arr1, V<U> arr2) : Matrix<T, U>() {

            // get the sizes of the arrays and set them as the matrix sizes
            this->size_1 = arr1.size(); this->size_2 = arr2.size();

            // define the iterators of the arrays
            typename V<U>::iterator a1 = arr1.begin(), a2 = arr2.begin();

            // loop over array 1
            for (int i = 0; i < this->size_1; i++) {

                // define a row array
                T<U> row;

                // loop over array 2
                for (int j = 0; j < this->size_2; j++) {

                    // add the multiple of the two elements
                    row.push_back((*next(a1, i)) * (*next(a2, j)));
                }

                // add the row
                _matrix.push_back(row);
            }
        };

        /**
         * @brief Copy Constructor of a Nested_Array class, copying a matrix
         * 
         * @param _new Nested_Array to copy;
         */
        Nested_Array(const Nested_Array<T, U>& _new) : Matrix<T, U>() {

            // get the sizes of the matrix and set them for this one
            tie(this->size_1, this->size_2) = _new.get_size();

            // loop over axis 1
            for (int i = 0; i < this->size_1; i++) {

                // make a row array
                T<U> row;

                // loop over axis 2 and add the element to the row
                for (int j = 0; j < this->size_2; j++) {
                    row.push_back(_new.get_element(i, j));
                }

                // add the row to the matrix
                _matrix.push_back(row);
            }
        }

        /**
         * @brief Constructor of a Nested_Array class from sizes
         * 
         * @param n1 size 1;
         * @param n2 size 2;
         */
        Nested_Array(int n1, int n2) : Matrix<T, U>() {

            // set the sizes of the matrix
            this->size_1 = n1; this->size_2 = n2;

            // loop over axis 1
            for (int i = 0; i < this->size_1; i++) {

                // create a row array
                T<U> row;

                // loop over axis 2 and initialize to 0
                for (int j = 0; j < this->size_2; j++) {
                    row.push_back(0);
                }

                // add the row to the matrix
                _matrix.push_back(row);
            }
        }
        ~Nested_Array() {};

    private:
        // define the internal matrix
        T<T<U>> _matrix;

    public:
        /**
         * @brief Get the element of the matrix at (i, j)
         * 
         * @param i index
         * @param j index
         * @return the element at (i, j)
         */
        U& get_element(int i, int j) {

            // make sure the indices are within the sizes of the matrix
            assert(0 <= i && i < this->size_1 && 0 <= j && j < this->size_2);

            // define the iterators
            typename T<T<U>>::iterator row = _matrix.begin();
            typename T<U>::iterator element = (*next(row, i)).begin();

            // return the element
            return *next(element, j);
        }

        /**
         * @brief Get the element of the matrix at (i, j) (constant)
         * 
         * @param i index
         * @param j index
         * @return the element at (i, j)
         */
        U get_element(int i, int j) const {
            
            // make sure the indices are within the sizes of the matrix
            assert(0 <= i && i < this->size_1 && 0 <= j && j < this->size_2);

            // define the constant iterators
            typename T<T<U>>::const_iterator row = _matrix.begin();
            typename T<U>::const_iterator element = (*next(row, i)).begin();

            // return the element
            return *next(element, j);
        }

        /**
         * @brief Set the element of the matrix at (i, j)
         * 
         * @param i index
         * @param j index
         * @param val the value to set the element to
         */
        void set_element(int i, int j, U val) { get_element(i, j) = val; };

        /**
         * @brief Get the internal matrix (constant)
         * 
         * @return return the internal matrix
         */
        T<T<U>> get_matrix() const { return _matrix; };

        /**
         * @brief addition operator of the matrix
         * 
         * @param matrix1 matrix to add to
         * @param matrix2 matrix to be added 
         * @return result of the summation
         */
        friend Nested_Array operator+(Nested_Array& matrix1,
                Nested_Array& matrix2) {
            // get the sizes of the matrices
            auto [size_1_1, size_2_1] = matrix1.get_size();
            auto [size_1_2, size_2_2] = matrix2.get_size();

            // make sure the matrices are the same size
            assert(size_1_1 == size_1_2 && size_2_1 == size_2_2);

            // copy matrix one to a new matrix
            Nested_Array matrix_out(matrix1);

            // get the sizes of the new matrix
            auto [size_1, size_2] = matrix_out.get_size();

            // loop over each element and add the correct element of the other
            for (int i = 0; i < size_1; i++) {
                for (int j = 0; j < size_2; j++) {
                    matrix_out.get_element(i, j) += matrix2.get_element(i, j);
                }
            }

            // return the result
            return matrix_out;
        };
};

/**
 * @brief Describes a flat array matrix that is based on the base matrix
 * 
 * @tparam T is the container type
 * @tparam U is the value type
 */
template <template <typename...> class T, typename U>
class Array : public Matrix<T, U> {

    public:
        /**
         * @brief Constructor for a Array class, creating a matrix from
         *        multiplying the values of 2 arrays.
         * 
         * @tparam V is the container type
         * @param arr1 first array;
         * @param arr2 second arrray;
         */
        template <template <typename...> class V>
        Array(V<U> arr1, V<U> arr2) : Matrix<T, U>() {

            // get the sizes of the arrays and set them as the matrix sizes
            this->size_1 = arr1.size(); this->size_2 = arr2.size();

            // define the iterators of the arrays
            typename V<U>::iterator a1 = arr1.begin(), a2 = arr2.begin();

            // loop over each element and add the mutliple of them
            for (int i = 0; i < this->size_1; i++) {
                for (int j = 0; j < this->size_2; j++) {
                    _matrix.push_back((*next(a1, i)) * (*next(a2, j)));
                }
            }
        };

        /**
         * @brief Copy Constructor of a Array class, copying a matrix
         * 
         * @param _new Nested_Array to copy;
         */
        Array(const Array<T, U>& _new) : Matrix<T, U>() {

            // get the sizes of the matrix and set them for this one
            tie(this->size_1, this->size_2) = _new.get_size();

            // loop over elements and set the element to the corrisponding one
            for (int i = 0; i < this->size_1; i++) {
                for (int j = 0; j < this->size_2; j++) {
                    _matrix.push_back(_new.get_element(i, j));
                }
            }
        };

        /**
         * @brief Constructor of a Array class from sizes
         * 
         * @param n1 size 1;
         * @param n2 size 2;
         */
        Array(int n1, int n2) : Matrix<T, U>() {

            // set the sizes of the matrix
            this->size_1 = n1; this->size_2 = n2;

            // initialize all the elements to 0
            for (int i = 0; i < this->size_1; i++) {
                for (int j = 0; j < this->size_2; j++) {
                    _matrix.push_back(0);
                }
            }
        };
        ~Array() {};

    private:
        // define the internal matrix
        T<U> _matrix;

    public:
        /**
         * @brief Get the element of the matrix at (i, j)
         * 
         * @param i index
         * @param j index
         * @return the element at (i, j)
         */
        U& get_element(int i, int j) {

            // make sure the indices are within the sizes of the matrix
            assert(0 <= i && i < this->size_1 && 0 <= j && j < this->size_2);

            // define the iterator
            typename T<U>::iterator element = _matrix.begin();

            // return the element
            return *next(element, j + i * this->size_2);
        };

        /**
         * @brief Get the element of the matrix at (i, j) (constant)
         * 
         * @param i index
         * @param j index
         * @return the element at (i, j)
         */
        U get_element(int i, int j) const {

            // mae sure the indices are within the sizes of the matrix
            assert(0 <= i && i < this->size_1 && 0 <= j && j < this->size_2);

            // define the constant iterator
            typename T<U>::const_iterator element = _matrix.begin();

            // return the element
            return *next(element, j + i * this->size_2);
        };

        /**
         * @brief Set the element of the matrix at (i, j)
         * 
         * @param i index
         * @param j index
         * @param val the value to set the element to
         */
        void set_element(int i, int j, U val) { get_element(i, j) = val; };

        /**
         * @brief Get the pointer of the element of the matrix at (i, j)
         * 
         * @param i index
         * @param j index
         * @return pointer of the element at (i, j)
         */
        U* get_element_ptr(int i, int j) { return &get_element(i, j); };

        /**
         * @brief Get the pointer of the element of the matrix at (i, j)
         *        (constant)
         * 
         * @param i index
         * @param j index
         * @return pointer of the element at (i, j)
         */
        U* get_element_ptr(int i, int j) const { return &get_element(i, j); };

        /**
         * @brief Get the internal matrix
         * 
         * @return return the internal matrix
         */
        T<U>& get_matrix() { return _matrix; };

        /**
         * @brief Get the matrix as a nested array
         * 
         * @return the matrix as a nested array
         */
        T<T<U>> get_matrix_unfold() const {
            // get the sizes
            auto [size1, size2] = this->get_size();

            // define the output array
            T<T<U>> out;

            // loop over axis 1
            for (int i = 0; i < size1; i++) {

                // define a row
                T<U> out_;

                // loop over axis 2
                for (int j = 0; j < size2; j++) {
                    
                    // get the element
                    auto val = this->get_element(i, j);

                    // add the element to the row
                    out_.push_back(val);
                }

                out.push_back(out_);
            }

            // return the output array
            return out;
        }

        /**
         * @brief load a raw array into the matrix
         * 
         * @param arr array to be loaded in
         */
        void load_raw(T<U> arr) {
            // get the sizes of the matrix
            auto [size1, size2] = this->get_size();
            
            // make sure the size of the array and sizes of the matrix can work
            assert (size1 * size2 == arr.size());

            // set the internal matrix to the array to be loaded
            _matrix = arr;
        }

        /**
         * @brief load a raw nested_array into the matrix
         * 
         * @param arr array to be loaded in
         */
        void load_raw(T<T<U>> arr) {
            //get the sizes of the matrix
            auto [size1, size2] = this->get_size();

            // make sure the sizes of the array and sizes of the matirx can work
            assert (size1 == arr.size() && size2 == arr[0].size());
            
            // set the index to 0
            int index = 0;
            
            // loop over each row
            for (auto row : arr) {
                // copy the row to the matrix at the right index
                copy(row.begin(), row.end(), _matrix.begin() + index);

                // increase the index with the size of the row
                index += row.size();
            }
        }

        /**
         * @brief addition operator of the matrix
         * 
         * @param matrix1 matrix to add to
         * @param matrix2 matrix to be added
         * @return result of the summation
         */
        friend Array operator+(Array& matrix1, Array& matrix2) {
            // get the sizes of the matrices
            auto [size_1_1, size_2_1] = matrix1.get_size();
            auto [size_1_2, size_2_2] = matrix2.get_size();

            // make sure the matrices are the same size
            assert(size_1_1 == size_1_2 && size_2_1 == size_2_2);

            // copy matrix one to a new matrix
            Array matrix_out(matrix1);

            // get the sizes of the new matrix
            auto [size_1, size_2] = matrix_out.get_size();

            // loop over each element and add the correct element of the other
            for (int i = 0; i < size_1; i++) {
                for (int j = 0; j < size_2; j++) {
                    matrix_out.get_element(i, j) += matrix2.get_element(i, j);
                }
            }

            // return the result
            return matrix_out;
        };
};

#endif
