#include "includes.h"

#ifndef array_H
#define array_H

template <class T, class U>
class Matrix {
    public:
        Matrix() {
            _last_i = 0; _last_j = 0;
            size_1 = 0; size_2 = 0;
        };
        ~Matrix() {};

    protected:
        int _last_i, _last_j;
        int size_1, size_2;

    public:
        tuple<int, int> get_size() const { return make_tuple(size_1, size_2); };
        virtual U get_element(int i, int j) const { return 0; };
    
        friend ostream& operator<<(ostream& os, const Matrix& matrix) {
            auto [size_1, size_2] = matrix.get_size();

            for (int i = 0; i < size_1; i++) {
                for (int j = 0; j < size_2; j++) {
                    os << matrix.get_element(i, j) << " ";
                }

                os << endl;
            }

            return os;
        };
};

template <class T_T, class T, class U>
class Nested_Array : public Matrix<T, U> {

    public:
        Nested_Array(T arr1, T arr2) : Matrix<T, U>() {
            this->size_1 = arr1.size(); this->size_2 = arr2.size();

            typename T::iterator a1 = arr1.begin(), a2 = arr2.begin();

            for (int i = 0; i < this->size_1; i++) {
                T row;
                for (int j = 0; j < this->size_2; j++) {
                    row.push_back((*next(a1, i)) * (*next(a2, j)));
                }

                _matrix.push_back(row);
            }
        };
        Nested_Array(const Nested_Array<T, T_T, U>& _new) : Matrix<T, U>() {
            tie(this->size_1, this->size_2) = _new.get_size();

            for (int i = 0; i < this->size_1; i++) {
                T row;
                for (int j = 0; j < this->size_2; j++) {
                    row.push_back(_new.get_element(i, j));
                }

                _matrix.push_back(row);
            }
        }
        Nested_Array(int n1, int n2) : Matrix<T, U>() {
            this->size_1 = n1; this->size_2 = n2;

            for (int i = 0; i < this->size_1; i++) {
                T row;
                for (int j = 0; j < this->size_2; j++) {
                    row.push_back(0);
                }

                _matrix.push_back(row);
            }
        }
        ~Nested_Array() {};

    private:
        T_T _matrix;

    public:
        U& get_element(int i, int j) {
            assert(0 <= i && i < this->size_1 && 0 <= j && j < this->size_2);

            this->_last_i = i; this->_last_j;

            typename T_T::iterator row = _matrix.begin();
            typename T::iterator element = (*next(row, i)).begin();

            return *next(element, j);
        }

        U get_element(int i, int j) const {
            assert(0 <= i && i < this->size_1 && 0 <= j && j < this->size_2);

            typename T_T::const_iterator row = _matrix.begin();
            typename T::const_iterator element = (*next(row, i)).begin();

            return *next(element, j);
        }

        T_T get_matrix() const { return _matrix; };

        friend Nested_Array operator+(Nested_Array& matrix1, Nested_Array& matrix2) {
            auto [size_1_1, size_2_1] = matrix1.get_size();
            auto [size_1_2, size_2_2] = matrix2.get_size();

            assert(size_1_1 == size_1_2 && size_2_1 == size_2_2);

            Nested_Array matrix_out(matrix1);

            auto [size_1, size_2] = matrix_out.get_size();

            for (int i = 0; i < size_1; i++) {
                for (int j = 0; j < size_2; j++) {
                    matrix_out.get_element(i, j) += matrix2.get_element(i, j);
                }
            }

            return matrix_out;
        };
};

template <class T, class U>
class Array : public Matrix<T, U> {

    public:
        Array(T arr1, T arr2) : Matrix<T, U>() {
            this->size_1 = arr1.size(); this->size_2 = arr2.size();

            typename T::iterator a1 = arr1.begin(), a2 = arr2.begin();

            for (int i = 0; i < this->size_1; i++) {
                for (int j = 0; j < this->size_2; j++) {
                    _matrix.push_back((*next(a1, i)) * (*next(a2, j)));
                }
            }
        };
        Array(const Array<T, U>& _new) : Matrix<T, U>() {
            tie(this->size_1, this->size_2) = _new.get_size();

            for (int i = 0; i < this->size_1; i++) {
                for (int j = 0; j < this->size_2; j++) {
                    _matrix.push_back(_new.get_element(i, j));
                }
            }
        };
        Array(int n1, int n2) : Matrix<T, U>() {
            this->size_1 = n1; this->size_2 = n2;

            for (int i = 0; i < this->size_1; i++) {
                for (int j = 0; j < this->size_2; j++) {
                    _matrix.push_back(0);
                }
            }
        }
        ~Array() {};

    private:
        T _matrix;

    public:
        U& get_element(int i, int j) {
            assert(0 <= i && i < this->size_1 && 0 <= j && j < this->size_2);

            this->_last_i = i; this->_last_j;

            typename T::iterator element = _matrix.begin();

            return *next(element, j + i * this->size_2);
        };

        U get_element(int i, int j) const {
            assert(0 <= i && i < this->size_1 && 0 <= j && j < this->size_2);

            typename T::const_iterator element = _matrix.begin();

            return *next(element, j + i * this->size_2);
        };

        list<T> get_matrix() const {
            auto [size1, size2] = this->get_size();

            list<T> out;

            for (int i = 0; i < size1; i++) {
                T out_;
                for (int j = 0; j < size2; j++) {
                    auto val = this->get_element(i, j);

                    out_.push_back(val);
                }

                out.push_back(out_);
            }

            return out;
        }

        friend Array operator+(Array& matrix1, Array& matrix2) {
            auto [size_1_1, size_2_1] = matrix1.get_size();
            auto [size_1_2, size_2_2] = matrix2.get_size();

            assert(size_1_1 == size_1_2 && size_2_1 == size_2_2);

            Array matrix_out(matrix1);

            auto [size_1, size_2] = matrix_out.get_size();

            for (int i = 0; i < size_1; i++) {
                for (int j = 0; j < size_2; j++) {
                    matrix_out.get_element(i, j) += matrix2.get_element(i, j);
                }
            }

            return matrix_out;
        };
};

#endif
