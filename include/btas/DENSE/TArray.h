#ifndef __BTAS_DENSE_TARRAY_H
#define __BTAS_DENSE_TARRAY_H 1

#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>

namespace btas {

   template<typename T, size_t N> class TArray;

   /*! \class TArray
    *  \brief Dense array class
    *
    *  Fixed-rank array class implemented in terms of std::array and std::vector
    *  Since using C++11 features, not compatible for C++03 compiler
    *
    *  \param T value type
    *  \param N array rank
    */
   template<typename T, size_t N>
      class TArray {
         private:

            friend class boost::serialization::access;

            //! Any TArray classes being friend of TArray<T, N>
            template<typename U, size_t M>
               friend class TArray;

            //! Enables to use boost serialization
            template<class Archive>
               void serialize(Archive& ar, const unsigned int version) { ar & m_shape & m_stride & m_store; }

         public:

            //! TArray<T, N>::iterator
            typedef typename std::vector<T>::iterator       iterator;

            //! TArray<T, N>::const_iterator
            typedef typename std::vector<T>::const_iterator const_iterator;

            //####################################################################################################
            // Member Functions
            //####################################################################################################

            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Constructor, Destructor, and Assignment operators
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            //! default constructor
            TArray() : m_shape(uniform<int, N>(0)), m_stride(uniform<int, N>(0)), m_store(new std::vector<T>()) { }

            //! destructor
            ~TArray() { }

            //! copy constructor
            TArray(const TArray& other) : m_store(new std::vector<T>()) {
               copy(other);
            }

            //! copy assignment operator
            TArray& operator= (const TArray& other) {
               copy(other);
               return *this;
            }

            //! addition assignment operator
            TArray& operator+=(const TArray& other) {
               add (other);
               return *this;
            }

            //! copying from other to this
            void copy(const TArray& other) {
               m_shape  = other.m_shape;
               m_stride = other.m_stride;
               Copy(*other.m_store, *m_store);
            }

            //! return copy of this
            TArray copy() const {
               TArray _cpy;
               _cpy.copy(*this);
               return std::move(_cpy);
            }

            //! scale by const value
            void scale(const T& alpha) {
               Scal(alpha, *m_store);
            }

            //! adding  from other to this
            void add (const TArray& other) {
               assert(m_shape  == other.m_shape);
               assert(m_stride == other.m_stride);
               Axpy(static_cast<T>(1), *other.m_store, *m_store);
            }

            //! move constructor
            explicit TArray(TArray&& other)
               : m_shape(std::move(other.m_shape)), m_stride(std::move(other.m_stride)), m_store(std::move(other.m_store) ) {

                  //make sure the other still point to something, else it will give errors when going out of scope.
                  other.m_store = shared_ptr< std::vector<T> >(new std::vector<T>());
                  other.m_shape = uniform<int, N>(0);
                  other.m_stride = uniform<int, N>(0);

               }

            //! move assignment
            TArray& operator= (TArray&& other) {

               this->swap(other);

               return *this;
            }

            //! take data reference from other
            void reference(const TArray& other) {
               m_shape  = other.m_shape;
               m_stride = other.m_stride;
               m_store  = other.m_store;
            }

            //! return data reference of this
            TArray reference() const {
               TArray _ref;
               _ref.reference(*this);
               return std::move(_ref);
            }

            //! convenient constructor with array shape, for N = 1
            explicit TArray(int n01) : m_store(new std::vector<T>()) {
               resize(n01);
            }

            //! convenient constructor with array shape, for N = 2
            TArray(int n01, int n02) : m_store(new std::vector<T>()) {
               resize(n01, n02);
            }

            //! convenient constructor with array shape, for N = 3
            TArray(int n01, int n02, int n03) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03);
            }

            //! convenient constructor with array shape, for N = 4
            TArray(int n01, int n02, int n03, int n04) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04);
            }

            //! convenient constructor with array shape, for N = 5
            TArray(int n01, int n02, int n03, int n04, int n05) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04, n05);
            }

            //! convenient constructor with array shape, for N = 6
            TArray(int n01, int n02, int n03, int n04, int n05, int n06) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04, n05, n06);
            }

            //! convenient constructor with array shape, for N = 7
            TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04, n05, n06, n07);
            }

            //! convenient constructor with array shape, for N = 8
            TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04, n05, n06, n07, n08);
            }

            //! convenient constructor with array shape, for N = 9
            TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04, n05, n06, n07, n08, n09);
            }

            //! convenient constructor with array shape, for N = 10
            TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04, n05, n06, n07, n08, n09, n10);
            }

            //! convenient constructor with array shape, for N = 11
            TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11);
            }

            //! convenient constructor with array shape, for N = 12
            TArray(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11, int n12) : m_store(new std::vector<T>()) {
               resize(n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11, n12);
            }

            //! convenient constructor with array shape, for arbitrary N
            TArray(const IVector<N>& _shape) : m_store(new std::vector<T>()) {
               resize(_shape);
            }

            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Resizing functions
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            //! resize array shape, for N = 1
            void resize(int n01) {
               IVector< 1> _shape = { n01 };
               resize(_shape);
            }

            //! resize array shape, for N = 2
            void resize(int n01, int n02) {
               IVector< 2> _shape = { n01, n02 };
               resize(_shape);
            }

            //! resize array shape, for N = 3
            void resize(int n01, int n02, int n03) {
               IVector< 3> _shape = { n01, n02, n03 };
               resize(_shape);
            }

            //! resize array shape, for N = 4
            void resize(int n01, int n02, int n03, int n04) {
               IVector< 4> _shape = { n01, n02, n03, n04 };
               resize(_shape);
            }

            //! resize array shape, for N = 5
            void resize(int n01, int n02, int n03, int n04, int n05) {
               IVector< 5> _shape = { n01, n02, n03, n04, n05 };
               resize(_shape);
            }

            //! resize array shape, for N = 6
            void resize(int n01, int n02, int n03, int n04, int n05, int n06) {
               IVector< 6> _shape = { n01, n02, n03, n04, n05, n06 };
               resize(_shape);
            }

            //! resize array shape, for N = 7
            void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07) {
               IVector< 7> _shape = { n01, n02, n03, n04, n05, n06, n07 };
               resize(_shape);
            }

            //! resize array shape, for N = 8
            void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08) {
               IVector< 8> _shape = { n01, n02, n03, n04, n05, n06, n07, n08 };
               resize(_shape);
            }

            //! resize array shape, for N = 9
            void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09) {
               IVector< 9> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09 };
               resize(_shape);
            }

            //! resize array shape, for N = 10
            void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10) {
               IVector<10> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10 };
               resize(_shape);
            }

            //! resize array shape, for N = 11
            void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11) {
               IVector<11> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11 };
               resize(_shape);
            }

            //! resize array shape, for N = 12
            void resize(int n01, int n02, int n03, int n04, int n05, int n06, int n07, int n08, int n09, int n10, int n11, int n12) {
               IVector<12> _shape = { n01, n02, n03, n04, n05, n06, n07, n08, n09, n10, n11, n12 };
               resize(_shape);
            }

            //! resize array by _shape, for arbitrary N
            /*! detects rank-mismatching error at compilation time */
            void resize(const IVector<N>& _shape) {
               m_shape = _shape;
               // calculate stride
               size_t stride = 1;
               for(int i = N-1; i >= 0; --i) {
                  m_stride[i] = stride;
                  stride *= m_shape[i];
               }
               // allocate memory
               m_store->resize(stride);
               return;
            }

            /// return reference reshaped: this object is cleared
            template<size_t M>
               TArray<T, M> reshape_clear(const IVector<M>& shape_)
               {
                  TArray<T, M> x; 

                  x.m_shape = shape_;

                  // calculate stride
                  size_t stride = 1;

                  for(int i = M-1; i >= 0; --i) {
                     x.m_stride[i] = stride;
                     stride *= x.m_shape[i];
                  }

                  assert(stride == this->size());

                  x.m_store = std::move(this->m_store);

                  this->m_store = shared_ptr< std::vector<T> >(new std::vector<T>());
                  this->m_shape = uniform<int, N>(0);
                  this->m_stride = uniform<int, N>(0);

                  return x;

               }

            /// return shared reference reshaped
            /// implaced reshape can be done as
            /// `y.swap(x.reshape(shape(n1,n2,...)));`
            template<size_t M>
               TArray<T, M> reshape (const IVector<M>& shape_)
               {
                  TArray<T, M> x(shape_);

                  assert(x.size() == this->size());

                  x.m_store = this->m_store; // shallow copy

                  return x;
               }

            /// return deep copy reshaped
            template<size_t M>
               TArray<T, M> reshape (const IVector<M>& shape_) const
               {
                  TArray<T, M> x(shape_);

                  assert(x.size() == this->size());

                  Copy(*this->m_store, *x.m_store); // deep copy

                  return x;
               }

            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // Data Accessing
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

            //! returns first iterator position (const)
            const_iterator begin() const { return m_store->begin(); }

            //! returns first iterator position
            iterator begin()       { return m_store->begin(); }

            //! returns last iterator position (const)
            const_iterator end() const { return m_store->end(); }

            //! returns last iterator position
            iterator end()       { return m_store->end(); }

            //! returns array shape
            const IVector<N>& shape() const { return m_shape; }

            //! returns array shape for rank i
            int shape(int i) const { return m_shape[i]; }

            //! returns array stride
            const IVector<N>& stride() const { return m_stride; }

            //! returns array stride for rank i
            int stride(int i) const { return m_stride[i]; }

            //! returns allocated size
            size_t size() const { return m_store->size(); }

            //! returns array element (N = 1) without range check
            const T& operator() (int i01) const {
               IVector< 1> _index = { i01 };
               return operator()(_index);
            }

            //! returns array element (N = 2) without range check
            const T& operator() (int i01, int i02) const {
               IVector< 2> _index = { i01, i02 };
               return operator()(_index);
            }

            //! returns array element (N = 3) without range check
            const T& operator() (int i01, int i02, int i03) const {
               IVector< 3> _index = { i01, i02, i03 };
               return operator()(_index);
            }

            //! returns array element (N = 4) without range check
            const T& operator() (int i01, int i02, int i03, int i04) const {
               IVector< 4> _index = { i01, i02, i03, i04 };
               return operator()(_index);
            }

            //! returns array element (N = 5) without range check
            const T& operator() (int i01, int i02, int i03, int i04, int i05) const {
               IVector< 5> _index = { i01, i02, i03, i04, i05 };
               return operator()(_index);
            }

            //! returns array element (N = 6) without range check
            const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06) const {
               IVector< 6> _index = { i01, i02, i03, i04, i05, i06 };
               return operator()(_index);
            }

            //! returns array element (N = 7) without range check
            const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07) const {
               IVector< 7> _index = { i01, i02, i03, i04, i05, i06, i07 };
               return operator()(_index);
            }

            //! returns array element (N = 8) without range check
            const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08) const {
               IVector< 8> _index = { i01, i02, i03, i04, i05, i06, i07, i08 };
               return operator()(_index);
            }

            //! returns array element (N = 9) without range check
            const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09) const {
               IVector< 9> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09 };
               return operator()(_index);
            }

            //! returns array element (N = 10) without range check
            const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10) const {
               IVector<10> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10 };
               return operator()(_index);
            }

            //! returns array element (N = 11) without range check
            const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11) const {
               IVector<11> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11 };
               return operator()(_index);
            }

            //! returns array element (N = 12) without range check
            const T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12) const {
               IVector<12> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12 };
               return operator()(_index);
            }

            //! returns array element (arbitrary N) without range check
            const T& operator() (const IVector<N>& _index) const {
               return (*m_store)[dot(_index, m_stride)];
            }

            //! returns array element (N = 1) without range check
            T& operator() (int i01) {
               IVector< 1> _index = { i01 };
               return operator()(_index);
            }

            //! returns array element (N = 2) without range check
            T& operator() (int i01, int i02) {
               IVector< 2> _index = { i01, i02 };
               return operator()(_index);
            }

            //! returns array element (N = 3) without range check
            T& operator() (int i01, int i02, int i03) {
               IVector< 3> _index = { i01, i02, i03 };
               return operator()(_index);
            }

            //! returns array element (N = 4) without range check
            T& operator() (int i01, int i02, int i03, int i04) {
               IVector< 4> _index = { i01, i02, i03, i04 };
               return operator()(_index);
            }

            //! returns array element (N = 5) without range check
            T& operator() (int i01, int i02, int i03, int i04, int i05) {
               IVector< 5> _index = { i01, i02, i03, i04, i05 };
               return operator()(_index);
            }

            //! returns array element (N = 6) without range check
            T& operator() (int i01, int i02, int i03, int i04, int i05, int i06) {
               IVector< 6> _index = { i01, i02, i03, i04, i05, i06 };
               return operator()(_index);
            }

            //! returns array element (N = 7) without range check
            T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07) {
               IVector< 7> _index = { i01, i02, i03, i04, i05, i06, i07 };
               return operator()(_index);
            }

            //! returns array element (N = 8) without range check
            T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08) {
               IVector< 8> _index = { i01, i02, i03, i04, i05, i06, i07, i08 };
               return operator()(_index);
            }

            //! returns array element (N = 9) without range check
            T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09) {
               IVector< 9> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09 };
               return operator()(_index);
            }

            //! returns array element (N = 10) without range check
            T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10) {
               IVector<10> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10 };
               return operator()(_index);
            }

            //! returns array element (N = 11) without range check
            T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11) {
               IVector<11> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11 };
               return operator()(_index);
            }

            //! returns array element (N = 12) without range check
            T& operator() (int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12) {
               IVector<12> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12 };
               return operator()(_index);
            }

            //! returns array element (arbitrary N) without range check
            T& operator() (const IVector<N>& _index) {
               return (*m_store)[dot(_index, m_stride)];
            }

            //! returns array element (N = 1) with range check
            const T& at(int i01) const {
               IVector< 1> _index = { i01 };
               return at(_index);
            }

            //! returns array element (N = 2) with range check
            const T& at(int i01, int i02) const {
               IVector< 2> _index = { i01, i02 };
               return at(_index);
            }

            //! returns array element (N = 3) with range check
            const T& at(int i01, int i02, int i03) const {
               IVector< 3> _index = { i01, i02, i03 };
               return at(_index);
            }

            //! returns array element (N = 4) with range check
            const T& at(int i01, int i02, int i03, int i04) const {
               IVector< 4> _index = { i01, i02, i03, i04 };
               return at(_index);
            }

            //! returns array element (N = 5) with range check
            const T& at(int i01, int i02, int i03, int i04, int i05) const {
               IVector< 5> _index = { i01, i02, i03, i04, i05 };
               return at(_index);
            }

            //! returns array element (N = 6) with range check
            const T& at(int i01, int i02, int i03, int i04, int i05, int i06) const {
               IVector< 6> _index = { i01, i02, i03, i04, i05, i06 };
               return at(_index);
            }

            //! returns array element (N = 7) with range check
            const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07) const {
               IVector< 7> _index = { i01, i02, i03, i04, i05, i06, i07 };
               return at(_index);
            }

            //! returns array element (N = 8) with range check
            const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08) const {
               IVector< 8> _index = { i01, i02, i03, i04, i05, i06, i07, i08 };
               return at(_index);
            }

            //! returns array element (N = 9) with range check
            const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09) const {
               IVector< 9> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09 };
               return at(_index);
            }

            //! returns array element (N = 10) with range check
            const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10) const {
               IVector<10> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10 };
               return at(_index);
            }

            //! returns array element (N = 11) with range check
            const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11) const {
               IVector<11> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11 };
               return at(_index);
            }

            //! returns array element (N = 12) with range check
            const T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12) const {
               IVector<12> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12 };
               return at(_index);
            }

            //! returns array element (arbitrary N) with range check
            const T& at(const IVector<N>& _index) const {
               return m_store->at(dot(_index, m_stride));
            }

            //! returns array element (N = 1) with range check
            T& at(int i01) {
               IVector< 1> _index = { i01 };
               return at(_index);
            }

            //! returns array element (N = 2) with range check
            T& at(int i01, int i02) {
               IVector< 2> _index = { i01, i02 };
               return at(_index);
            }

            //! returns array element (N = 3) with range check
            T& at(int i01, int i02, int i03) {
               IVector< 3> _index = { i01, i02, i03 };
               return at(_index);
            }

            //! returns array element (N = 4) with range check
            T& at(int i01, int i02, int i03, int i04) {
               IVector< 4> _index = { i01, i02, i03, i04 };
               return at(_index);
            }

            //! returns array element (N = 5) with range check
            T& at(int i01, int i02, int i03, int i04, int i05) {
               IVector< 5> _index = { i01, i02, i03, i04, i05 };
               return at(_index);
            }

            //! returns array element (N = 6) with range check
            T& at(int i01, int i02, int i03, int i04, int i05, int i06) {
               IVector< 6> _index = { i01, i02, i03, i04, i05, i06 };
               return at(_index);
            }

            //! returns array element (N = 7) with range check
            T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07) {
               IVector< 7> _index = { i01, i02, i03, i04, i05, i06, i07 };
               return at(_index);
            }

            //! returns array element (N = 8) with range check
            T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08) {
               IVector< 8> _index = { i01, i02, i03, i04, i05, i06, i07, i08 };
               return at(_index);
            }

            //! returns array element (N = 9) with range check
            T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09) {
               IVector< 9> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09 };
               return at(_index);
            }

            //! returns array element (N = 10) with range check
            T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10) {
               IVector<10> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10 };
               return at(_index);
            }

            //! returns array element (N = 11) with range check
            T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11) {
               IVector<11> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11 };
               return at(_index);
            }

            //! returns array element (N = 12) with range check
            T& at(int i01, int i02, int i03, int i04, int i05, int i06, int i07, int i08, int i09, int i10, int i11, int i12) {
               IVector<12> _index = { i01, i02, i03, i04, i05, i06, i07, i08, i09, i10, i11, i12 };
               return at(_index);
            }

            //! returns array element (arbitrary N) with range check
            T& at(const IVector<N>& _index) {
               return m_store->at(dot(_index, m_stride));
            }

            //! returns the first pointer of array elements
            const T* data() const { return m_store->data(); }

            //! returns the first pointer of array elements
            T* data()       { return m_store->data(); }


            //! fills elements by constant value
            void fill(const T& val) {
               std::fill(m_store->begin(), m_store->end(), val);
            }

            //! fills elements by constant value
            void operator= (const T& val) { fill(val); }

            //! generates array elements by function gen
            /*! Generator is either default constructible class or function pointer, which can be called by gen() */
            template<class Generator>
               void generate(Generator gen) {
                  std::generate(m_store->begin(), m_store->end(), gen);
               }

            //! deallocate storage
            void clear() {
               m_shape = uniform<int, N>(0);
               m_stride = uniform<int, N>(0);
               m_store->clear();
            }

            /// swap object
            void swap (TArray& x)
            {
               std::swap(this->m_shape,  x.m_shape);
               std::swap(this->m_stride, x.m_stride);
               this->m_store.swap(x.m_store);
            }

            int use_count(){

               return m_store.use_count();

            }

         private:

            //####################################################################################################
            // Member Variables
            //####################################################################################################

            //! array shape
            IVector<N>
               m_shape;

            //! array stride
            IVector<N>
               m_stride;

            //! array storage
            shared_ptr<std::vector<T>>
               m_store;

      }; // class TArray

}; // namespace btas

#endif // __BTAS_DENSE_TARRAY_H
