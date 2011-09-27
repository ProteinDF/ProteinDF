#ifndef TLSIMPLEVECTOR_H
#define TLSIMPLEVECTOR_H

#include <cassert>
#include <algorithm>

template<typename T>
class TlSimpleVector {
public:
    explicit TlSimpleVector(size_t nSize =0, const T& value = T()) : m_nSize(nSize), m_pData(NULL) {
        this->m_pData = new T[m_nSize];
        //     const size_t nEnd = m_nSize;
        //     for (size_t i = 0; i < nEnd; ++i){
        //       this->m_pData[i] = value;
        //     }
        std::fill(this->m_pData, (this->m_pData + this->m_nSize), value);
    }

    TlSimpleVector(const TlSimpleVector<T>& rhs) : m_nSize(rhs.m_nSize), m_pData(NULL) {
        this->m_pData = new T[m_nSize];
        //     const size_t nEnd = m_nSize;
        //     for (size_t i = 0; i < nEnd; ++i){
        //       pData[i] = rhs.m_pData[i];
        //     }
        std::copy(rhs.m_pData, (rhs.m_pData +rhs.m_nSize), this->m_pData);
    }

    ~TlSimpleVector() {
        delete [] this->m_pData;
        this->m_pData = NULL;
    }

    TlSimpleVector& operator=(const TlSimpleVector<T>& rhs) {
        if (this != &rhs) {
            if (this->m_pData != NULL) {
                delete [] this->m_pData;
            }

            m_nSize = rhs.m_nSize;
            this->m_pData = new T[m_nSize];
            //       const size_t nEnd = m_nSize;
            //       for (size_t i = 0; i < nEnd; ++i){
            //    pData[i] = rhs.pData[i];
            //       }
            std::copy(rhs.m_pData, (rhs.m_pData +rhs.m_nSize), this->m_pData);
        }

        return *this;
    }

    TlSimpleVector& operator+=(const TlSimpleVector<T>& rhs) {
        assert(m_nSize == rhs.m_nSize);

        const size_t nEnd = m_nSize;
        for (size_t i = 0; i < nEnd; ++i) {
            this->m_pData[i] += rhs.pData[i];
        }

        return *this;
    }


    T& operator[](size_t index) {
        assert(index < m_nSize);
        return this->m_pData[index];
    }

    const T& operator[](size_t index) const {
        assert(index < m_nSize);
        return this->m_pData[index];
    }

    size_t size() const {
        return this->m_nSize;
    }

    T max() const {
        assert(m_nSize > 0);
        T answer = this->m_pData[0];
        const size_t nEnd = m_nSize;
        for (size_t i = 1; i < nEnd; ++i) {
            answer = std::max(answer, this->m_pData[i]);
        }

        return answer;
    }

public:
    friend TlSimpleVector<T>& operator+(const TlSimpleVector<T>& x, TlSimpleVector<T>& y) {
        TlSimpleVector<T> answer = x;
        answer += y;

        return answer;
    }

private:
    size_t m_nSize;
    T* m_pData;
};

#endif // TLSIMPLEVECTOR_H

