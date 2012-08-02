#ifndef TLSHAREDPOINTER_H
#define TLSHAREDPOINTER_H

#include <cassert>

/** ��ե���󥹥�����ȼ�Smart Pointer
 *
 *  - new�Ǻ��������ݥ��󥿤Τ��ݻ������뤳�ȡ�new �ʳ������������ݥ��󥿤ϲ����Ǥ��ʤ���
 *  - new�Ǻ��������ݥ��󥿤�ʣ����TlSharedPoinnter���Ϥ��ʤ�����(2�Ų����ˤʤ�)��
 *  - ����Υݥ��󥿤��ݻ������ʤ�����(delete []��Ԥ�ʤ�����)��
 */
template<typename T>
class TlSharedPointer {
public:
    TlSharedPointer(T* pObject);
    TlSharedPointer(const TlSharedPointer<T>& psp);

    ~TlSharedPointer();

    TlSharedPointer<T>& operator=(const TlSharedPointer<T>& sp);

    bool operator==(const T* p);
    bool operator!=(const T* p);

    T* operator->() const;
    T& operator*() const;

private:
    void set(const TlSharedPointer<T>& psp);
    void release();

private:
    T* m_pObject;
    int* m_pCount;

    const static T* m_pNullObject; // NULL�ݥ���
};

////////////////////////////////////////////////////////////////////////////////

template<typename T>
const T* TlSharedPointer<T>::m_pNullObject = NULL;

template<typename T>
TlSharedPointer<T>::TlSharedPointer(T* pObject) {
    this->m_pObject = pObject;

    this->m_pCount = new int;
    *(this->m_pCount) = 1;
}

template<typename T>
TlSharedPointer<T>::TlSharedPointer(const TlSharedPointer<T>& psp) {
    this->m_pObject = NULL;
    this->m_pCount = NULL;

    this->set(psp);
}

template<typename T>
TlSharedPointer<T>::~TlSharedPointer() {
    this->release();
}

template<typename T>
TlSharedPointer<T>& TlSharedPointer<T>::operator=(const TlSharedPointer<T>& psp) {
    this->set(psp);
    return (*this);
}

template<typename T>
bool TlSharedPointer<T>::operator==(const T* p) {
    if (this->m_pObject == p) {
        return true;
    } else {
        return false;
    }
}

template<typename T>
bool TlSharedPointer<T>::operator!=(const T* p) {
    return !(this->operator==(p));
}

template<typename T>
T* TlSharedPointer<T>::operator->() const {
    return this->m_pObject;
}

template<typename T>
T& TlSharedPointer<T>::operator*() const {
    return *(this->m_pObject);
}

template<typename T>
void TlSharedPointer<T>::set(const TlSharedPointer<T>& psp){
    if (this != &psp){
        this->release();

        this->m_pObject = psp.m_pObject;
        this->m_pCount = psp.m_pCount;
        ++(*(this->m_pCount));
    }
}

template<typename T>
void TlSharedPointer<T>::release(){
    if (this->m_pCount != NULL){
        --(*(this->m_pCount));

        if (*(this->m_pCount) == 0){
            delete this->m_pObject;
            delete this->m_pCount;

            this->m_pObject = NULL;
            this->m_pCount = NULL;
        }
    }
}

#endif // TLSHAREDPOINTER_H
