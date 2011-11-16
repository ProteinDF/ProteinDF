#ifndef TLLISTMATRIX_H
#define TLLISTMATRIX_H

#include <deque>

class TlListMatrix {
public:
    struct Element {
public:
        Element(int r = 0, int c = 0, double v = 0.0) : row(r), col(c), value(v) {
        }

public:
        int row;
        int col;
        double value;
    };

public:
    TlListMatrix(std::size_t reserveSize = 0);
    ~TlListMatrix();

public:
    void clear();

    std::size_t size() const {
        return this->size_;
    };

    void add(int row, int col, double value);

    Element pop();

private:
    std::size_t size_;
    std::deque<Element> elements_;
};

#endif // TLLISTMATRIX_H
