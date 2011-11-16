#include <cassert>
#include <TlListMatrix.h>

TlListMatrix::TlListMatrix(const std::size_t reserveSize) : size_(0)
{
    this->elements_.clear();
}

TlListMatrix::~TlListMatrix()
{
}

void TlListMatrix::clear()
{
    this->size_ = 0;
    this->elements_.clear();
}

void TlListMatrix::add(const int row, const int col, const double value)
{
    this->elements_.push_back(Element(row, col, value));
    ++this->size_;
}

TlListMatrix::Element TlListMatrix::pop()
{
    Element answer(-1, -1, 0.0);
    if (this->size() > 0) {
        answer = this->elements_.front();
        this->elements_.pop_front();
        --this->size_;
    }

    return answer;
}
