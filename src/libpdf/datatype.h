#ifndef DATATYPE_H
#define DATATYPE_H

#include <vector>

typedef int index_type;

class Index2 {
   public:
    explicit Index2(index_type i1 = 0, index_type i2 = 0)
        : index1_(i1), index2_(i2) {}

    bool operator<(const Index2& rhs) const {
        if (this->index1_ < rhs.index1_) {
            return true;
        } else if (this->index1_ == rhs.index1_) {
            return (this->index2_ < rhs.index2_);
        }

        return false;
    }

    bool operator==(const Index2& rhs) const {
        return ((this->index1_ == rhs.index1_) &&
                (this->index2_ == rhs.index2_));
    }

    index_type index1() const { return this->index1_; }

    index_type index2() const { return this->index2_; }

   protected:
    index_type index1_;
    index_type index2_;
};

typedef std::vector<Index2> PQ_PairArray;

#endif  // DATATYPE_H
