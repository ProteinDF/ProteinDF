#include <cassert>
#include "TlPosition.h"
#include "TlMatrix.h"

const double TlPosition::TOO_SMALL = 1.0E-20;
const double TlPosition::EQUIV_POSITION_RANGE = 1.0E-6;


TlPosition::TlPosition(const double x, const double y, const double z)
    : m_bHasDistance(false), m_distance(0.0)
{
    this->v_[0] = x;
    this->v_[1] = y;
    this->v_[2] = z;
}


TlPosition::TlPosition(const TlPosition& p)
    : m_bHasDistance(false), m_distance(0.0)
{
    this->v_[0] = p.x();
    this->v_[1] = p.y();
    this->v_[2] = p.z();
}


TlPosition::~TlPosition()
{
}


TlPosition& TlPosition::operator=(const TlPosition& rhs)
{
    if (this != &rhs) {
        this->v_[0] = rhs.x();
        this->v_[1] = rhs.y();
        this->v_[2] = rhs.z();
        this->m_bHasDistance = rhs.m_bHasDistance;
        this->m_distance     = rhs.m_distance;
    }

    return *this;
}


void TlPosition::clear()
{
    this->v_[0] = 0.0;
    this->v_[1] = 0.0;
    this->v_[2] = 0.0;
    this->m_bHasDistance = false;
}


TlPosition& TlPosition::rotateBy(const TlMatrix& rot)
{
    assert(rot.getNumOfRows() == 3);
    assert(rot.getNumOfCols() == 3);

    const double x = rot(0, 0) * this->x() + rot(0, 1) * this->y() + rot(0, 2) * this->z();
    const double y = rot(1, 0) * this->x() + rot(1, 1) * this->y() + rot(1, 2) * this->z();
    const double z = rot(2, 0) * this->x() + rot(2, 1) * this->y() + rot(2, 2) * this->z();

    this->v_[0] = x;
    this->v_[1] = y;
    this->v_[2] = z;
    this->m_bHasDistance = false;

    return *this;
}


TlPosition& TlPosition::operator+=(const TlPosition& rhs)
{
    this->shiftBy(rhs);

    return *this;
}


TlPosition& TlPosition::operator-=(const TlPosition& rhs)
{
    this->shiftBy(-1.0 * rhs);

    return *this;
}


TlPosition& TlPosition::operator*=(const double a)
{
    this->v_[0] *= a;
    this->v_[1] *= a;
    this->v_[2] *= a;
    this->m_bHasDistance = false;

    return *this;
}


TlPosition& TlPosition::operator/=(const double a)
{
    assert(a != 0.0);

    const double div_a = 1.0 / a;
    this->v_[0] *= div_a;
    this->v_[1] *= div_a;
    this->v_[2] *= div_a;
    this->m_bHasDistance = false;

    return *this;
}


bool TlPosition::operator==(const TlPosition& rhs) const
{
    return this->isWithinRange(rhs, EQUIV_POSITION_RANGE);
}


bool TlPosition::isWithinRange(const TlPosition& pos, double range) const
{
    const double dx = std::fabs(this->x() - pos.x());
    if (dx > range) {
        return false;
    }

    const double dy = std::fabs(this->y() - pos.y());
    if (dy > range) {
        return  false;
    }

    const double dz = std::fabs(this->z() - pos.z());
    if (dz > range) {
        return false;
    }

    const double dd = dx*dx + dy*dy + dz*dz;
    if (dd > range) {
        return false;
    }

    return true;
}


void TlPosition::unit()
{
    const double s = this->distanceFrom();
    if (s > 0) {
        *this /= s;
    }
}

