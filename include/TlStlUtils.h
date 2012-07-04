#ifndef TLSTLUTILS_H
#define TLSTLUTILS_H

class TlStlUtils {
public:
    template<typename MapType, typename KeyArgType, typename ValueArgType>
    static typename MapType::iterator efficientAddorUpdate(MapType& m,
                                                           const KeyArgType& k,
                                                           const ValueArgType& v);
};


// see Effective STL 3-24
template<typename MapType, typename KeyArgType, typename ValueArgType>
typename MapType::iterator efficientAddorUpdate(MapType& m,
                                                const KeyArgType& k,
                                                const ValueArgType& v)
{
    typename MapType::iterator lb = m.lower_bound(k);
    if (lb != m.end() &&
        !(m.key_comp()(k, lb->first))) {
        lb->second = v;
        return lb;
    } else {
        typedef typename MapType::value_type MVT;
        return m.insert(lb, MVT(k, v));
    }
}

#endif // TLSTLUTILS_H
