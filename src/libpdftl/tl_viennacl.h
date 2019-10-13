#ifndef TL_VIENNACL_H
#define TL_VIENNACL_H

#include <string>

class TlViennaCL {
   public:
    void setupAllAvailableDevices();
    std::string listDevices();

    void switchDevice(const int id);
    std::string listCurrentDevice();
};

#endif  // TL_VIENNACL_H
