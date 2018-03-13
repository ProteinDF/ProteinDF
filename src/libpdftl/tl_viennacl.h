#ifndef TL_VIENNACL_H
#define TL_VIENNACL_H

class TlViennaCL {
 public:
    void setupAllAvailableDevices();
    void showDevices();

    void switchDevice(const int id);
    void showCurrentDevice();
};

#endif  // TL_VIENNACL_H
