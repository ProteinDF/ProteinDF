#include "config.h"

#include "TlGetopt.h"
#include "gtest/gtest.h"

#ifdef HAVE_VIENNACL
#include "tl_viennacl.h"
#endif  // HAVE_VIENNACL

int main(int argc, char **argv) {
    TlGetopt opt(argc, argv, "d:");
    ::testing::InitGoogleTest(&argc, argv);

// ViennaCL
#ifdef HAVE_VIENNACL
    {
        int deviceId = 0;
        if (!opt["d"].empty()) {
            deviceId = std::atoi(opt["d"].c_str());
        }

        TlViennaCL vcl;
        vcl.setupAllAvailableDevices();
        std::cout << vcl.listDevices() << std::endl;

        vcl.switchDevice(deviceId);
        std::cout << vcl.listCurrentDevice() << std::endl;
    }
#endif  // HAVE_VIENNACL

    return RUN_ALL_TESTS();
}
