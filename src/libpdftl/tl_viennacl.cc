#ifdef HAVE_CONFIG_H
#include "config.h"
#endif  // HAVE_CONFIG_H

#include <iostream>
#include <vector>
#include "viennacl/ocl/backend.hpp"

#include "TlUtils.h"
#include "tl_viennacl.h"

void TlViennaCL::setupAllAvailableDevices() {
    viennacl::ocl::platform pf;
    std::cout << "Platform info: " << pf.info() << std::endl;

    // const std::vector<viennacl::ocl::device> devices =
    //     pf.devices(CL_DEVICE_TYPE_DEFAULT);
    const std::vector<viennacl::ocl::device> devices =
        pf.devices(CL_DEVICE_TYPE_ALL);
    std::cout << "Number of devices: " << devices.size() << std::endl;

    std::vector<cl_device_id> device_id_array;
    for (std::size_t i = 0; i < devices.size(); ++i) {
        device_id_array.push_back(devices[i].id());
    }

    std::cout << "Creating context..." << std::endl;
    cl_int err;
    const cl_context my_context =
        clCreateContext(0, cl_uint(device_id_array.size()),
                        &(device_id_array[0]), NULL, NULL, &err);
    VIENNACL_ERR_CHECK(err);

    std::vector<cl_command_queue> queues(devices.size());
    for (std::size_t i = 0; i < devices.size(); ++i) {
#if OpenCL_VERSION_MAJOR >= 2
        {
            cl_queue_properties properties[] = {0};
            queues[i] = clCreateCommandQueueWithProperties(
                my_context, devices[i].id(), properties, &err);
        }
#else
        queues[i] = clCreateCommandQueue(my_context, devices[i].id(), 0, &err);
#endif
        VIENNACL_ERR_CHECK(err);
    }

    viennacl::ocl::setup_context(0, my_context, device_id_array, queues);
    viennacl::ocl::switch_context(
        0);  // activate the new context (only mandatory
             // with context-id not equal to zero)
}

std::string TlViennaCL::listDevices() {
    std::string answer = "";

    const std::vector<viennacl::ocl::device> devices =
        viennacl::ocl::current_context().devices();

    answer += TlUtils::format("the number of devices: %d\n", devices.size());
    const int numOfDevices = devices.size();
    for (int i = 0; i < numOfDevices; ++i) {
        answer += TlUtils::format("[%d] %s\n", i, devices[i].name().c_str());
    }

    return answer;
}

void TlViennaCL::switchDevice(const int id) {
    viennacl::ocl::current_context().switch_device(id);
}

std::string TlViennaCL::listCurrentDevice() {
    std::string answer = "";
    answer += TlUtils::format("current Device Name: %s\n",
                              viennacl::ocl::current_device().name().c_str());

    if (viennacl::ocl::current_device().double_support()) {
        answer += "this device supports double precision.\n";
    } else {
        answer += "this device does NOT support double precision.\n";
    }

    return answer;
}
