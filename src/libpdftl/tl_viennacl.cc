#include "tl_viennacl.h"
#include <iostream>
#include <vector>
#include "config.h"
#include "viennacl/ocl/backend.hpp"

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
      clCreateContext(0, cl_uint(device_id_array.size()), &(device_id_array[0]),
                      NULL, NULL, &err);
  VIENNACL_ERR_CHECK(err);

  cl_queue_properties properties[] = {0};
  std::vector<cl_command_queue> queues(devices.size());
  for (std::size_t i = 0; i < devices.size(); ++i) {
    //queues[i] = clCreateCommandQueue(my_context, devices[i].id(), 0, &err);
    queues[i] = clCreateCommandQueueWithProperties(my_context, devices[i].id(), properties, &err);
    VIENNACL_ERR_CHECK(err);
  }

  viennacl::ocl::setup_context(0, my_context, device_id_array, queues);
  viennacl::ocl::switch_context(0);  // activate the new context (only mandatory
                                     // with context-id not equal to zero)
}

void TlViennaCL::showDevices() {
  const std::vector<viennacl::ocl::device> devices =
      viennacl::ocl::current_context().devices();

  std::cout << "# devices: " << devices.size() << std::endl;
  for (int i = 0; i < devices.size(); ++i) {
    std::cout << i << ":" << devices[i].name() << std::endl;
  }
}

void TlViennaCL::switchDevice(const int id) {
  viennacl::ocl::current_context().switch_device(id);
}

void TlViennaCL::showCurrentDevice() {
  std::cout << "current Device Name: " << viennacl::ocl::current_device().name()
            << std::endl;

  if (viennacl::ocl::current_device().double_support()) {
      std::cout << "this device supports double precision." << std::endl;
  } else {
      std::cout << "this device does NOT support double precision." << std::endl;
  }
}
