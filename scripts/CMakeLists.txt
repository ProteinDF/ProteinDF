cmake_minimum_required(VERSION 3.5)

set(INSTALL_SCRIPTS
  pdf-clean.sh
  pdf-setup.sh
  pdf-printenv.sh
  pdf-serial.sh
  pdf-parallel.sh
  )

install(PROGRAMS ${INSTALL_SCRIPTS}
  DESTINATION bin
  )

