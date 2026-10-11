# cmake_minimum_required(VERSION 3.5)
add_subdirectory(unit_test)

find_program(YML2MPAC_EXECUTABLE NAMES yml2mpac yml2mpac.py)
if(NOT YML2MPAC_EXECUTABLE)
  message(FATAL_ERROR
    "yml2mpac (or yml2mpac.py) not found. "
    "Please install ProteinDF_bridge in your Python environment or specify -DYML2MPAC_EXECUTABLE=/path/to/yml2mpac."
  )
endif()

set(PDFKWD_YAML "${CMAKE_CURRENT_SOURCE_DIR}/pdfkwd.yaml")
set(PDFKWD_MPAC "${CMAKE_CURRENT_BINARY_DIR}/pdfkwd.mpac")

add_custom_command(
  OUTPUT "${PDFKWD_MPAC}"
  COMMAND "${YML2MPAC_EXECUTABLE}" "${PDFKWD_YAML}" "${PDFKWD_MPAC}"
  DEPENDS "${PDFKWD_YAML}"
  COMMENT "Generating ${PDFKWD_MPAC} from ${PDFKWD_YAML}"
  VERBATIM
)

add_custom_target(generate_pdfkwd_mpac ALL DEPENDS "${PDFKWD_MPAC}")

set(DATA_FILES
  basis2
  alias.mpac
  harris.mpac
  )

install(FILES ${DATA_FILES} DESTINATION data)
install(FILES "${PDFKWD_MPAC}" DESTINATION data)
