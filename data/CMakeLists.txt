cmake_minimum_required(VERSION 3.5)

set(DATA_FILES
  basis2
  pdfkwd.mpac
  alias.mpac
  harris.mpac
  )

install(FILES ${DATA_FILES} DESTINATION data)

