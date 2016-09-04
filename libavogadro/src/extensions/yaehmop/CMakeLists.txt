include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(LINK_LIBS avogadro)
set(PLUGIN_LABEL extensions)
set(PLUGIN_TARGET extensions)

set(yaehmopextension_SRCS
  specialkpoints.cpp
  yaehmopbanddialog.cpp
  yaehmopextension.cpp
  yaehmopout.cpp
)

set(yaehmopextension_UIS
  yaehmopbanddialog.ui
)

avogadro_plugin_nogl(yaehmopextension
  "${yaehmopextension_SRCS}";
  "${yaehmopextension_UIS}"
)

# Download Yaehmop executable
# Let's set the name. Windows likes to add '.exe' at the end
if(WIN32)
  set(YAEHMOP_NAME "yaehmop.exe")
else(WIN32)
  set(YAEHMOP_NAME "yaehmop")
endif(WIN32)

if(NOT EXISTS ${CMAKE_BINARY_DIR}/bin/${YAEHMOP_NAME})
  # Linux
  if(UNIX AND NOT APPLE)
    set(YAEHMOP_DOWNLOAD_LOCATION "https://github.com/psavery/yaehmop/releases/download/3.0/linux64-yaehmop.tgz")
    set(MD5 "704ae01241b77fb08651d0ac786180d0")

  # Apple
  elseif(APPLE)
    set(YAEHMOP_DOWNLOAD_LOCATION "https://github.com/psavery/yaehmop/releases/download/3.0/mac64-yaehmop.tgz")
    set(MD5 "df37ca47554fb5c068eeac59dc5f8797")

  # Windows
  elseif(WIN32 AND NOT CYGWIN)
    set(YAEHMOP_DOWNLOAD_LOCATION "https://github.com/psavery/yaehmop/releases/download/3.0/win32-yaehmop.exe.tgz")
    set(MD5 "cea021178e82e26dadbe26453b8ada02")

  else()
    message(FATAL_ERROR "Yaehmop is not supported with the current OS type!")
  endif()

  message("-- Downloading Yaehmop executable from ${YAEHMOP_DOWNLOAD_LOCATION}")

  file(DOWNLOAD
       ${YAEHMOP_DOWNLOAD_LOCATION}
       "${CMAKE_BINARY_DIR}/bin/${YAEHMOP_NAME}.tgz"
       SHOW_PROGRESS
       EXPECTED_MD5 ${MD5})

  execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzvf ${YAEHMOP_NAME}.tgz
                  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

  file(REMOVE "${CMAKE_BINARY_DIR}/bin/${YAEHMOP_NAME}.tgz")

  install(FILES "${CMAKE_BINARY_DIR}/bin/${YAEHMOP_NAME}"
          DESTINATION ${BIN_INSTALL_DIR}
          PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                      GROUP_READ GROUP_EXECUTE
                      WORLD_READ WORLD_EXECUTE)

endif(NOT EXISTS ${CMAKE_BINARY_DIR}/bin/${YAEHMOP_NAME})
