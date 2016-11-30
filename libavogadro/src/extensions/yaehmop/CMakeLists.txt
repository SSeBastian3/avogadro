include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(LINK_LIBS avogadro)
set(PLUGIN_LABEL extensions)
set(PLUGIN_TARGET extensions)

set(yaehmopextension_SRCS
  specialkpoints.cpp
  yaehmopbanddialog.cpp
  yaehmopextension.cpp
  yaehmopprojecteddosdialog.cpp
  yaehmoptotaldosdialog.cpp
  yaehmopout.cpp
)

set(yaehmopextension_UIS
  yaehmopbanddialog.ui
  yaehmopprojecteddosdialog.ui
  yaehmoptotaldosdialog.ui
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

option(USE_SYSTEM_YAEHMOP "Use system YAeHMOP" OFF)

if(NOT ${USE_SYSTEM_YAEHMOP} AND
   NOT EXISTS ${CMAKE_BINARY_DIR}/bin/${YAEHMOP_NAME})
   set(YAEHMOP_V "v3.0.1-bug-fix")
  # Linux
  if(UNIX AND NOT APPLE)
    set(YAEHMOP_DOWNLOAD_LOCATION "https://github.com/psavery/yaehmop/releases/download/${YAEHMOP_V}/linux64-yaehmop.tgz")
    set(MD5 "92a2738f36b9362a2153a6798574de0e")

  # Apple
  elseif(APPLE)
    set(YAEHMOP_DOWNLOAD_LOCATION "https://github.com/psavery/yaehmop/releases/download/${YAEHMOP_V}/mac64-yaehmop.tgz")
    set(MD5 "06a07e17e9d50251c292b18048eff7d0")

  # Windows
  elseif(WIN32 AND NOT CYGWIN)
    set(YAEHMOP_DOWNLOAD_LOCATION "https://github.com/psavery/yaehmop/releases/download/${YAEHMOP_V}/win32-yaehmop.exe.tgz")
    set(MD5 "a0d0ecd7453f58ac8ff5f898de208ace")

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
