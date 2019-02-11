include(ExternalProject)
find_program(GIT_EXECUTABLE git)

ExternalProject_Add(
  project_spdlog
  PREFIX ${CMAKE_BINARY_DIR}/spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.2.1
  TIMEOUT 01
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
  )

ExternalProject_Get_Property(project_spdlog source_dir)
set(SPDLOG_INCLUDE_DIR ${source_dir}/include)
