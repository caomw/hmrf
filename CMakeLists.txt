cmake_minimum_required(VERSION 2.8)
PROJECT(hiermrf)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

ADD_DEFINITIONS(-std=gnu++0x)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_49_0)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon-1.2.3)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon_build)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

LINK_DIRECTORIES (/home/sci/weiliu/packages/lemon_build/lemon)
LINK_DIRECTORIES (/home/sci/weiliu/packages/boost_1_49_0/build/lib)

ADD_EXECUTABLE(groupmrf
  groupmrfmain.cxx
  utility.cxx
  )

ADD_EXECUTABLE(label2ind
  label2ind.cxx
  )
ADD_EXECUTABLE(ind2label
  ind2label.cxx
  )

ADD_EXECUTABLE(sample2label
  sample2label.cxx
  )

ADD_EXECUTABLE(copyheader
  copyheader.cxx
  )

ADD_EXECUTABLE(applymask
  applymask.cxx
  )

ADD_EXECUTABLE(seedcorr
  seedcorr.cxx
  )

ADD_EXECUTABLE(compareclusterings
  compareclusterings.cxx
  )

ADD_EXECUTABLE(alignlabelmaps
  alignlabelmaps.cxx
  )


TARGET_LINK_LIBRARIES(groupmrf emon ${ITK_LIBRARIES} boost_program_options boost_filesystem)
TARGET_LINK_LIBRARIES(label2ind ${ITK_LIBRARIES} boost_program_options boost_filesystem)
TARGET_LINK_LIBRARIES(ind2label ${ITK_LIBRARIES} boost_program_options boost_filesystem)
TARGET_LINK_LIBRARIES(sample2label ${ITK_LIBRARIES} boost_program_options boost_filesystem)
TARGET_LINK_LIBRARIES(copyheader ${ITK_LIBRARIES} boost_program_options boost_filesystem)

TARGET_LINK_LIBRARIES(applymask ${ITK_LIBRARIES} boost_program_options boost_filesystem)
TARGET_LINK_LIBRARIES(seedcorr ${ITK_LIBRARIES} boost_program_options boost_filesystem)

TARGET_LINK_LIBRARIES(compareclusterings ${ITK_LIBRARIES} boost_program_options boost_filesystem)

TARGET_LINK_LIBRARIES(alignlabelmaps ${ITK_LIBRARIES} boost_program_options boost_filesystem)




