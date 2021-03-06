cmake_minimum_required(VERSION 2.8)
PROJECT(groupmrf)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE OFF)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/boost_1_53_0)
# INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon-1.2.3)
# INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon-1.3)
# INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon_build)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon-1.3)
INCLUDE_DIRECTORIES (/home/sci/weiliu/packages/lemon-1.3/build/wukong)
INCLUDE_DIRECTORIES (${PROJECT_SOURCE_DIR})

# LINK_DIRECTORIES (/home/sci/weiliu/packages/lemon_build/lemon)
LINK_DIRECTORIES (/home/sci/weiliu/packages/lemon-1.3/build/wukong/lemon)

# boost library
set(BOOST_LIB "" CACHE FILEPATH "Boost library path")
LINK_DIRECTORIES( ${BOOST_LIB} )

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

ADD_EXECUTABLE(geticc
  geticc.cxx
  )

ADD_EXECUTABLE(geticcs
  geticcs.cxx
  )

ADD_EXECUTABLE(get012
  get012.cxx
  )

ADD_EXECUTABLE(fmriresample
  fmriresample.cxx
  utility.cxx
  )

ADD_EXECUTABLE(gentestlabel
  gentestlabel_base.cxx
  gentestlabel.cxx
  utility.cxx
  )

ADD_EXECUTABLE(gentestlabel2
  gentestlabel_base.cxx
  gentestlabel2.cxx
  utility.cxx
  )

ADD_EXECUTABLE(margll
  margll.cxx
  utility.cxx
  )
ADD_EXECUTABLE(computesnr
  compute_snr.cxx
  utility.cxx
  )


TARGET_LINK_LIBRARIES(groupmrf emon ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(label2ind ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(ind2label ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(sample2label ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(copyheader ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)

TARGET_LINK_LIBRARIES(applymask ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(seedcorr ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)

TARGET_LINK_LIBRARIES(compareclusterings ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)

TARGET_LINK_LIBRARIES(alignlabelmaps ${ITK_LIBRARIES} boost_program_options boost_filesystem boost_system)

TARGET_LINK_LIBRARIES(geticc ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(geticcs ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(get012 ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(fmriresample ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(gentestlabel ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(gentestlabel2 ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(margll ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)
TARGET_LINK_LIBRARIES(computesnr ${ITK_LIBRARIES} boost_program_options boost_filesystem  boost_system)




