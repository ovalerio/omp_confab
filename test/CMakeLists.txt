# define TESTDATADIR for tests that need input files
add_definitions(-DTESTDATADIR="\\"${CMAKE_SOURCE_DIR}/test/files/\\"")

# define FORMATDIR for location of format plugin binaries
add_definitions(-DFORMATDIR="\\"${openbabel_BINARY_DIR}/lib/\\"")

###########################################################
#  new tests using obtest.h
############################################################

set(tests
    automorphism
    builder
    canonconsistent
    canonstable
    cistrans
    confsearch
    divposes
    dummy
    graphsym
    implicitH
    lssr
    isomorphism
    rotor
    shuffle
    smiles
    spectrophore
    squareplanar
    stereo
    stereoperception
    tetrahedral
    tetranonplanar
    tetraplanar
    uniqueid
  )

if (EIGEN2_FOUND)
set(tests
        align
        ${tests}
)
endif (EIGEN2_FOUND)

set(TEST_PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

foreach(test ${tests})
  add_executable(test_${test} ${test}test.cpp)
  target_link_libraries(test_${test} ${BABEL_LIBRARY})
  if(MSVC)
    add_test(${test}_Test ${CMAKE_BINARY_DIR}/test_${test})
  else(MSVC)
    add_test(${test}_Test ${TEST_PATH}/test_${test})
  endif(MSVC)
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed")
  if(MSVC)
    set_target_properties(test_${test} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
endforeach(test ${tests})


############################################################
#  old-style tests using "not ok"
############################################################

set(old_tests
    aromatest
    atom
    bond
    cansmi
    charge-mmff94
    charge-gasteiger
    cmlreadfile
    conversion
    datatest
    ffgaff
    ffghemical
    ffmmff94
    ffuff
    formalcharge
    format
    formula
    internalcoord
    invalidsmarts
    invalidsmiles
    iterators
    logp_psa
    math
    mol
    phmodel
    residue
    ringtest
    smartstest
    smartsparse
    smilesmatch
    unitcell
   )

# TODO: forgot that inchiwrite requires arguments -- it compiles now, but
#   doesn't run a test
if(WITH_INCHI)
  add_executable(test_inchiwrite inchiwrite.cpp ../src/formats/getinchi.cpp)
  target_link_libraries(test_inchiwrite ${BABEL_LIBRARY})
  if(MSVC)
    set_target_properties(test_inchiwrite PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
  # files in test/inchi -- both .sdf and .txt
  set(inchitests Samples SamplesTechMan)
  set(inchidata ${CMAKE_SOURCE_DIR}/test/inchi)
  foreach(test ${inchitests})
    if(MSVC)
      add_test(inchi${test}_Test 
               ${CMAKE_BINARY_DIR}/test_inchiwrite ${inchidata}/${test}.sdf ${inchidata}/${test}.txt)
    else(MSVC)
      add_test(inchi${test}_Test 
              ${TEST_PATH}/test_inchiwrite ${inchidata}/${test}.sdf ${inchidata}/${test}.txt)
    endif(MSVC)

    set_tests_properties(inchi${test}_Test PROPERTIES
                         FAIL_REGULAR_EXPRESSION "Not ok"
                         ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data"
	 	        )
  endforeach(test ${inchitests})
endif(WITH_INCHI)

foreach(test ${old_tests})
  add_executable(test_${test} ${test}.cpp)
  target_link_libraries(test_${test} ${BABEL_LIBRARY})
  if(MSVC)
    add_test(${test}_Test ${CMAKE_BINARY_DIR}/test_${test})
  else(MSVC)
    add_test(${test}_Test ${TEST_PATH}/test_${test})
  endif(MSVC)
  set_tests_properties(${test}_Test PROPERTIES
    FAIL_REGULAR_EXPRESSION "not ok"
    ENVIRONMENT "BABEL_DATADIR=${CMAKE_SOURCE_DIR}/data"
 )
  if(MSVC)
    set_target_properties(test_${test} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  endif(MSVC)
endforeach(test ${old_tests})

############################################################  
#  roundtrip 
############################################################  

add_executable(roundtrip roundtrip.cpp)
target_link_libraries(roundtrip ${BABEL_LIBRARY})
install(TARGETS roundtrip
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib
)

###########################
# Tests wrapped in Python #
###########################

include(UsePythonTest)
set(pytests
    babel sym fastsearch unique kekule)
foreach(pytest ${pytests})
  ADD_PYTHON_TEST(pytest_${pytest} test${pytest}.py)
  set_tests_properties(pytest_${pytest} PROPERTIES
       FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
       ENVIRONMENT "BABEL_LIBDIR=${CMAKE_BINARY_DIR}/lib"
  )
endforeach(pytest ${pytests})

###############################
# Tests using Python bindings #
###############################

if (PYTHON_BINDINGS)
  set(pybindtests
      bindings _pybel example)
  foreach(pybindtest ${pybindtests})
    ADD_PYTHON_TEST(pybindtest_${pybindtest} test${pybindtest}.py)
    set_tests_properties(pybindtest_${pybindtest} PROPERTIES
         FAIL_REGULAR_EXPRESSION "ERROR;FAIL;Test failed"
         ENVIRONMENT "BABEL_LIBDIR=${CMAKE_BINARY_DIR}/lib"
    )
  endforeach(pybindtest ${pybindtests})
endif (PYTHON_BINDINGS)

