set(examples
#  obconformersearch_default
#  obconversion_readstring
#  obforcefield_energy
  obenergy_mmff94
  obminimize_mmff94
  obconformer_mmff94
  results_mmff94
  validate_mmff94
  confab_mmff94
)

add_library(StopWatch
    hr_time.cpp
)

foreach(example ${examples})
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} ${BABEL_LIBRARY} StopWatch)
  install(TARGETS ${example}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
  )
endforeach(example)
