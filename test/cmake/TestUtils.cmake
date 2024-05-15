################################################################################
# Utility functions for creating tests

if(CCPP_ENABLE_MEMCHECK)
  find_program(MEMORYCHECK_COMMAND "valgrind")
endif()

################################################################################
# Runs a test with memory checking if enabled

function(add_memory_check_test test_name test_binary test_args working_dir)
  if(CCPP_ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name}
      COMMAND mpirun -v -np 1 ${MEMORYCHECK_COMMAND} --leak-check=full --error-exitcode=1 --trace-children=yes
      ${test_binary} ${test_args}
      WORKING_DIRECTORY ${working_dir})
  endif()
endfunction(add_memory_check_test)

################################################################################