function(liknorm_add_test name)
  add_executable(test_${name} src/${name}.c)
  target_compile_features(test_${name} PRIVATE c_std_11)
  target_compile_definitions(
    test_${name} PRIVATE $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>)
  target_link_libraries(test_${name} PRIVATE LIKNORM::liknorm)
  add_test(NAME ${name} COMMAND test_${name} -E environment)

  # Ugly hacky to make the tests work on Windows.
  file(TO_CMAKE_PATH "$ENV{PATH}" PATH)
  list(APPEND PATH $<TARGET_FILE_DIR:LIKNORM::liknorm>)
  string(REPLACE ";" "\\;" PATH "${PATH}")
  set_tests_properties(${name} PROPERTIES ENVIRONMENT "PATH=${PATH}")
endfunction()

liknorm_add_test(empty)
liknorm_add_test(linking)
liknorm_add_test(integration)
liknorm_add_test(specific)
liknorm_add_test(logprod)

file(COPY data/table.csv DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)
