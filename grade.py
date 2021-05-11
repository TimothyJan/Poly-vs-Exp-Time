
from saddl import *

a = Assignment(18, 9)
s = State(a)

horizontal_rule()

s.reject_if_missing_contributors()
s.reject_if_starter_contributors()

s.reject_unless_files_exist(['poly_exp.hpp',
                             'poly_exp_test.cpp'])

s.reject_if_file_unchanged('poly_exp.hpp',
                           '953ed73434f4cae54b9161b48da2f25a2622522198a655c00de571bb596b16df')

s.reject_if_file_changed('poly_exp_test.cpp',
                         '126e0d8ff80a4abcd57f6655f6e28c410bdbfe983226f2ca44115a161565a9af')

s.reject_unless_command_succeeds(['make', 'clean', 'poly_exp_test'])

s.string_removed_test('TODO comments removed', 3, 'TODO', ['poly_exp.hpp'])

s.gtest_run('poly_exp_test')
s.gtest_suite_test('max_subarray_exh_SmallCases', 3)
s.gtest_suite_test('max_subarray_exh_LargeCases', 3)
s.gtest_suite_test('max_subarray_dbh_SmallCases', 3)
s.gtest_suite_test('max_subarray_dbh_LargeCases', 3)
s.gtest_suite_test('subset_sum_exh', 3)

s.summarize()
