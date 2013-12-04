# from partfinder import options, threadpool
# import shlex

#def pytest_generate_tests(metafunc):
#    # This function feeds the output of the above function into the tests below
#    metafunc.parametrize(
#        'good_args', [
#            '-p 3 -v folder',
#            'folder',
#            ]
#    )

# def test_processes():
#     o = options.Options()
#     o.parse_args(shlex.split('folder'))
#     assert o.processes == threadpool.get_cpu_count()
#     assert o.folder_path == 'folder'
#
#     o = options.Options()
#     o.parse_args(shlex.split('-p 3 folder'))
#     assert o.processes == 3
#

