"""
Pytest configuration file to run slow test only if pytest is started with option --runslow or if
environment variable RUN_SLOW_TESTS is set.
"""

import pytest, os

def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False, help="run slow tests")

def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    import platform

    run_slow_tests = config.getoption("--runslow")

    # gitlab-ci: run slow tests only on Linux
    if platform.system()=='Linux' and 'RUN_SLOW_TESTS' in os.environ and os.environ['RUN_SLOW_TESTS']:
        run_slow_tests = True

    if not run_slow_tests:
        skip_slow = pytest.mark.skip(reason="need --runslow option to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)
