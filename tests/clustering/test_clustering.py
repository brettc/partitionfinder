import os
import pytest
from partfinder import main, util, analysis, config

# Get the Errors all into the local space
from partfinder.config import ConfigurationError
from partfinder.util import PartitionFinderError

HERE = os.path.abspath(os.path.dirname(__file__))

test_description = """
DNA_hcluster1  | success              |
DNA_hcluster2  | success              |
DNA_hcluster3  | PartitionFinderError |
DNA_hcluster4  | success              | --weights "1,1,1,1"
DNA_hcluster5  | PartitionFinderError | --weights "0,0,0,0"
DNA_hcluster6  | success              | --weights "1000,0.01,0.003,0"
DNA_hcluster7  | PartitionFinderError | --weights "-1000,0.01,0.003,0"
DNA_rcluster1  | success              |
DNA_rcluster2  | success              |

DNA_rcluster3  | success              |

DNA_rcluster4  | success              | --weights "1,1,1,1"
DNA_rcluster5  | PartitionFinderError | --weights "0,0,0,0"
DNA_rcluster6  | success              | --weights "1000, 0.01, 0.003, 0"
DNA_rcluster7  | ConfigurationError   | --weights "-1000, 0.01, 0.003, 0"
DNA_rcluster8  | success              | --rcluster-percent 0
DNA_rcluster9  | success              | --rcluster-percent 0.000001
DNA_rcluster10 | success              | --rcluster-percent 99.999999
DNA_rcluster11 | success              | --rcluster-percent 100.00
DNA_rcluster12 | ConfigurationError   | --rcluster-percent 100.001
DNA_rcluster13 | ConfigurationError   | --rcluster-percent -0.001
DNA_rcluster14 | success              | --rcluster-max 1000 --rcluster-percent 10 -p 1

prot_hcluster1 | success              |
prot_hcluster2 | success              |
prot_hcluster3 | success              |
prot_rcluster1 | ConfigurationError   | --weight  "1, 0, 100, egg"
prot_rcluster2 | success              | --rcluster-percent 12 --weight "1,0,0,0"
prot_rcluster3 | PartitionFinderError | --weight "0,0,0,0"
prot_rcluster4 | PartitionFinderError |
prot_rcluster5 | success              |
"""

# Turn the text above into a container of info about the tests
test_container = {}
for line in test_description.split('\n'):
    line = line.strip()
    if not line:
        continue
    name, res, cmdline = [x.strip() for x in line.split('|')]
    if name[:3] == "DNA":
        kind = "DNA"
    else:
        kind = "protein"
    if res == 'success':
        error = None
    elif res[:5] == 'xfail':
        error = 'xfail', res[5:].strip()
    else:
        # Get the Exception from the local space (must be imported)
        error = locals()[res]

    test_container[name] = error, kind, cmdline


def pytest_generate_tests(metafunc):
    # This function feeds the output of the above function into the tests below
    if 'folder_name' in metafunc.fixturenames:
        # Send it sorted order
        k = test_container.keys()
        k.sort()
        metafunc.parametrize('folder_name', k)

# The actual test function
def test_clustering(folder_name):
    full_path = os.path.join(HERE, folder_name)
    error, kind, cmdline = test_container[folder_name]
    if error is None:
        main.call_main(kind, '--no-ml-tree "%s" --raxml %s' % (full_path, cmdline))
    elif type(error) == type((0, 0)):
        pytest.xfail(error[1])
    else:
        with pytest.raises(error):
            main.call_main(kind, '--no-ml-tree "%s" --raxml %s' % (full_path, cmdline))
