"""Parse the environment variables for the use of the entire
package"""

import os

CUTS_USER = os.environ.get("USER", "yilung")
CUTS_DIR = os.environ.get("CUTS_DIR","/projects/ACT/yilung/cuts/")
CUTS_PYENV = os.environ.get("CUTS_PYENV", "myenv")
CUTS_DEPOT = os.environ.get("CUTS_DEPOT", "/projects/ACT/yilung/depot/")
