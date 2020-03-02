"""moby2 related recipes"""

import os

moby2_dir = "/home/yilung/software/moby2"
moby2_prefix = os.environ.get("MOBY2_PREFIX", None)

# update moby2
def update():
    cmd = f"cd {moby2_dir}; git pull origin master; MOBY2_PREFIX={moby2_prefix} make; MOBY2_PREFIX={moby2_prefix} make install"
    return [cmd]
