import moby2
from moby2.util.database import TODList
import os

input_list = "./2017_ar6_wide01hn.txt"
error_list = "./error_list.txt"

tod_list = TODList.from_file(input_list)
with open(error_list, "w") as f:
    for todn in tod_list[30:100]:
        ret = os.system("python test_radec_load.py %s" % todn)
        print("ret: %s" % ret)
        if ret != 0:
            f.write("%s\n" % todn)
