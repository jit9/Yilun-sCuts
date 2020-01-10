import glob, os

def combine(ver):
    """Combine splitted outputs into one. The splited files will
    be removed after this.
    Example:
        cuts results combine v0
    """
    tags = ["*.db", "done_list.txt", "error_list.txt"]
    rm_cmds = []
    for tag in tags:
        files = glob.glob("run_{}/{}.*".format(ver, tag))
        for f in files: print(f)
        files = [f for f in files if os.path.basename(f).split('.')[-1]!='db']
        if len(files) == 0:
            print("No files found with tag {}".format(tag))
            continue
        filename = '.'.join(files[0].split('.')[:-1])
        mode = "a" if os.path.exists(filename) else "w"
        if mode == "a":
            first = False
        else:
            first = True
        with open(filename, mode) as ff:
            for f in files:
                with open(f, "r") as tf:
                    lines = tf.readlines()
                    for l in lines:
                        if l[0] == '#':  # headers comment
                            if first:
                                ff.write(l)
                            else:
                                continue
                        else:
                            ff.write(l)
                first = False
        for f in files:
            rm_cmds.append("rm {}".format(f))
    return rm_cmds
