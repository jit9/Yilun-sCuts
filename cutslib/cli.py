"""This scripts define the command line interface of the cutslib"""
import importlib
import sys, os

def main():
    argv = sys.argv[1:]
    if len(argv) < 2:
        print("Too few arguments!")
        print("Format: cuts recipe method arg1 arg2 ...")
        sys.exit(-1)
    recipe = argv[0]
    method = argv[1]

    has_arguments = False
    if len(argv) > 2:
        has_arguments = True
        arguments = argv[2:]

    mod_name = "cutslib.recipes.%s" % recipe
    mod = importlib.import_module(mod_name)
    if has_arguments:
        ret = getattr(mod, method)(*arguments)
    else:
        ret = getattr(mod, method)()

    if isinstance(ret, list):
        print("======================")
        print("Recipe: %s" % mod_name)
        print("Run: \n")
        for cmd in ret:
            print(cmd)
            os.system(cmd)
        print("======================")
