"""This scripts define the command line interface of the cutslib"""
import importlib
import sys, os

def main():
    argv = sys.argv[1:]
    recipe = argv[0]
    method = argv[1]
    arguments = argv[2:]
    mod_name = "cutslib.recipes.%s" % recipe
    print("Recipe: %s" % mod_name)
    mod = importlib.import_module(mod_name)
    ret = getattr(mod, method)(*arguments)

    if isinstance(ret, list):
        print("======================")
        print("Run: \n")
        for cmd in ret:
            print(cmd)
            os.system(cmd)
        print("======================")
