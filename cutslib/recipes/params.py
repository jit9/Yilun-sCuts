"""This collection contains tools to perform quick
tweak to cutparams"""

def quickfix(var, val, cpar, flag=None):
    """quick fix for a unique identifier in cutparam, does not always work"""
    import os
    # allow adding an extra flag, specifically what i have in mind is -i
    # for inplace editing
    if not flag:
        cmd = fr'sed -r "s/({var}.\s*:\s*)(\S*)/\1{escape_quote(val)},/g" {cpar}'
    else:
        cmd = fr'sed {flag} -r "s/({var}.\s*:\s*)(\S*)/\1{escape_quote(val)},/g" {cpar}'
    os.system(cmd)


####################
# utility function #
####################
def escape_quote(string):
    if "\"" in string:
        return string.replace("\"", "\\\"")
    elif "'" in string:
        return string.replace("'","\\\"")
    else:
        return string
