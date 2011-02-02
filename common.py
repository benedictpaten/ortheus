
import os

def ortheusRootPath():
    """
    function for finding external location
    """
    import ortheus.Ortheus
    i = os.path.abspath(ortheus.Ortheus.__file__)
    return os.path.split(i)[0] #os.path.split(os.path.split(os.path.split(i)[0])[0])[0]
