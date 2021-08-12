import importlib.util

#################################################

def check_package(package_name):
    '''
    To check if the package is available in the
    present python environment

    Parameter :
    package_name : string

    For example, check_package('matplotlib')
    '''
    spec = importlib.util.find_spec(package_name)
    if spec is None:
        return 0
    else :
        return 1

#################################################

req = ['numpy', 'scipy', 'matplotlib']

for i in req:
    found = check_package(i)
    if (found != 1):
        print("\t\tError: Required package not found.  ",i)
        print("\t\tExiting.  ",i)
        quit()


#################################################
print ("\t\tRequired packages available.")
