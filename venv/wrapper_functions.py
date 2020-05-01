import os

# Catch number of bins >0, not string, not list etc
# Catch distance limits, have to be integers
# Catch plots needs to be False / True only
# Patterns needs to only contain +/-
# Threshold p-value, how do people insert very small p-values? How do we explain them how to do it?
# Expected bias needs to be between 0 and 1

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message

def path_checker(paths):
    """:param paths Input given by the user as a string of paths, comma seperated
       :return a list of all paths given by the user (a list with one element if only one path was given)

       This functions takes the input paths given by the user, verifies that they exist and returns a list of paths
    """
    pathsL = [x.strip() for x in paths.split(',')]
    for path in pathsL:
        if not os.path.exists(path):
            err= "File \"" + path + "\" was not found"
            raise IOError(err)
    return pathsL

def name_splitter(names):
    """
    This function splits the optional argument names into a list of the different names

    :param names: Input given by the user as a string of names, corresponding to the paths in the path argument, comma seperated
    :return: a list of all names given by the user (a list with one element if only one path and name was given
    """
    namesL = [x.strip() for x in names.split(',')]
    return namesL

def sanitize(paths, orientation, names):
    """
    All three parameters are entered by the user in string format. If multiple paths are given, they are in "path1, path2",
    comma-seperated format. If orientation or names is given as an optional argument by the user, there needs to be one file
    for each corresponding path.

    This function takes the three inputs, converts them to lists of paths and verifies that each each orientation or name
    file, corresponds to one path file. If orientation or names was not given, they are set to none

    :param paths: String of paths entered by the user. If multiple, as comma-seperated string.
    :param orientation: String of orientation files entered by the user
    :param names: String of names entered by the user. Each name must correspond to a file entered in --paths
    :return: A tuple consisting of the three input arguments in list form
    """
    paths = path_checker(paths)  # Converts the input to a list of paths. List can include only one element, if one path is given by the user
    # Converts the input to a list of path for the orientation files. If this optional argument was not given, the variable is set to None
    try:
        orientation_paths = path_checker(orientation)
    except AttributeError:
        orientation_paths = None
    # Converts the input to a list of names. If this optional argument was not given, the variable is set to none
    try:
        names = name_splitter(names)
    except AttributeError:
        names = None
    print(paths, orientation_paths, names)
    if names is not None:
        if len(paths) != len(names):
            err = "Please enter the same number of arguments for paths and names"
            raise InputError(err)
    if orientation_paths is not None:
        if len(paths)  != len(orientation_paths):
            err = "Please enter the same number of arguments for paths and orientation_paths"
            raise InputError(err)
    return paths, orientation_paths, names
