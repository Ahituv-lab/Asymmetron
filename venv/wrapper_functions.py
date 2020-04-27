import os


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
    return names

def sanitize(paths, orientation, names):
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
    return paths, orientation_paths, names