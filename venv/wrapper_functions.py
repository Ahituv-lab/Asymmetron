import os
import time
import argparse

# Catch number of bins negative or 0 --> check_positive_int
# Catch distance limits have to be >=0 --> check_positive_int
# Patterns needs to only contain +/- --> check_valid_pattern
# Recognize if one of the input files is using a header for column names
# Threshold p-value (with Bonferoni correction) how do people insert very small p-values? How do we explain them how to do it?
# Expected bias needs to be between 0 and 1 --> check_valid_probability QUESTION does 0 or 1 make sense?
# If user uses Scores it should be float or integer.
# Set patterns to default if not included in user-input
# Small file with inputs. Trunctuate path if too long

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

def output_path(fun_name, extension, *args):
    """
    Creates the path to be used for saving output files. File type needs to be added manually depending on the function

    :param fun_name: Function name to append to the path
    :param args: Additional arguments to append to the path
    :param extension: file extension for the file
    :return: Path under which to save the output. No file type is selected
    """

    time_stamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())  # To add timestamp to output file names
    if not os.path.exists("Asymmetron_output"):
        os.makedirs("Asymmetron_output/")
    if not os.path.exists("Asymmetron_output/"+fun_name):
        os.makedirs("Asymmetron_output/" + fun_name)

    return "Asymmetron_output/" + fun_name + "/" + time_stamp+ "_" + fun_name + "_" + "_".join(args)+ "." + extension



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
    return paths, orientation_paths, names


def check_positive_int(value):
    """
    Checks if the value is a non-zero integer. Raises an argparser error otherwise. Returns the integer

    """
    try:
        my_value = int(value)
        if my_value < 0:
            raise
    except:
        msg = "{} is an invalid value. Please enter a positive integer.".format(value)
        raise argparse.ArgumentTypeError(msg)
    return my_value


def check_valid_pattern(value):
    """
    Checks if the value is a valid pattern consisting of + and -. Raises an argparser error otherwise.
    Returns the pattern.
    """
    for char in value:
        if char not in ("+", "-", ",", " "):
            msg = "{} is an invalid pattern. Please enter a pattern consisting of + or - only.".format(value)
            raise argparse.ArgumentTypeError(msg)
    return value


def check_valid_probability(value):
    """
    Checks if the value is a valid probability between 0 and 1. Raises an argparser error otherwise.
    Returns the pattern.
    """
    try:
        my_value = float(value)
        if my_value < 0 or my_value > 1:
            raise
    except:
        msg = "{} is an invalid value. Please enter a probability between 0 and 1.".format(value)
        raise argparse.ArgumentTypeError(msg)



if __name__ == "__name__":
    print(output_path("consecutive_patterns", "bed", "test", "test2"))





