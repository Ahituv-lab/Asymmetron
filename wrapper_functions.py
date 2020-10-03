import os
import time
import argparse
import configparser
import glob


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
    config = configparser.ConfigParser()
    config.read('config.txt')
    time_stamp = time.strftime("%Y%m%d_%H%M%S", time.localtime()) +"_"  # To add timestamp to output file names
    # Remove the time stamp if chosen by the user

    if config['DEFAULT']['time_stamp'] == "False":
        time_stamp = ""
    if not os.path.exists("Asymmetron_output"):
        os.makedirs("Asymmetron_output/")
    if not os.path.exists("Asymmetron_output/"+fun_name):
        os.makedirs("Asymmetron_output/" + fun_name)

    return "Asymmetron_output/" + fun_name + "/" + time_stamp + fun_name + "_" + "_".join(args)+ "." + extension


def bed_file_validate(paths):
    for path in paths:
        with open(path) as my_bed_file:
            for line in my_bed_file.readlines():
                if line[0] == "-":
                    pass
                elif len(line.split("\t")) < 6:
                    raise InputError("BED file not compatible")
    return



def path_checker(paths):
    """:param paths Input given by the user as a string of paths, comma seperated
       :return a list of all paths given by the user (a list with one element if only one path was given)

       This functions takes the input paths given by the user, verifies that they exist and returns a list of paths
    """
    pathsL = [glob.glob(x.strip()) for x in paths.split(',')]
    pathsL = [path for glob_paths in pathsL for path in glob_paths]  # Convert paths to simple list
    for path in pathsL:
        if not os.path.exists(path):
            err= "File \"" + path + "\" was not found"
            raise InputError(err)
    return pathsL

def name_splitter(names, paths):
    """
    This function splits the optional argument names into a list of the different names

    :param names: Input given by the user as a string of names, corresponding to the paths in the path argument, comma seperated
    :param paths: List of paths used to extract the names, if the names argument is None
    :return: a list of all names given by the user (a list with one element if only one path and name was given
    """
    if names != None:
        namesL = [x.strip() for x in names.split(',')]
        if len(namesL) != len(paths):
            err = "If the names argument is used, a name for every corresponding path must be provided"
            raise InputError(err)
    else:
        namesL = [os.path.basename(x) for x in paths]
    return namesL

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
    chars = set(value)
    if not (chars.issubset(("+", "-", ",", " ", ".")) or value == "basic"):
        msg = "{} is an invalid pattern. Please enter a pattern consisting of + or - only. Insert \"basic\" for " \
              "same / opposite orientations " \
              "".format(value)
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


if __name__ == "__main__":
    print(output_path("consecutive_patterns", "bed", "test", "test2"))





