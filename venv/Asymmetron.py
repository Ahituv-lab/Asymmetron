import sys
#import functions
import argparse



def fun1():
    return


def fun2():
    return


def fun3():
    return

if __name__ == "__main__":
    """
    # Determine which of the three subroutines to run based on the first argument
    # """
    parser = argparse.ArgumentParser()
    parser.parse_args()
    print(test3)
    # if sys.argv[1] == "FUN1":
    #     # Can modify later to allow for multiple paths
    #     path = sys.argv[2]
    #     min_dist = sys.argv[3]
    #     max_dist = sys.argv[4]
    #
    #     fun1()
    # if sys.argv[1] == "FUN2":
    #     fun2()
    # if sys.argv[1] == 'FUN3':
    #     fun3()
