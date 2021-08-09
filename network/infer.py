import torch
import argparse
import os
from utils.detectionRemoval import test_single_highlight_image

def get_parameters():
    parameters = argparse.ArgumentParser()
    parameters.add_argument('-i', '--image_path', type=str, default=None, help='input image path')
    parameters.add_argument('-r', '--results_path', type=str, default='./output', help='folder for saving result')
    parameters.add_argument('-c', '--checkpoint', type=str, default='jshdr', help='my checkpoint')
    return parameters

if __name__ == "__main__":
    parameters = get_parameters().parse_args()
    print("running ... ...")
    test_single_highlight_image(parameters)
    print("done!")
