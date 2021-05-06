# coding=utf-8

import os.path
import time

import numpy as np


def check_opencv_module():
    try:
        import cv2
    except ImportError:
        return False
    else:
        return True


def sequence2image(sequence, imageXY, colorMap=True):
    if (isinstance(imageXY, tuple) and len(imageXY) == 2 and isinstance(imageXY[0], int) and isinstance(imageXY[1], int)) and (isinstance(sequence, np.ndarray) and sequence.ndim == 1 and sequence.size == imageXY[0] * imageXY[1]):
        sequence.resize(imageXY)
        array2image(sequence, imageXY, colorMap)
    else:
        raise ValueError("Incompatible size of sequence and target image.")


def array2image(array, imageXY, colorMap=True):
    if (check_opencv_module() is True):
        if (isinstance(imageXY, tuple) and len(imageXY) == 2 and isinstance(imageXY[0], int) and isinstance(imageXY[1], int)) and (isinstance(array, np.ndarray) and array.ndim == 2 and array.shape == imageXY):
            imgGrey = cv2.normalize(
                src=array, dst=None, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX, dtype=cv2.CV_8UC1)
            if (colorMap is True):
                imgColor = cv2.applyColorMap(imgGrey, cv2.COLORMAP_JET)
                return imgColor
            else:
                return imgGrey
        else:
            raise ValueError("Incompatible size of array and target image.")
    else:
        raise ImportError(
            "Can not find cv2 module. Please install it first with command \'pip3 install opencv-python\'.")


def writeImage(image, name, dir="."):
    if (check_opencv_module() is True):
        # yyyymmdd-HHMMSS
        timeFormat = time.strftime("%Y%m%d_%H%M%S")
        fileName = name + "_" + timeFormat + ".png"
        cv2.imwrite(os.path.join(dir, fileName), image)
    else:
        raise ImportError(
            "Can not find cv2 module. Please install it first with command \'pip3 install opencv-python\'.")
