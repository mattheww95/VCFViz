"""
Createing a simple logging class to provide better and cleaner information
on the vcf creation

2022-05-19: Matthew Wells
"""
import logging

class VCFLogger(object):
    logger = logging.getLogger()
    FORMAT = "[%(levelname)s %(asctime)s] %(message)s"
    logger.setLevel(logging.DEBUG)
    logging.basicConfig(format=FORMAT)


