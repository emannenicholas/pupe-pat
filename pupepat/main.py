import matplotlib
matplotlib.use('Agg')

import logging
import logging.config
from lcogt_logging import LCOGTFormatter
import traceback
import sys

logging.captureWarnings(True)

logConf = { "formatters": { "default": {"()": LCOGTFormatter}},
            "handlers": { "console": { "class": "logging.StreamHandler", "formatter": "default",
                                       "stream": "ext://sys.stdout"}},
            "loggers": {"pupepat": { "handlers": ["console"], "level": logging.INFO, 'propagate': False}},
            "version": 1 }

logging.config.dictConfig(logConf)
logger = logging.getLogger('pupepat')

import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from glob import glob
import argparse
import os

from pupepat.fitting import fit_defocused_image
from pupepat.utils import save_results
import tempfile



def run_watcher():
    parser = argparse.ArgumentParser(description='Run the PUPE-PAT analysis on a directory')
    parser.add_argument('--input-dir', dest='input_dir', required=True, help='Input directory where the new files will appear.')
    parser.add_argument('--output-dir', dest='output_dir', required=True, help='Directory to store output files.')

    args = parser.parse_args()
    observer = Observer()
    observer.schedule(Handler(args.output_dir), args.input_dir, recursive=True)
    observer.start()
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
        logger.info('Stopping PUPE-PAT watcher because of keyboard interrupt.')

    observer.join()


def analyze_directory():
    parser = argparse.ArgumentParser(description='Run the PUPE-PAT analysis on a directory')
    parser.add_argument('--input-dir', dest='input_dir', required=True, help='Input directory where the new files will appear.')
    parser.add_argument('--output-dir', dest='output_dir', required=True, help='Directory to store output files.')

    args = parser.parse_args()
    images_to_analyze = glob(os.path.join(args.input_dir, '*x??.fits*'))
    output_table = None
    for image_filename in images_to_analyze:
        output_table = analyze_image(image_filename, output_table, os.path.join(args.output_dir, 'pupe-pat.dat'), args.output_dir)


def analyze_image(filename, output_table, output_filename, output_directory):
    try:
        # Run the analysis on the new frame
        logger.info('Fitting pupil model', extra={'tags': {'filename': os.path.basename(filename)}})
        file_basename, file_extension = os.path.splitext(os.path.basename(filename))

        with tempfile.TemporaryDirectory() as temp_dir:
            if file_extension == '.fz':
                temp_filename = os.path.join(temp_dir, os.path.basename(file_basename))
                os.system('funpack -O {temp_filename} {filename}'.format(temp_filename=temp_filename,
                                                                         filename=filename))
                filename = temp_filename

            plot_basename = os.path.join(output_directory, os.path.splitext(os.path.basename(filename))[0])

            best_fit_models = fit_defocused_image(filename, plot_basename)
            best_fit_models = [best_fit_model for best_fit_model in best_fit_models if best_fit_model is not None]
            if len(best_fit_models) > 0:
                output_table = save_results(filename, best_fit_models, output_table, output_filename)
    except Exception as e:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tb_msg = traceback.format_exception(exc_type, exc_value, exc_tb)

        logger.error('Error processing: {tb_msg}'.format(tb_msg=tb_msg), extra={'tags': {'filename': os.path.basename(os.path.basename(filename)),
                                                                                'error': str(e)}})
    return output_table


class Handler(FileSystemEventHandler):
    def __init__(self, output_directory):
        self.output_directory = output_directory
        self.output_table = None
        self.output_filename = os.path.join(output_directory, 'pupe-pat.dat')
        super(Handler, self).__init__()

    def on_created(self, event):
        if 'x00.fits.fz' in event.src_path:
            self.output_table = analyze_image(event.src_path, self.output_table, self.output_filename, self.output_directory)
