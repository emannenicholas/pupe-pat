import matplotlib
matplotlib.use('Agg')

import logging
import logging.config
from lcogt_logging import LCOGTFormatter


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


class Handler(FileSystemEventHandler):
    def __init__(self, output_directory):
        self.output_directory = output_directory
        self.output_table = None
        self.output_filename = os.path.join(output_directory, 'pupe-pat.dat')
        super(Handler, self).__init__()

    def on_created(self, event):
        if 'x00.fits' in event.src_path:
            try:
                # Run the analysis on the new frame
                logger.info('Fitting pupil model', extra={'tags': {'filename': os.path.basename(event.src_path)}})
                file_basename, file_extension = os.path.splitext(os.path.basename(event.src_path))

                with tempfile.TemporaryDirectory() as temp_dir:
                    if file_extension == '.fz':
                        temp_filename = os.path.join(temp_dir, os.path.basename(file_basename))
                        os.system('funpack -O {temp_filename} {filename}'.format(temp_filename=temp_filename, filename=event.src_path))
                        filename = temp_filename
                    else:
                        filename = event.src_path
                        
                    plot_basename = os.path.join(self.output_directory, os.path.splitext(os.path.basename(filename))[0])

                    best_fit_models = fit_defocused_image(filename, plot_basename)
                    self.output_table = save_results(os.path.basename(filename), best_fit_models, self.output_table,
                                                     self.output_filename)
            except Exception as e:
                logger.error('Error processing', extra={'tags': {'filename': os.path.basename(event.src_path),
                                                                 'error': e}})
