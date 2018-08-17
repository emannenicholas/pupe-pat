"""
pupepat.main - console scripts for running pupe-pat

Author
    Curtis McCully (cmccully@lco.global)

License
    GPL v3.0

February 2018
"""
import matplotlib
matplotlib.use('Agg')

import logging.config
from lcogt_logging import LCOGTFormatter
import traceback
import sys

logging.captureWarnings(True)

logConf = {"formatters": {"default": {"()": LCOGTFormatter}},
           "handlers": {"console": {"class": "logging.StreamHandler", "formatter": "default",
                                    "stream": "ext://sys.stdout"}},
           "loggers": {"pupepat": {"handlers": ["console"], "level": logging.DEBUG, 'propagate': False}},
           "version": 1}

logging.config.dictConfig(logConf)
logger = logging.getLogger('pupepat')

import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from glob import glob
import argparse
import os
import yaml

from pupepat.quiver import make_quiver_plot
from pupepat.fitting import fit_defocused_image
from pupepat.utils import save_results, merge_pdfs, should_process_image, config
import tempfile


def handle_args_and_config():
    '''Parse and return the command line arguments. 
    '''
    parser = argparse.ArgumentParser(description='Run the PUPE-PAT analysis on a directory',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--input-dir', dest='input_dir', required=True,
                        help='Input directory where the new files will appear.')
    parser.add_argument('--output-dir', dest='output_dir', required=True,
                        help='Directory to store output files.')
    parser.add_argument('--output-table', dest='output_table', default='pupe-pat.dat',
                        help='Filename of the table of fit results.')
    parser.add_argument('--output-plot', dest='output_plot', default='pupe-pat',
                        help='Filename of the quiver plot of fit results.')
    parser.add_argument('--proposal-id', dest='proposal_id', default='LCOEngineering',
                        help='Proposal ID used to take the pupil plate images')
    parser.add_argument('--config-file', dest='config_filename', default=None,
                        help='Filename of configuration YAML. See config.yaml.example, for example.')

    args = parser.parse_args()

    # update utils.config dict with custom values if supplied, and
    # supplied-or-not, dump the config to output-dir
    config_dump_filename = 'pupe-pat.yaml.out'
    
    if args.config_filename:
        with open(args.config_filename, "r") as yaml_file:
            custom_config = yaml.load(yaml_file)
        config.update(custom_config)
        config_dump_filename = args.config_filename + '.out'
    
    with open(os.path.join(args.output_dir, config_dump_filename), 'w') as output_file:
        yaml.dump(config, output_file, default_flow_style=False)

    return args


def run_watcher():
    args = handle_args_and_config()
    observer = Observer()
    observer.schedule(Handler(args.output_dir, args.output_table, args.proposal_id), args.input_dir, recursive=True)
    observer.start()
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
        logger.info('Stopping PUPE-PAT watcher because of keyboard interrupt.')

    merge_pdfs(args.output_dir, args.output_table)
    make_quiver_plot(args.output_dir, args.output_table, output_plot=args.output_plot)
    observer.join()


def analyze_directory():
    args = handle_args_and_config()
    images_to_analyze = [image for image in glob(os.path.join(args.input_dir, '*fits*'))
                         if should_process_image(image, args.proposal_id)]

    if len(images_to_analyze) == 0:
        logger.error('No images to analyze in directory: {}'.format(args.input_dir))
        sys.exit(1)

    output_table = None
    for image_filename in images_to_analyze:
        output_table = analyze_image(image_filename, output_table, os.path.join(args.output_dir, args.output_table),
                                     args.output_dir)
    merge_pdfs(args.output_dir, args.output_table)
    make_quiver_plot(args.output_dir, args.output_table, output_plot=args.output_plot)


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

        logger.error('Error processing: {tb_msg}'.format(tb_msg=tb_msg),
                     extra={'tags': {'filename': os.path.basename(os.path.basename(filename)), 'error': str(e)}})
    return output_table


class Handler(FileSystemEventHandler):
    def __init__(self, output_directory, output_filename, proposal_id):
        self.output_directory = output_directory
        self.output_table = None
        self.output_filename = os.path.join(output_directory, output_filename)
        self.proposal_id = proposal_id
        super(Handler, self).__init__()

    def on_created(self, event):
        if should_process_image(event.src_path, self.proposal_id):
            self.output_table = analyze_image(event.src_path, self.output_table,
                                              self.output_filename, self.output_directory)
